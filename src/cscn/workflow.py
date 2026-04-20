from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from .aggregate import write_consensus_csv
from .config import CSCNConfig, serialize_config
from .io import load_dataset
from .layout import RunLayout
from .prep import PreparedGroupData, PreparedRunData, prepare_run_inputs


@dataclass(frozen=True)
class RunSummary:
    run_dir: Path
    groups: dict[str, int]
    gene_count: int


def _log(message: str) -> None:
    print(f"[cscn] {message}")


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    import csv

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_snapshot(layout: RunLayout, config: CSCNConfig) -> None:
    layout.run_dir.mkdir(parents=True, exist_ok=True)
    layout.config_snapshot_path.write_text(
        yaml.safe_dump(serialize_config(config), sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )


def _write_prepare_outputs(layout: RunLayout, prepared: PreparedRunData) -> RunSummary:
    layout.ensure_dirs()
    gene_rows = [
        {"node_index": index, "gene_name": gene_name}
        for index, gene_name in enumerate(prepared.gene_names)
    ]
    _write_csv(layout.genes_path, ["node_index", "gene_name"], gene_rows)

    group_rows = []
    for group_key, group_data in prepared.groups.items():
        np.save(layout.matrix_path(group_key), group_data.matrix.astype(np.float32))
        _write_csv(
            layout.group_cells_path(group_key),
            ["cell_id"],
            [{"cell_id": cell_id} for cell_id in group_data.cell_ids],
        )
        group_rows.append(
            {
                "group_key": group_key,
                "group_label": group_key,
                "n_cells": len(group_data.cell_ids),
            }
        )
    _write_csv(layout.group_manifest_path, ["group_key", "group_label", "n_cells"], group_rows)
    prepared.sampled_metadata.to_csv(layout.cell_metadata_path, index=False)

    summary = {
        "run_name": layout.run_name,
        "gene_count": len(prepared.gene_names),
        "groups": {group: len(data.cell_ids) for group, data in prepared.groups.items()},
    }
    layout.summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return RunSummary(
        run_dir=layout.run_dir,
        groups=summary["groups"],
        gene_count=summary["gene_count"],
    )


def prepare_run(config: CSCNConfig) -> RunSummary:
    layout = RunLayout.from_config(config)
    _log(f"loading dataset for run {layout.run_name}")
    dataset = load_dataset(config)
    prepared = prepare_run_inputs(dataset, config)
    _write_snapshot(layout, config)
    summary = _write_prepare_outputs(layout, prepared)
    _log(f"prepared run directory: {summary.run_dir}")
    return summary


def load_prepared_run(layout: RunLayout) -> PreparedRunData:
    if not layout.genes_path.is_file():
        raise FileNotFoundError(
            f"Prepared inputs are missing for run {layout.run_name}. Run `cscn prepare` first."
        )
    genes_frame = pd.read_csv(layout.genes_path)
    gene_names = genes_frame["gene_name"].astype(str).tolist()

    groups_frame = pd.read_csv(layout.group_manifest_path)
    groups: dict[str, PreparedGroupData] = {}
    for row in groups_frame.to_dict(orient="records"):
        group_key = str(row["group_key"])
        matrix = np.load(layout.matrix_path(group_key)).astype(np.float32)
        cells_frame = pd.read_csv(layout.group_cells_path(group_key))
        cell_ids = cells_frame["cell_id"].astype(str).tolist()
        groups[group_key] = PreparedGroupData(
            group_key=group_key,
            cell_ids=cell_ids,
            matrix=matrix,
        )
    sampled_metadata = pd.read_csv(layout.cell_metadata_path)
    return PreparedRunData(
        gene_names=gene_names,
        groups=groups,
        sampled_metadata=sampled_metadata,
    )


def run_cscn(config: CSCNConfig) -> RunSummary:
    from .core import CSCN

    layout = RunLayout.from_config(config)
    prepared = load_prepared_run(layout)
    max_workers = config.run.max_workers or min(8, os.cpu_count() or 1)

    for group_key, group_data in prepared.groups.items():
        dag_dir = layout.dag_group_dir(group_key)
        dag_dir.mkdir(parents=True, exist_ok=True)
        cscn = CSCN(
            output_dir=str(dag_dir),
            sigmoid_score=config.run.sigmoid_score,
            significance_level=config.run.significance_level,
            max_cond_vars=config.run.max_cond_vars,
            use_bitmap=config.run.use_bitmap,
            show_progress=config.run.show_progress,
            progress_interval=config.run.progress_interval,
        )
        _log(
            f"running CSCN for group={group_key} cells={len(group_data.cell_ids)} genes={group_data.matrix.shape[1]}"
        )
        cscn.run_core(group_data.matrix, usingNMF=config.run.using_nmf)
        cscn.run_pc_concurrently(
            max_workers=max_workers,
            progress_interval=config.run.progress_interval,
            progress_label=group_key,
        )
        CSCN.save_to_file(cscn, layout.cscn_object_path(group_key))

    return RunSummary(
        run_dir=layout.run_dir,
        groups={group: len(data.cell_ids) for group, data in prepared.groups.items()},
        gene_count=len(prepared.gene_names),
    )


def aggregate_run(config: CSCNConfig) -> RunSummary:
    layout = RunLayout.from_config(config)
    prepared = load_prepared_run(layout)
    summary_rows = []
    if config.aggregate.consensus:
        from .aggregate import load_node_map

        node_map = load_node_map(layout.run_dir)
        for group_key in prepared.groups:
            metadata = write_consensus_csv(
                group_dir=layout.dag_group_dir(group_key),
                output_path=layout.consensus_csv_path(group_key),
                node_map=node_map,
                threshold_mode=config.aggregate.consensus_threshold_mode,
            )
            summary_rows.append({"group_key": group_key, **metadata})
    summary_path = layout.consensus_dir / "summary.json"
    summary_path.write_text(
        json.dumps(summary_rows, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return RunSummary(
        run_dir=layout.run_dir,
        groups={group: len(data.cell_ids) for group, data in prepared.groups.items()},
        gene_count=len(prepared.gene_names),
    )


def run_biomarker_workflow(config: CSCNConfig) -> Path:
    from .postprocess.biomarker import run_biomarker

    layout = RunLayout.from_config(config)
    prepared = load_prepared_run(layout)
    return run_biomarker(config=config, layout=layout, prepared_run=prepared)


def run_all(config: CSCNConfig) -> RunSummary:
    summary = prepare_run(config)
    run_cscn(config)
    aggregate_run(config)
    if config.biomarker.enabled:
        run_biomarker_workflow(config)
    return summary
