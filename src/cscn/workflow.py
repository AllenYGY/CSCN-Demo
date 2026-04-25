from __future__ import annotations

import json
import math
import os
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from .aggregate import load_group_dags, write_consensus_csv
from .config import CSCNConfig, serialize_config
from .io import load_dataset
from .layout import RunLayout, safe_component
from .prep import (
    PreparedAdaptiveBlockData,
    PreparedGroupData,
    PreparedRunData,
    prepare_run_inputs,
)
from .spatial_blocks import CellBlockAssignment, SpatialBlock


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


def _summary_payload(layout: RunLayout, prepared: PreparedRunData) -> dict[str, Any]:
    return {
        "run_name": layout.run_name,
        "gene_count": len(prepared.gene_names),
        "groups": {group: len(data.cell_ids) for group, data in prepared.groups.items()},
        "adaptive_blocks": {
            group: (
                len(data.adaptive_blocks.blocks)
                if data.adaptive_blocks is not None
                else 0
            )
            for group, data in prepared.groups.items()
        },
    }


def _write_summary(layout: RunLayout, payload: dict[str, Any]) -> None:
    layout.summary_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def _read_summary(layout: RunLayout) -> dict[str, Any]:
    if not layout.summary_path.is_file():
        return {}
    return json.loads(layout.summary_path.read_text(encoding="utf-8"))


def _write_adaptive_block_outputs(layout: RunLayout, group_key: str, adaptive: PreparedAdaptiveBlockData) -> None:
    manifest_rows = []
    block_cell_rows = []
    for block_index, block in enumerate(adaptive.blocks):
        manifest_rows.append(
            {
                "block_index": block_index,
                "block_id": block.block_id,
                "cluster_id": block.cluster_id,
                "center_x": block.center_x,
                "center_y": block.center_y,
                "min_x": block.min_x,
                "max_x": block.max_x,
                "min_y": block.min_y,
                "max_y": block.max_y,
                "n_cells": len(block.cell_ids),
            }
        )
        for cell_id in block.cell_ids:
            block_cell_rows.append({"block_id": block.block_id, "cell_id": cell_id})
    _write_csv(
        layout.block_manifest_path(group_key),
        [
            "block_index",
            "block_id",
            "cluster_id",
            "center_x",
            "center_y",
            "min_x",
            "max_x",
            "min_y",
            "max_y",
            "n_cells",
        ],
        manifest_rows,
    )
    _write_csv(
        layout.block_cells_path(group_key),
        ["block_id", "cell_id"],
        block_cell_rows,
    )
    assignment_rows = []
    for assignment in adaptive.assignments.values():
        assignment_rows.append(
            {
                "cell_id": assignment.cell_id,
                "primary_block_id": assignment.primary_block_id or "",
                "covering_block_ids": "|".join(assignment.covering_block_ids),
                "halo_block_ids": "|".join(assignment.halo_block_ids),
            }
        )
    _write_csv(
        layout.cell_block_assignment_path(group_key),
        ["cell_id", "primary_block_id", "covering_block_ids", "halo_block_ids"],
        assignment_rows,
    )
    _write_csv(
        layout.block_genes_path(group_key),
        ["node_index", "gene_name"],
        [
            {"node_index": index, "gene_name": gene_name}
            for index, gene_name in enumerate(adaptive.block_gene_names)
        ],
    )
    np.save(layout.block_matrix_path(group_key, "sum_then_normalize"), adaptive.matrix_sum.astype(np.float32))
    np.save(layout.block_matrix_path(group_key, "mean"), adaptive.matrix_mean.astype(np.float32))


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
        np.save(layout.raw_matrix_path(group_key), group_data.raw_matrix.astype(np.float32))
        _write_csv(
            layout.group_cells_path(group_key),
            ["cell_id"],
            [{"cell_id": cell_id} for cell_id in group_data.cell_ids],
        )
        if group_data.spatial_coords is not None:
            coord_columns = ["spatial_x", "spatial_y"]
            if group_data.spatial_coords.shape[1] == 3:
                coord_columns.append("spatial_z")
            spatial_rows = []
            for idx, cell_id in enumerate(group_data.cell_ids):
                row = {"cell_id": cell_id}
                for col_idx, col_name in enumerate(coord_columns):
                    row[col_name] = float(group_data.spatial_coords[idx, col_idx])
                spatial_rows.append(row)
            _write_csv(
                layout.spatial_coords_path(group_key),
                ["cell_id", *coord_columns],
                spatial_rows,
            )
        if group_data.adaptive_blocks is not None:
            _write_adaptive_block_outputs(layout, group_key, group_data.adaptive_blocks)
        group_rows.append(
            {
                "group_key": group_key,
                "group_label": group_key,
                "n_cells": len(group_data.cell_ids),
            }
        )
    _write_csv(layout.group_manifest_path, ["group_key", "group_label", "n_cells"], group_rows)
    prepared.sampled_metadata.to_csv(layout.cell_metadata_path, index=False)

    summary = _summary_payload(layout, prepared)
    _write_summary(layout, summary)
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


def _load_adaptive_blocks(layout: RunLayout, group_key: str) -> PreparedAdaptiveBlockData | None:
    manifest_path = layout.block_manifest_path(group_key)
    if not manifest_path.is_file():
        return None

    manifest_frame = pd.read_csv(manifest_path)
    block_cells_frame = pd.read_csv(layout.block_cells_path(group_key))
    block_genes_frame = pd.read_csv(layout.block_genes_path(group_key))
    assignment_frame = pd.read_csv(layout.cell_block_assignment_path(group_key))
    block_to_cells: dict[str, list[str]] = {}
    for row in block_cells_frame.to_dict(orient="records"):
        block_to_cells.setdefault(str(row["block_id"]), []).append(str(row["cell_id"]))

    blocks: list[SpatialBlock] = []
    for row in manifest_frame.to_dict(orient="records"):
        block_id = str(row["block_id"])
        blocks.append(
            SpatialBlock(
                block_id=block_id,
                cluster_id=int(row["cluster_id"]),
                center_x=float(row["center_x"]),
                center_y=float(row["center_y"]),
                min_x=float(row["min_x"]),
                max_x=float(row["max_x"]),
                min_y=float(row["min_y"]),
                max_y=float(row["max_y"]),
                cell_ids=tuple(block_to_cells.get(block_id, [])),
            )
        )
    assignments: dict[str, CellBlockAssignment] = {}
    for row in assignment_frame.to_dict(orient="records"):
        assignments[str(row["cell_id"])] = CellBlockAssignment(
            cell_id=str(row["cell_id"]),
            primary_block_id=(
                str(row["primary_block_id"]).strip()
                if str(row.get("primary_block_id") or "").strip()
                else None
            ),
            covering_block_ids=tuple(
                item for item in str(row.get("covering_block_ids") or "").split("|") if item
            ),
            halo_block_ids=tuple(
                item for item in str(row.get("halo_block_ids") or "").split("|") if item
            ),
        )
    block_gene_names = block_genes_frame["gene_name"].astype(str).tolist()
    return PreparedAdaptiveBlockData(
        blocks=blocks,
        assignments=assignments,
        block_gene_names=block_gene_names,
        matrix_sum=np.load(layout.block_matrix_path(group_key, "sum_then_normalize")).astype(np.float32),
        matrix_mean=np.load(layout.block_matrix_path(group_key, "mean")).astype(np.float32),
    )


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
        raw_matrix_path = layout.raw_matrix_path(group_key)
        raw_matrix = (
            np.load(raw_matrix_path).astype(np.float32)
            if raw_matrix_path.is_file()
            else matrix.astype(np.float32)
        )
        cells_frame = pd.read_csv(layout.group_cells_path(group_key))
        cell_ids = cells_frame["cell_id"].astype(str).tolist()
        spatial_coords = None
        spatial_path = layout.spatial_coords_path(group_key)
        if spatial_path.is_file():
            spatial_frame = pd.read_csv(spatial_path)
            expected = [
                "cell_id",
                *[
                    f"spatial_{axis}"
                    for axis in ("x", "y", "z")
                    if f"spatial_{axis}" in spatial_frame.columns
                ],
            ]
            spatial_frame = spatial_frame[expected]
            spatial_frame["cell_id"] = spatial_frame["cell_id"].astype(str)
            if spatial_frame["cell_id"].tolist() != cell_ids:
                raise ValueError(
                    f"Spatial coordinate order does not match cell order for group {group_key}."
                )
            spatial_coords = spatial_frame.drop(columns=["cell_id"]).to_numpy(
                dtype=np.float32,
                copy=True,
            )
        groups[group_key] = PreparedGroupData(
            group_key=group_key,
            cell_ids=cell_ids,
            matrix=matrix,
            raw_matrix=raw_matrix,
            spatial_coords=spatial_coords,
            adaptive_blocks=_load_adaptive_blocks(layout, group_key),
        )
    sampled_metadata = pd.read_csv(layout.cell_metadata_path)
    return PreparedRunData(
        gene_names=gene_names,
        groups=groups,
        sampled_metadata=sampled_metadata,
    )


def _resolve_prior_threshold(mode: str | int, num_sources: int) -> int:
    if num_sources <= 0:
        return 1
    if isinstance(mode, int):
        return max(1, mode)
    if str(mode).lower() == "auto":
        return max(1, int(math.ceil(num_sources / 2.0)))
    return max(1, int(mode))


def _write_block_prior_files(
    layout: RunLayout,
    group_key: str,
    adaptive: PreparedAdaptiveBlockData,
) -> dict[str, set[tuple[str, str]]]:
    prior_dir = layout.block_prior_group_dir(group_key)
    prior_dir.mkdir(parents=True, exist_ok=True)
    block_priors: dict[str, set[tuple[str, str]]] = {}
    gene_names = adaptive.block_gene_names
    node_map = {index: gene for index, gene in enumerate(gene_names)}
    for dag_id, dag in load_group_dags(layout.block_dag_group_dir(group_key)):
        if dag_id < 0 or dag_id >= len(adaptive.blocks):
            continue
        block_id = adaptive.blocks[dag_id].block_id
        edges: set[tuple[str, str]] = set()
        for source, target in dag.edges():
            source_gene = node_map.get(int(source))
            target_gene = node_map.get(int(target))
            if not source_gene or not target_gene or source_gene == target_gene:
                continue
            edges.add(tuple(sorted((source_gene, target_gene))))
        block_priors[block_id] = edges
        _write_csv(
            prior_dir / f"{safe_component(block_id)}.csv",
            ["gene_a", "gene_b"],
            [{"gene_a": source, "gene_b": target} for source, target in sorted(edges)],
        )
    return block_priors


def _run_block_level_cscn(
    config: CSCNConfig,
    layout: RunLayout,
    group_key: str,
    adaptive: PreparedAdaptiveBlockData,
    max_workers: int,
):
    from .core import CSCN

    if not adaptive.blocks:
        return {}
    dag_dir = layout.block_dag_group_dir(group_key)
    dag_dir.mkdir(parents=True, exist_ok=True)
    aggregation = config.run.spatial.default_aggregation
    matrix = (
        adaptive.matrix_sum
        if aggregation == "sum_then_normalize"
        else adaptive.matrix_mean
    )
    block_cscn = CSCN(
        output_dir=str(dag_dir),
        sigmoid_score=config.run.sigmoid_score,
        significance_level=config.run.significance_level,
        max_cond_vars=config.run.max_cond_vars,
        use_bitmap=config.run.use_bitmap,
        show_progress=config.run.show_progress,
        progress_interval=config.run.progress_interval,
        spatial_enabled=False,
        spatial_strategy="weighted_counts",
        spatial_mode=config.run.spatial.mode,
        spatial_k=config.run.spatial.k,
        spatial_radius=config.run.spatial.radius,
        spatial_kernel=config.run.spatial.kernel,
        spatial_bandwidth=config.run.spatial.bandwidth,
        spatial_lambda_expr=config.run.spatial.lambda_expr,
        spatial_min_effective_neighbors=config.run.spatial.min_effective_neighbors,
    )
    block_cscn.run_core(matrix, usingNMF=False)
    block_cscn.run_pc_concurrently(
        max_workers=max_workers,
        progress_interval=config.run.progress_interval,
        progress_label=f"{group_key}-blocks",
    )
    return _write_block_prior_files(layout, group_key, adaptive)


def _build_adaptive_group_runtime(
    config: CSCNConfig,
    group_data: PreparedGroupData,
    prepared: PreparedRunData,
    block_priors: dict[str, set[tuple[str, str]]],
) -> tuple[dict[int, np.ndarray], dict[int, set[tuple[int, int]]], dict[str, Any]]:
    assert group_data.adaptive_blocks is not None
    adaptive = group_data.adaptive_blocks
    gene_index = {gene: idx for idx, gene in enumerate(prepared.gene_names)}
    cell_index = {cell_id: idx for idx, cell_id in enumerate(group_data.cell_ids)}
    block_cell_index_sets = {
        block.block_id: {cell_index[cell_id] for cell_id in block.cell_ids if cell_id in cell_index}
        for block in adaptive.blocks
    }
    local_subsets: dict[int, np.ndarray] = {}
    allowed_edges_by_cell: dict[int, set[tuple[int, int]]] = {}
    prior_empty_cells = 0
    fallback_cells = 0
    local_sizes: list[int] = []
    prior_sizes: list[int] = []

    for cell_id, assignment in adaptive.assignments.items():
        if cell_id not in cell_index:
            continue
        idx = cell_index[cell_id]
        subset_blocks = [
            block_id
            for block_id in [assignment.primary_block_id, *assignment.halo_block_ids]
            if block_id
        ]
        subset_indices = sorted(
            {
                item
                for block_id in subset_blocks
                for item in block_cell_index_sets.get(block_id, set())
            }
        )
        if not subset_indices:
            subset_indices = [idx]
        local_subsets[idx] = np.asarray(subset_indices, dtype=np.int64)
        local_sizes.append(len(subset_indices))

        source_blocks = list(assignment.covering_block_ids)
        if not source_blocks and assignment.primary_block_id:
            source_blocks = [assignment.primary_block_id]
        edge_counter: Counter[tuple[str, str]] = Counter()
        for block_id in source_blocks:
            for edge in block_priors.get(block_id, set()):
                edge_counter[edge] += 1
        threshold = _resolve_prior_threshold(
            config.run.spatial.prior_consensus_threshold,
            len(source_blocks),
        )
        allowed_edges = {
            tuple(sorted((gene_index[source], gene_index[target])))
            for (source, target), count in edge_counter.items()
            if count >= threshold and source in gene_index and target in gene_index
        }
        if not allowed_edges:
            prior_empty_cells += 1
            if config.run.spatial.fallback_to_local_knn_subset:
                fallback_cells += 1
        else:
            allowed_edges_by_cell[idx] = allowed_edges
        prior_sizes.append(len(allowed_edges))

    metrics = {
        "block_count": len(adaptive.blocks),
        "valid_block_count": len(block_priors),
        "fallback_to_local_knn_subset_cells": fallback_cells,
        "empty_prior_cells": prior_empty_cells,
        "avg_local_subset_size": float(np.mean(local_sizes)) if local_sizes else 0.0,
        "avg_prior_edge_count": float(np.mean(prior_sizes)) if prior_sizes else 0.0,
    }
    return local_subsets, allowed_edges_by_cell, metrics


def run_cscn(config: CSCNConfig) -> RunSummary:
    from .core import CSCN

    layout = RunLayout.from_config(config)
    prepared = load_prepared_run(layout)
    max_workers = config.run.max_workers or min(8, os.cpu_count() or 1)
    summary_payload = _read_summary(layout)
    adaptive_summary: dict[str, Any] = {}

    for group_key, group_data in prepared.groups.items():
        dag_dir = layout.dag_group_dir(group_key)
        dag_dir.mkdir(parents=True, exist_ok=True)
        adaptive_local_subsets: dict[int, np.ndarray] | None = None
        adaptive_allowed_edges: dict[int, set[tuple[int, int]]] | None = None

        if (
            config.run.spatial.enabled
            and config.run.spatial.strategy == "adaptive_block_prior"
            and group_data.adaptive_blocks is not None
        ):
            _log(
                f"running block-level CSCN for group={group_key} blocks="
                f"{len(group_data.adaptive_blocks.blocks)}"
            )
            block_priors = _run_block_level_cscn(
                config=config,
                layout=layout,
                group_key=group_key,
                adaptive=group_data.adaptive_blocks,
                max_workers=max_workers,
            )
            adaptive_local_subsets, adaptive_allowed_edges, metrics = _build_adaptive_group_runtime(
                config=config,
                group_data=group_data,
                prepared=prepared,
                block_priors=block_priors,
            )
            adaptive_summary[group_key] = metrics

        cscn = CSCN(
            output_dir=str(dag_dir),
            sigmoid_score=config.run.sigmoid_score,
            significance_level=config.run.significance_level,
            max_cond_vars=config.run.max_cond_vars,
            use_bitmap=config.run.use_bitmap,
            show_progress=config.run.show_progress,
            progress_interval=config.run.progress_interval,
            spatial_enabled=config.run.spatial.enabled,
            spatial_strategy=config.run.spatial.strategy,
            spatial_mode=config.run.spatial.mode,
            spatial_k=config.run.spatial.k,
            spatial_radius=config.run.spatial.radius,
            spatial_kernel=config.run.spatial.kernel,
            spatial_bandwidth=config.run.spatial.bandwidth,
            spatial_lambda_expr=config.run.spatial.lambda_expr,
            spatial_min_effective_neighbors=config.run.spatial.min_effective_neighbors,
            adaptive_local_subsets=adaptive_local_subsets,
            adaptive_allowed_edges=adaptive_allowed_edges,
            adaptive_fallback_to_local_knn_subset=config.run.spatial.fallback_to_local_knn_subset,
        )
        _log(
            f"running CSCN for group={group_key} cells={len(group_data.cell_ids)} "
            f"genes={group_data.matrix.shape[1]}"
        )
        cscn.run_core(
            group_data.matrix,
            usingNMF=config.run.using_nmf,
            spatial_coords=group_data.spatial_coords,
        )
        cscn.run_pc_concurrently(
            max_workers=max_workers,
            progress_interval=config.run.progress_interval,
            progress_label=group_key,
        )
        CSCN.save_to_file(cscn, layout.cscn_object_path(group_key))

    if adaptive_summary:
        summary_payload["adaptive_block_prior_runtime"] = adaptive_summary
        _write_summary(layout, summary_payload)

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
