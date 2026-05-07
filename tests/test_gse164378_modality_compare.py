from __future__ import annotations

import csv
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from scripts.analysis.compare_gse164378_modalities import run_analysis
from scripts.prep.prepare_GSE164378 import sample_shared_cells
from tests.support.fake_graph import FakeGraph


def _write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def _make_run_dir(base: Path, name: str, genes: list[str], ckm: np.ndarray, missing_last_dag: bool = False) -> Path:
    run_dir = base / name
    (run_dir / "inputs").mkdir(parents=True, exist_ok=True)
    (run_dir / "matrices").mkdir(parents=True, exist_ok=True)
    (run_dir / "ckm").mkdir(parents=True, exist_ok=True)
    (run_dir / "dags" / "all").mkdir(parents=True, exist_ok=True)

    _write_csv(
        run_dir / "inputs" / "genes.csv",
        ["node_index", "gene_name"],
        [[idx, gene] for idx, gene in enumerate(genes)],
    )
    _write_csv(
        run_dir / "inputs" / "groups.csv",
        ["group_key", "group_label", "n_cells"],
        [["all", "all", 4]],
    )
    _write_csv(
        run_dir / "inputs" / "cell_metadata.csv",
        ["cell_id", "cscn_group", "celltype.l1"],
        [
            ["c1", "all", "A"],
            ["c2", "all", "A"],
            ["c3", "all", "B"],
            ["c4", "all", "B"],
        ],
    )
    _write_csv(
        run_dir / "matrices" / "all_cells.csv",
        ["cell_id"],
        [["c1"], ["c2"], ["c3"], ["c4"]],
    )
    np.save(run_dir / "ckm" / "all_ckm.npy", ckm)

    last = 3 if missing_last_dag else 4
    for idx in range(last):
        with (run_dir / "dags" / "all" / f"result_{idx}.pkl").open("wb") as handle:
            pickle.dump(FakeGraph(nodes=[0, 1], edges=[(0, 1)]), handle)
    return run_dir


def test_compare_gse164378_modalities_drops_missing_joint_dag_and_writes_metrics(tmp_path):
    expr_path = tmp_path / "rna_expr.csv"
    adt_expr_path = tmp_path / "adt_expr.csv"
    joint_expr_path = tmp_path / "joint_expr.csv"
    _write_csv(
        expr_path,
        ["cell_id", "G1", "G2"],
        [
            ["c1", 0.0, 0.1],
            ["c2", 0.2, 0.0],
            ["c3", 5.0, 5.1],
            ["c4", 5.2, 5.0],
        ],
    )
    _write_csv(
        adt_expr_path,
        ["cell_id", "ADT_A", "ADT_B"],
        [
            ["c1", 0.1, 0.0],
            ["c2", 0.2, 0.1],
            ["c3", 4.9, 5.0],
            ["c4", 5.1, 4.9],
        ],
    )
    _write_csv(
        joint_expr_path,
        ["cell_id", "G1", "ADT_A"],
        [
            ["c1", 0.0, 0.1],
            ["c2", 0.1, 0.2],
            ["c3", 5.0, 5.0],
            ["c4", 5.2, 5.1],
        ],
    )

    rna_run = _make_run_dir(
        tmp_path,
        "rna",
        ["G1", "G2"],
        np.asarray([[0.0, 0.1], [0.2, 0.0], [5.0, 5.1], [5.2, 5.0]], dtype=float),
    )
    adt_run = _make_run_dir(
        tmp_path,
        "adt",
        ["ADT_A", "ADT_B"],
        np.asarray([[0.1, 0.0], [0.2, 0.1], [4.9, 5.0], [5.1, 4.9]], dtype=float),
    )
    joint_run = _make_run_dir(
        tmp_path,
        "joint",
        ["G1", "ADT_A"],
        np.asarray([[0.0, 0.1], [0.1, 0.2], [5.0, 5.0], [5.2, 5.1]], dtype=float),
        missing_last_dag=True,
    )

    output_dir = tmp_path / "analysis"
    run_analysis(
        rna_expr_path=expr_path,
        adt_expr_path=adt_expr_path,
        joint_expr_path=joint_expr_path,
        rna_run_dir=rna_run,
        adt_run_dir=adt_run,
        joint_run_dir=joint_run,
        output_dir=output_dir,
        label_column="celltype.l1",
        enable_umap=False,
    )

    metrics = pd.read_csv(output_dir / "clustering_metrics.csv")
    shared = pd.read_csv(output_dir / "shared_cells.csv")
    assignments = pd.read_csv(output_dir / "cell_assignments.csv")

    assert metrics["representation"].tolist() == [
        "expr_rna",
        "expr_adt",
        "expr_joint",
        "ckm_rna",
        "ckm_adt",
        "ckm_joint",
    ]
    assert metrics["n_cells"].tolist() == [3, 3, 3, 3, 3, 3]
    assert shared["cell_id"].tolist() == ["c1", "c2", "c3"]
    assert assignments["cell_id"].tolist() == ["c1", "c2", "c3"]


def test_sample_shared_cells_stratifies_evenly_and_is_deterministic():
    header = ["cell_id", "celltype.l1", "dummy"]
    rows: list[list[str]] = []
    for label in ["A", "B", "C", "D"]:
        for idx in range(10):
            rows.append([f"{label}_{idx}", label, "x"])

    sampled_ids, sampled_rows, counts = sample_shared_cells(
        header,
        rows,
        total_cells=8,
        stratify_key="celltype.l1",
        random_seed=7,
    )

    assert len(sampled_ids) == 8
    assert len(sampled_rows) == 8
    assert counts == {"A": 2, "B": 2, "C": 2, "D": 2}
    sampled_labels = [row[1] for row in sampled_rows]
    assert sampled_labels.count("A") == 2
    assert sampled_labels.count("B") == 2
    assert sampled_labels.count("C") == 2
    assert sampled_labels.count("D") == 2

    sampled_ids_repeat, _, counts_repeat = sample_shared_cells(
        header,
        rows,
        total_cells=8,
        stratify_key="celltype.l1",
        random_seed=7,
    )
    assert sampled_ids_repeat == sampled_ids
    assert counts_repeat == counts
