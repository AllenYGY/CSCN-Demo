from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from cscn.aggregate import load_group_dags
from cscn.core import CSCN
from tests.support.fake_graph import FakeGraph  # noqa: F401  # ensure pickle compat in tests

from scripts.analysis.compare_ckm_clustering import (
    _group_cells_path,
    _group_object_path,
    _load_gene_names,
    _log,
    analyze_representation,
)


def _ckm_path(run_dir: Path, group_key: str) -> Path:
    return run_dir / "ckm" / f"{group_key}_ckm.npy"


def _matrix_path(run_dir: Path, group_key: str) -> Path:
    return run_dir / "matrices" / f"{group_key}.npy"


def _load_run_metadata(run_dir: Path) -> pd.DataFrame:
    metadata_path = run_dir / "inputs" / "cell_metadata.csv"
    metadata = pd.read_csv(metadata_path)
    metadata["cell_id"] = metadata["cell_id"].astype(str)
    return metadata.set_index("cell_id")


def _load_group_cells(run_dir: Path, group_key: str) -> list[str]:
    return pd.read_csv(_group_cells_path(run_dir, group_key))["cell_id"].astype(str).tolist()


def _available_dag_ids(run_dir: Path, group_key: str) -> set[int]:
    dag_dir = run_dir / "dags" / group_key
    return {int(path.stem.split("_", 1)[1]) for path in dag_dir.glob("result_*.pkl")}


def _available_cell_ids_by_run(run_dir: Path, group_key: str) -> tuple[list[str], dict[str, int]]:
    cells = _load_group_cells(run_dir, group_key)
    available_ids = _available_dag_ids(run_dir, group_key)
    available_pairs = [(cell_id, dag_idx) for dag_idx, cell_id in enumerate(cells) if dag_idx in available_ids]
    ordered_cell_ids = [cell_id for cell_id, _ in available_pairs]
    dag_index_by_cell = {cell_id: dag_idx for cell_id, dag_idx in available_pairs}
    return ordered_cell_ids, dag_index_by_cell


def _load_or_compute_ckm(
    run_dir: Path,
    group_key: str,
    *,
    beta_transform: str = "auto",
    strict: bool = False,
) -> pd.DataFrame:
    gene_names = _load_gene_names(run_dir)
    cells = _load_group_cells(run_dir, group_key)
    ckm_path = _ckm_path(run_dir, group_key)

    if ckm_path.is_file():
        return pd.DataFrame(__import__("numpy").load(ckm_path), index=cells, columns=gene_names)

    object_path = _group_object_path(run_dir, group_key)
    if object_path.is_file():
        cscn = CSCN.load_from_file(object_path)
    else:
        cscn = CSCN(output_dir=str(run_dir / "dags" / group_key))
        matrix = __import__("numpy").load(_matrix_path(run_dir, group_key))
        cscn.run_core(matrix, usingNMF=False)

    dags = load_group_dags(run_dir / "dags" / group_key)
    ckm = cscn.compute_ckm(
        dags=dags,
        beta_transform=beta_transform,
        save_path=str(ckm_path),
        strict=strict,
    )
    return pd.DataFrame(ckm, index=cells, columns=gene_names)


def _load_full_expression_subset(expr_path: Path, cell_ids: list[str]) -> pd.DataFrame:
    expr = pd.read_csv(expr_path, sep="\t")
    if "feature_name" not in expr.columns:
        raise ValueError("Expression baseline requires a `feature_name` column.")
    expr["feature_name"] = expr["feature_name"].astype(str)
    expr = expr.set_index("feature_name").T
    expr.index = expr.index.astype(str)

    missing_cells = [cell_id for cell_id in cell_ids if cell_id not in expr.index]
    if missing_cells:
        raise ValueError(f"Expression baseline is missing shared cells: {missing_cells[:10]}")

    return expr.loc[cell_ids].astype(float).copy()


def run_analysis(
    *,
    rna_expr_path: Path,
    adt_expr_path: Path,
    joint_expr_path: Path,
    rna_run_dir: Path,
    adt_run_dir: Path,
    joint_run_dir: Path,
    output_dir: Path,
    label_column: str = "hto_best_label",
    group_key: str = "all",
    random_seed: int = 42,
    enable_umap: bool = True,
    n_pcs: int = 2,
    ckm_beta_transform: str = "auto",
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    canonical_metadata = _load_run_metadata(rna_run_dir)
    if label_column not in canonical_metadata.columns:
        raise ValueError(f"Missing label column `{label_column}` in {rna_run_dir}/inputs/cell_metadata.csv")

    rna_cell_ids, rna_dag_by_cell = _available_cell_ids_by_run(rna_run_dir, group_key)
    adt_cell_ids, adt_dag_by_cell = _available_cell_ids_by_run(adt_run_dir, group_key)
    joint_cell_ids, joint_dag_by_cell = _available_cell_ids_by_run(joint_run_dir, group_key)

    shared_cell_ids = set(rna_cell_ids) & set(adt_cell_ids) & set(joint_cell_ids)
    if not shared_cell_ids:
        raise ValueError("No shared cell ids across the provided runs.")

    common_cell_ids = [cell_id for cell_id in rna_cell_ids if cell_id in shared_cell_ids]
    metadata = canonical_metadata.loc[common_cell_ids].copy()
    truth_labels = metadata[label_column].astype(str)

    expr_rna = _load_full_expression_subset(rna_expr_path, common_cell_ids)
    expr_adt = _load_full_expression_subset(adt_expr_path, common_cell_ids)
    expr_joint = _load_full_expression_subset(joint_expr_path, common_cell_ids)

    ckm_rna = _load_or_compute_ckm(
        rna_run_dir,
        group_key,
        beta_transform=ckm_beta_transform,
        strict=False,
    ).loc[common_cell_ids]
    ckm_adt = _load_or_compute_ckm(
        adt_run_dir,
        group_key,
        beta_transform=ckm_beta_transform,
        strict=False,
    ).loc[common_cell_ids]
    ckm_joint = _load_or_compute_ckm(
        joint_run_dir,
        group_key,
        beta_transform=ckm_beta_transform,
        strict=False,
    ).loc[common_cell_ids]

    representations = {
        "expr_rna": expr_rna,
        "expr_adt": expr_adt,
        "expr_joint": expr_joint,
        "ckm_rna": ckm_rna,
        "ckm_adt": ckm_adt,
        "ckm_joint": ckm_joint,
    }

    metrics_rows: list[dict[str, float | int | str]] = []
    assignment_frame = pd.DataFrame(
        {
            "cell_id": metadata.index,
            label_column: truth_labels.values,
        }
    ).set_index("cell_id")

    for name, frame in representations.items():
        metrics, clusters = analyze_representation(
            name,
            frame,
            truth_labels,
            output_dir,
            random_seed=random_seed,
            n_pcs=n_pcs,
            enable_umap=enable_umap,
        )
        metrics_rows.append(metrics)
        assignment_frame[clusters.name] = clusters

    metrics_df = pd.DataFrame(metrics_rows)
    metrics_df.to_csv(output_dir / "clustering_metrics.csv", index=False)
    assignment_frame.reset_index().to_csv(output_dir / "cell_assignments.csv", index=False)
    pd.DataFrame(
        {
            "cell_id": common_cell_ids,
            "rna_dag_index": [rna_dag_by_cell[cell_id] for cell_id in common_cell_ids],
            "adt_dag_index": [adt_dag_by_cell[cell_id] for cell_id in common_cell_ids],
            "joint_dag_index": [joint_dag_by_cell[cell_id] for cell_id in common_cell_ids],
        }
    ).to_csv(output_dir / "shared_cells.csv", index=False)
    _log(f"saved GSE128639 modality comparison outputs to {output_dir}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare GSE128639 expression baselines and CSCN CKM representations across modalities."
    )
    parser.add_argument(
        "--rna-expr-path",
        type=Path,
        default=REPO_ROOT / "data" / "GSE128639" / "cscn_inputs" / "gse128639_mnc_shared3000_rna_only_expression.tsv.gz",
    )
    parser.add_argument(
        "--adt-expr-path",
        type=Path,
        default=REPO_ROOT / "data" / "GSE128639" / "cscn_inputs" / "gse128639_mnc_shared3000_adt_only_expression.tsv.gz",
    )
    parser.add_argument(
        "--joint-expr-path",
        type=Path,
        default=REPO_ROOT / "data" / "GSE128639" / "cscn_inputs" / "gse128639_mnc_shared3000_rna_adt_joint_expression.tsv.gz",
    )
    parser.add_argument(
        "--rna-run-dir",
        type=Path,
        default=REPO_ROOT / "runs" / "gse128639_mnc_shared3000_rna_only",
    )
    parser.add_argument(
        "--adt-run-dir",
        type=Path,
        default=REPO_ROOT / "runs" / "gse128639_mnc_shared3000_adt_only",
    )
    parser.add_argument(
        "--joint-run-dir",
        type=Path,
        default=REPO_ROOT / "runs" / "gse128639_mnc_shared3000_rna_adt_joint",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "results" / "gse128639_modality_clustering_shared3000",
    )
    parser.add_argument("--label-column", default="hto_best_label")
    parser.add_argument("--group-key", default="all")
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--n-pcs", type=int, default=2)
    parser.add_argument(
        "--ckm-beta-transform",
        choices=("auto", "log1p", "identity"),
        default="auto",
    )
    parser.add_argument("--umap", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    run_analysis(
        rna_expr_path=args.rna_expr_path,
        adt_expr_path=args.adt_expr_path,
        joint_expr_path=args.joint_expr_path,
        rna_run_dir=args.rna_run_dir,
        adt_run_dir=args.adt_run_dir,
        joint_run_dir=args.joint_run_dir,
        output_dir=args.output_dir,
        label_column=args.label_column,
        group_key=args.group_key,
        random_seed=args.random_seed,
        enable_umap=args.umap,
        n_pcs=args.n_pcs,
        ckm_beta_transform=args.ckm_beta_transform,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
