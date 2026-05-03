from __future__ import annotations

import argparse
from pathlib import Path

from .compare_ckm_clustering import run_analysis as run_generic_analysis


def run_analysis(
    expr_path: Path,
    metadata_path: Path,
    weighted_run_dir: Path,
    local_knn_run_dir: Path,
    output_dir: Path,
    *,
    random_seed: int = 42,
    enable_umap: bool = True,
    n_pcs: int = 2,
    ckm_beta_transform: str = "auto",
) -> None:
    run_generic_analysis(
        expr_path=expr_path,
        metadata_path=metadata_path,
        weighted_run_dir=weighted_run_dir,
        local_knn_run_dir=local_knn_run_dir,
        output_dir=output_dir,
        label_column="cell_class_name",
        expr_orientation="cells_by_genes",
        expr_cell_id_column="cell_id",
        metadata_cell_id_column="cell_id",
        random_seed=random_seed,
        enable_umap=enable_umap,
        n_pcs=n_pcs,
        ckm_beta_transform=ckm_beta_transform,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compare seqFISH expression and CKM clustering.")
    parser.add_argument("--expr-path", type=Path, required=True)
    parser.add_argument("--metadata-path", type=Path, required=True)
    parser.add_argument("--weighted-run-dir", type=Path, required=True)
    parser.add_argument("--local-knn-run-dir", type=Path, required=True)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("config/seqfish/analysis_ckm_compare"),
    )
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
        expr_path=args.expr_path,
        metadata_path=args.metadata_path,
        weighted_run_dir=args.weighted_run_dir,
        local_knn_run_dir=args.local_knn_run_dir,
        output_dir=args.output_dir,
        random_seed=args.random_seed,
        enable_umap=args.umap,
        n_pcs=args.n_pcs,
        ckm_beta_transform=args.ckm_beta_transform,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
