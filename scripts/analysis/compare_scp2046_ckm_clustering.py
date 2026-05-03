from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .compare_ckm_clustering import run_analysis as run_generic_analysis


def _resolve_scp2046_metadata(metadata_path: Path, output_dir: Path) -> Path:
    metadata = pd.read_csv(metadata_path)
    if "NAME" in metadata.columns:
        metadata = metadata.loc[metadata["NAME"].astype(str) != "TYPE"].copy()
    if "cell_id" in metadata.columns:
        metadata = metadata.loc[metadata["cell_id"].astype(str) != "TYPE"].copy()
    if "cell_id" in metadata.columns and "zones" in metadata.columns:
        return metadata_path

    if "NAME" in metadata.columns and "zones" in metadata.columns:
        resolved = metadata.rename(columns={"NAME": "cell_id"})
        output_dir.mkdir(parents=True, exist_ok=True)
        resolved_path = output_dir / "_resolved_scp2046_metadata.csv"
        resolved.to_csv(resolved_path, index=False)
        return resolved_path

    if "cell_id" not in metadata.columns:
        if "NAME" in metadata.columns:
            metadata = metadata.rename(columns={"NAME": "cell_id"})
        else:
            raise ValueError(
                "SCP2046 comparison metadata must contain either `cell_id` or `NAME`."
            )

    repo_root = Path(__file__).resolve().parents[2]
    canonical_path = repo_root / "data" / "SCP2046" / "metadata" / "meta_data.csv"
    canonical = pd.read_csv(canonical_path, skiprows=[1]).rename(columns={"NAME": "cell_id"})
    canonical["cell_id"] = canonical["cell_id"].astype(str)
    metadata["cell_id"] = metadata["cell_id"].astype(str)

    merged = metadata.merge(
        canonical[["cell_id", "biosample_id", "zones", "disease", "donor_id"]],
        on="cell_id",
        how="left",
        validate="one_to_one",
    )
    if merged["zones"].isnull().any():
        missing = merged.loc[merged["zones"].isnull(), "cell_id"].head(10).tolist()
        raise ValueError(
            f"Could not recover `zones` for some SCP2046 cells from canonical metadata: {missing}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    resolved_path = output_dir / "_resolved_scp2046_metadata.csv"
    merged.to_csv(resolved_path, index=False)
    return resolved_path


def run_analysis(
    expr_path: Path,
    metadata_path: Path,
    weighted_run_dir: Path,
    comparison_run_dir: Path,
    output_dir: Path,
    *,
    comparison_label: str = "ckm_local_knn",
    extra_run_dir: Path | None = None,
    extra_label: str | None = None,
    random_seed: int = 42,
    enable_umap: bool = True,
    n_pcs: int = 2,
    ckm_beta_transform: str = "auto",
) -> None:
    resolved_metadata_path = _resolve_scp2046_metadata(metadata_path, output_dir)
    run_generic_analysis(
        expr_path=expr_path,
        metadata_path=resolved_metadata_path,
        weighted_run_dir=weighted_run_dir,
        local_knn_run_dir=comparison_run_dir,
        output_dir=output_dir,
        label_column="zones",
        weighted_representation_name="ckm_weighted",
        local_knn_representation_name=comparison_label,
        extra_run_dir=extra_run_dir,
        extra_representation_name=extra_label,
        expr_orientation="genes_by_cells",
        expr_gene_key="GENE",
        expr_cell_id_column="cell_id",
        metadata_cell_id_column="cell_id",
        random_seed=random_seed,
        enable_umap=enable_umap,
        n_pcs=n_pcs,
        ckm_beta_transform=ckm_beta_transform,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compare SCP2046 expression and CKM clustering.")
    parser.add_argument("--expr-path", type=Path, required=True)
    parser.add_argument("--metadata-path", type=Path, required=True)
    parser.add_argument("--weighted-run-dir", type=Path, required=True)
    parser.add_argument("--comparison-run-dir", type=Path, required=True)
    parser.add_argument("--comparison-label", default="ckm_local_knn")
    parser.add_argument("--extra-run-dir", type=Path)
    parser.add_argument("--extra-label")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("runs/scp2046_sham1_compare"),
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
        comparison_run_dir=args.comparison_run_dir,
        output_dir=args.output_dir,
        comparison_label=args.comparison_label,
        extra_run_dir=args.extra_run_dir,
        extra_label=args.extra_label,
        random_seed=args.random_seed,
        enable_umap=args.umap,
        n_pcs=args.n_pcs,
        ckm_beta_transform=args.ckm_beta_transform,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
