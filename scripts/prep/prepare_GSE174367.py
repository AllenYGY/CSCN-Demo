from __future__ import annotations

import argparse
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from gse174367_utils import (
    DATA_SET,
    DEFAULT_CASE_GROUP,
    DEFAULT_CONTROL_GROUP,
    build_pseudobulk_count_matrix,
    build_run_slug,
    build_sample_metadata,
    filter_metadata_for_cell_type,
    load_gse174367_metadata,
)


def log(message: str) -> None:
    print(f"[{DATA_SET}] {message}", flush=True)


def log_stage(title: str) -> None:
    print(flush=True)
    print(f"=== {title} ===", flush=True)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Prepare donor-level pseudobulk DESeq2 inputs for GSE174367 "
            "within a selected cell type."
        )
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=REPO_ROOT / "data" / DATA_SET,
        help="Dataset directory containing the GSE174367 snRNA-seq files.",
    )
    parser.add_argument(
        "--cell-type",
        default="ODC",
        help="Cell.Type value to prepare, for example ODC, MG, ASC, EX, INH, or OPC.",
    )
    parser.add_argument(
        "--case-group",
        default=DEFAULT_CASE_GROUP,
        help="Case diagnosis label. Default: AD",
    )
    parser.add_argument(
        "--control-group",
        default=DEFAULT_CONTROL_GROUP,
        help="Control diagnosis label. Default: Control",
    )
    parser.add_argument(
        "--min-cells-per-sample",
        type=int,
        default=20,
        help="Minimum number of cells required per donor after cell-type filtering.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.min_cells_per_sample <= 0:
        raise ValueError("--min-cells-per-sample must be a positive integer")

    data_dir = args.data_dir.resolve()
    metadata_path = data_dir / "GSE174367_snRNA-seq_cell_meta.csv.gz"
    h5_path = data_dir / "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5"
    output_dir = data_dir / "output_deseq"
    run_slug = build_run_slug(
        cell_type=args.cell_type,
        case_group=args.case_group,
        control_group=args.control_group,
    )
    count_path = output_dir / f"count_matrix_{run_slug}_by_sample.csv"
    metadata_out_path = output_dir / f"metadata_{run_slug}_by_sample.csv"
    selected_cells_path = output_dir / f"covariates_{run_slug}_selected_cells.csv"

    log_stage("Configuration")
    log(f"dataset dir: {data_dir}")
    log(f"cell type: {args.cell_type}")
    log(f"case group: {args.case_group}")
    log(f"control group: {args.control_group}")
    log(f"min cells per sample: {args.min_cells_per_sample}")
    log(f"count output: {count_path}")
    log(f"metadata output: {metadata_out_path}")

    log_stage("Load Metadata")
    metadata_df = load_gse174367_metadata(metadata_path)
    selected_df = filter_metadata_for_cell_type(
        metadata_df,
        cell_type=args.cell_type,
        case_group=args.case_group,
        control_group=args.control_group,
        min_cells_per_sample=args.min_cells_per_sample,
    )
    sample_metadata_df = build_sample_metadata(selected_df)
    group_summary = sample_metadata_df.groupby("diagnosis")["sample"].count().to_dict()
    cell_summary = sample_metadata_df.groupby("diagnosis")["n_cells"].sum().to_dict()
    for diagnosis in (args.control_group, args.case_group):
        log(
            f"{diagnosis}: donors={group_summary.get(diagnosis, 0)}, "
            f"cells={cell_summary.get(diagnosis, 0)}"
        )
    if not all(group_summary.get(group, 0) > 0 for group in (args.control_group, args.case_group)):
        raise ValueError(
            "Both diagnosis groups must retain at least one donor after filtering."
        )

    log_stage("Build Pseudobulk Counts")
    count_df = build_pseudobulk_count_matrix(h5_path, selected_df, sample_metadata_df)
    output_dir.mkdir(parents=True, exist_ok=True)
    selected_df.to_csv(selected_cells_path, index=False)
    sample_metadata_df.to_csv(metadata_out_path, index=False)
    count_df.to_csv(count_path, index=False)
    log(f"saved selected-cell metadata: {selected_cells_path}")
    log(f"saved sample metadata: {metadata_out_path}")
    log(f"saved pseudobulk count matrix: {count_path}")
    log(f"genes written: {count_df.shape[0]}")
    log(f"donor columns written: {count_df.shape[1] - 1}")


if __name__ == "__main__":
    main()
