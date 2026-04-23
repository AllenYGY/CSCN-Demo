from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]

DATA_SET = "GSE131907"
RUN_SLUG = "nlung_epithelial_vs_tlung_malignant"

GROUP_CONFIG = {
    "normal": {
        "sample_origins": {"nLung"},
        "cell_types": {"Epithelial cells"},
        "cell_subtypes": {"AT1", "AT2", "Club", "Ciliated"},
        "disease": "normal",
    },
    "cancer": {
        "sample_origins": {"tLung"},
        "cell_types": {"Epithelial cells"},
        "cell_subtypes": {"tS1", "tS2", "tS3"},
        "disease": "cancer",
    },
}


def log(message):
    print(f"[{DATA_SET}] {message}")


def log_stage(title):
    print()
    print(f"=== {title} ===")


def validate_required_file(path: Path, description: str):
    if not path.exists():
        raise FileNotFoundError(f"Missing {description}: {path}")
    log(f"{description}: {path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Prepare GSE131907 pseudobulk DESeq2 inputs for "
            "nLung epithelial cells versus tLung malignant epithelial cells."
        )
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=REPO_ROOT / "data" / DATA_SET,
        help="Dataset directory containing GSE131907 GEO files.",
    )
    return parser.parse_args()


def load_selected_annotations(annotation_path: Path):
    selected_by_cell = {}
    summary_by_group = {}

    with gzip.open(annotation_path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Index", "Barcode", "Sample", "Sample_Origin", "Cell_type", "Cell_subtype"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"Missing required annotation columns in {annotation_path}: {sorted(missing)}"
            )

        for group_name in GROUP_CONFIG:
            summary_by_group[group_name] = {
                "total_cells": 0,
                "sample_counts": Counter(),
                "cell_subtype_counts": Counter(),
            }

        for row in reader:
            for group_name, config in GROUP_CONFIG.items():
                if row["Sample_Origin"] not in config["sample_origins"]:
                    continue
                if row["Cell_type"] not in config["cell_types"]:
                    continue
                if row["Cell_subtype"] not in config["cell_subtypes"]:
                    continue

                cell_id = row["Index"]
                selected_by_cell[cell_id] = {
                    "cell_id": cell_id,
                    "barcode": row["Barcode"],
                    "sample": row["Sample"],
                    "sample_origin": row["Sample_Origin"],
                    "cell_type": row["Cell_type"],
                    "cell_subtype": row["Cell_subtype"],
                    "group": group_name,
                    "disease": config["disease"],
                }
                summary = summary_by_group[group_name]
                summary["total_cells"] += 1
                summary["sample_counts"][row["Sample"]] += 1
                summary["cell_subtype_counts"][row["Cell_subtype"]] += 1
                break

    return selected_by_cell, summary_by_group


def build_covariates_in_matrix_order(counts_path: Path, selected_by_cell):
    with gzip.open(counts_path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)

    matrix_cells = header[1:]
    selected_rows = []
    selected_positions = []
    missing_cells = []
    for position, cell_id in enumerate(matrix_cells):
        row = selected_by_cell.get(cell_id)
        if row is not None:
            selected_rows.append(row)
            selected_positions.append(position)

    for cell_id in selected_by_cell:
        if cell_id not in set(cell["cell_id"] for cell in selected_rows):
            missing_cells.append(cell_id)

    if missing_cells:
        raise ValueError(
            f"Selected annotation cells missing from counts matrix: {missing_cells[:5]}"
        )

    return selected_rows, np.array(selected_positions, dtype=np.int64)


def build_sample_groups(covariates):
    groups = {}
    ordered_labels = []
    for index, row in enumerate(covariates):
        label = row["sample"]
        if label not in groups:
            groups[label] = {
                "label": label,
                "sample": row["sample"],
                "disease": row["disease"],
                "sample_origin": row["sample_origin"],
                "cell_subtypes": [],
                "indices": [],
                "n_cells": 0,
            }
            ordered_labels.append(label)

        group = groups[label]
        if group["disease"] != row["disease"]:
            raise ValueError(
                f"Sample {label} mixes disease labels {group['disease']} and {row['disease']}"
            )
        if row["cell_subtype"] not in group["cell_subtypes"]:
            group["cell_subtypes"].append(row["cell_subtype"])

        group["indices"].append(index)
        group["n_cells"] += 1

    ordered_groups = []
    for label in ordered_labels:
        group = groups[label]
        group["indices"] = np.array(group["indices"], dtype=np.int64)
        ordered_groups.append(group)
    return ordered_groups


def write_csv(path: Path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_gzip_csv(path: Path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def metadata_rows(groups):
    rows = []
    for group in groups:
        rows.append(
            {
                "sample": group["sample"],
                "disease": group["disease"],
                "sample_origin": group["sample_origin"],
                "cell_subtypes": "+".join(group["cell_subtypes"]),
                "n_cells": group["n_cells"],
            }
        )
    return rows


def write_pseudobulk_matrix(counts_path: Path, output_path: Path, selected_positions, groups):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    group_ids = np.empty(len(selected_positions), dtype=np.int64)
    for group_index, group in enumerate(groups):
        group_ids[group["indices"]] = group_index

    with open(output_path, "w", newline="") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["gene", *[group["label"] for group in groups]])

        gene_count = 0
        with gzip.open(counts_path, "rt", newline="") as in_handle:
            reader = csv.reader(in_handle, delimiter="\t")
            next(reader)
            n_selected = len(selected_positions)

            for row in reader:
                gene = row[0]
                selected_values = np.fromiter(
                    (int(row[position + 1]) for position in selected_positions),
                    dtype=np.int64,
                    count=n_selected,
                )
                sums = np.bincount(
                    group_ids,
                    weights=selected_values,
                    minlength=len(groups),
                ).astype(np.int64)
                writer.writerow([gene, *[int(value) for value in sums]])
                gene_count += 1

    return gene_count


def main():
    args = parse_args()
    data_dir = args.data_dir.resolve()
    output_dir = data_dir / "output_deseq"
    annotation_path = data_dir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    counts_path = data_dir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    covariates_path = data_dir / f"{DATA_SET}_{RUN_SLUG}_covariates.csv.gz"
    metadata_path = output_dir / f"metadata_{RUN_SLUG}_by_sample.csv"
    count_matrix_path = output_dir / f"count_matrix_{RUN_SLUG}_by_sample.csv"

    log_stage("Configuration")
    log(f"dataset dir: {data_dir}")
    log(f"output dir: {output_dir}")
    validate_required_file(annotation_path, "cell annotation")
    validate_required_file(counts_path, "raw UMI matrix")

    log_stage("Select Cells")
    selected_by_cell, summary_by_group = load_selected_annotations(annotation_path)
    for group_name, summary in summary_by_group.items():
        log(f"{group_name} selected cells: {summary['total_cells']}")
        log(f"{group_name} sample counts: {dict(summary['sample_counts'])}")
        log(f"{group_name} subtype counts: {dict(summary['cell_subtype_counts'])}")

    log_stage("Align To Counts Matrix")
    covariates, selected_positions = build_covariates_in_matrix_order(counts_path, selected_by_cell)
    log(f"selected cells present in matrix: {len(covariates)}")

    log_stage("Write Covariates")
    covariate_columns = [
        "cell_id",
        "barcode",
        "sample",
        "sample_origin",
        "cell_type",
        "cell_subtype",
        "group",
        "disease",
    ]
    write_gzip_csv(covariates_path, covariate_columns, covariates)
    log(f"saved covariates: {covariates_path}")

    log_stage("Build Sample-Level Pseudobulk")
    groups = build_sample_groups(covariates)
    rows = metadata_rows(groups)
    write_csv(metadata_path, ["sample", "disease", "sample_origin", "cell_subtypes", "n_cells"], rows)
    log(f"saved metadata: {metadata_path}")
    log(f"pseudo-samples: {len(groups)}")

    gene_count = write_pseudobulk_matrix(
        counts_path=counts_path,
        output_path=count_matrix_path,
        selected_positions=selected_positions,
        groups=groups,
    )
    log(f"saved count matrix: {count_matrix_path}")
    log(f"genes written: {gene_count}")


if __name__ == "__main__":
    main()
