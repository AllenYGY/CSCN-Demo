from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]

DATA_SET = "GSE132465"
RUN_SLUG = "tumor_epithelial_vs_normal_epithelial"

GROUP_CONFIG = {
    "normal": {
        "class_labels": {"Normal"},
        "cell_types": {"Epithelial cells"},
        "disease": "normal",
    },
    "tumor": {
        "class_labels": {"Tumor"},
        "cell_types": {"Epithelial cells"},
        "disease": "tumor",
    },
}


def log(message):
    print(f"[{DATA_SET}] {message}", flush=True)


def log_stage(title):
    print(flush=True)
    print(f"=== {title} ===", flush=True)


def resolve_default_data_dir():
    candidates = [
        REPO_ROOT / "data" / DATA_SET,
        REPO_ROOT / "data" / "GSE32465",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def validate_required_file(path: Path, description: str):
    if not path.exists():
        raise FileNotFoundError(f"Missing {description}: {path}")
    log(f"{description}: {path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Prepare GSE132465 paired pseudobulk DESeq2 inputs for tumor epithelial "
            "cells versus normal epithelial cells."
        )
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=resolve_default_data_dir(),
        help="Dataset directory containing GSE132465 GEO files.",
    )
    parser.add_argument(
        "--min-cells-per-sample",
        type=int,
        default=20,
        help="Minimum epithelial cells required per sample to keep a matched pair.",
    )
    return parser.parse_args()


def load_selected_annotations(annotation_path: Path):
    selected_by_cell = {}
    summary_by_group = {}

    with gzip.open(annotation_path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Index", "Patient", "Class", "Sample", "Cell_type", "Cell_subtype"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"Missing required annotation columns in {annotation_path}: {sorted(missing)}"
            )

        for group_name in GROUP_CONFIG:
            summary_by_group[group_name] = {
                "total_cells": 0,
                "patient_counts": Counter(),
                "sample_counts": Counter(),
                "cell_subtype_counts": Counter(),
            }

        for row in reader:
            for group_name, config in GROUP_CONFIG.items():
                if row["Class"] not in config["class_labels"]:
                    continue
                if row["Cell_type"] not in config["cell_types"]:
                    continue

                cell_id = row["Index"]
                selected_by_cell[cell_id] = {
                    "cell_id": cell_id,
                    "patient": row["Patient"],
                    "class_label": row["Class"],
                    "sample": row["Sample"],
                    "cell_type": row["Cell_type"],
                    "cell_subtype": row["Cell_subtype"],
                    "group": group_name,
                    "disease": config["disease"],
                }
                summary = summary_by_group[group_name]
                summary["total_cells"] += 1
                summary["patient_counts"][row["Patient"]] += 1
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
    selected_cells_set = set()
    for position, cell_id in enumerate(matrix_cells):
        row = selected_by_cell.get(cell_id)
        if row is not None:
            selected_rows.append(row)
            selected_positions.append(position)
            selected_cells_set.add(cell_id)

    missing_cells = [cell_id for cell_id in selected_by_cell if cell_id not in selected_cells_set]
    if missing_cells:
        raise ValueError(
            f"Selected annotation cells missing from counts matrix: {missing_cells[:5]}"
        )

    return selected_rows, np.array(selected_positions, dtype=np.int64)


def build_paired_sample_groups(covariates, min_cells_per_sample):
    groups = {}
    for index, row in enumerate(covariates):
        label = row["sample"]
        if label not in groups:
            groups[label] = {
                "label": label,
                "sample": row["sample"],
                "patient": row["patient"],
                "class_label": row["class_label"],
                "disease": row["disease"],
                "cell_subtypes": [],
                "indices": [],
                "n_cells": 0,
            }

        group = groups[label]
        for key in ("patient", "class_label", "disease"):
            if group[key] != row[key]:
                raise ValueError(
                    f"Sample {label} mixes {key} values {group[key]} and {row[key]}"
                )
        if row["cell_subtype"] not in group["cell_subtypes"]:
            group["cell_subtypes"].append(row["cell_subtype"])
        group["indices"].append(index)
        group["n_cells"] += 1

    patient_to_samples = defaultdict(dict)
    for sample in groups.values():
        if sample["n_cells"] < min_cells_per_sample:
            continue
        patient_to_samples[sample["patient"]][sample["disease"]] = sample

    kept_patients = []
    ordered_groups = []
    for patient in sorted(patient_to_samples):
        sample_map = patient_to_samples[patient]
        if set(sample_map) != {"normal", "tumor"}:
            continue
        kept_patients.append(patient)
        for disease in ("normal", "tumor"):
            group = sample_map[disease]
            group["indices"] = np.array(group["indices"], dtype=np.int64)
            ordered_groups.append(group)

    excluded_samples = []
    for sample in groups.values():
        if sample["patient"] in kept_patients:
            continue
        reason = []
        if sample["n_cells"] < min_cells_per_sample:
            reason.append(f"n_cells<{min_cells_per_sample}")
        if sample["patient"] not in kept_patients:
            reason.append("unpaired_after_filter")
        excluded_samples.append(
            {
                "sample": sample["sample"],
                "patient": sample["patient"],
                "disease": sample["disease"],
                "n_cells": sample["n_cells"],
                "reason": "+".join(reason),
            }
        )

    return ordered_groups, kept_patients, excluded_samples


def restrict_to_kept_positions(selected_positions, groups):
    kept_indices = np.sort(np.concatenate([group["indices"] for group in groups]))
    kept_positions = selected_positions[kept_indices]
    index_map = {old_index: new_index for new_index, old_index in enumerate(kept_indices.tolist())}

    for group in groups:
        group["indices"] = np.array(
            [index_map[int(old_index)] for old_index in group["indices"]],
            dtype=np.int64,
        )

    return kept_positions


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
                "patient": group["patient"],
                "disease": group["disease"],
                "class_label": group["class_label"],
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
    annotation_path = data_dir / "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
    counts_path = data_dir / "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
    covariates_path = data_dir / f"{DATA_SET}_{RUN_SLUG}_covariates.csv.gz"
    metadata_path = output_dir / f"metadata_{RUN_SLUG}_paired_by_sample.csv"
    count_matrix_path = output_dir / f"count_matrix_{RUN_SLUG}_paired_by_sample.csv"
    excluded_samples_path = output_dir / f"excluded_samples_{RUN_SLUG}_paired_by_sample.csv"

    log_stage("Configuration")
    log(f"dataset dir: {data_dir}")
    log(f"output dir: {output_dir}")
    log(f"min cells per sample: {args.min_cells_per_sample}")
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
        "patient",
        "class_label",
        "sample",
        "cell_type",
        "cell_subtype",
        "group",
        "disease",
    ]
    write_gzip_csv(covariates_path, covariate_columns, covariates)
    log(f"saved covariates: {covariates_path}")

    log_stage("Build Paired Sample-Level Pseudobulk")
    groups, kept_patients, excluded_samples = build_paired_sample_groups(
        covariates,
        min_cells_per_sample=args.min_cells_per_sample,
    )
    if not groups:
        raise ValueError("No matched paired samples remain after filtering.")
    write_csv(
        metadata_path,
        ["sample", "patient", "disease", "class_label", "cell_subtypes", "n_cells"],
        metadata_rows(groups),
    )
    write_csv(
        excluded_samples_path,
        ["sample", "patient", "disease", "n_cells", "reason"],
        excluded_samples,
    )
    selected_positions = restrict_to_kept_positions(selected_positions, groups)
    log(f"kept patients: {kept_patients}")
    log(f"paired pseudo-samples: {len(groups)}")
    log(f"cells retained for pseudobulk: {len(selected_positions)}")
    log(f"saved metadata: {metadata_path}")
    log(f"saved excluded samples table: {excluded_samples_path}")

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
