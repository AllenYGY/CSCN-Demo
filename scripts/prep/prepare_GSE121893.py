from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]

DATA_SET = "GSE121893"
MAIN_CASE_GROUP = "dHF"
MAIN_CONTROL_GROUP = "N"


def log(message):
    print(f"[{DATA_SET}] {message}")


def log_stage(title):
    print()
    print(f"=== {title} ===")


def validate_required_file(path: Path, description: str):
    if not path.exists():
        raise FileNotFoundError(f"Missing {description}: {path}")
    log(f"{description}: {path}")


def unique_preserve_order(values):
    seen = set()
    ordered = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Prepare GSE121893 cell covariates and pseudo-bulk DESeq2 inputs "
            "from the GEO raw matrix plus cluster annotations."
        )
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=REPO_ROOT / "data" / DATA_SET,
        help="Dataset directory containing the downloaded GSE121893 GEO files.",
    )
    return parser.parse_args()


def iter_gzip_tsv_rows(path: Path):
    with gzip.open(path, "rt", newline="") as handle:
        header = None
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            header = line.rstrip("\n").split("\t")
            break
        if header is None:
            raise ValueError(f"Could not locate a TSV header in {path}")
        reader = csv.DictReader(handle, fieldnames=header, delimiter="\t")
        yield from reader


def read_count_cells(counts_path: Path):
    with gzip.open(counts_path, "rt", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
    if len(header) <= 1:
        raise ValueError(f"Counts matrix header is empty: {counts_path}")
    return header[1:]


def parse_region_condition(value: str):
    disease, region = value.split("_", 1)
    return disease, region


def load_filtered_covariates(cluster_info_path: Path, count_cells):
    count_cell_set = set(count_cells)
    row_by_cell = {}

    for row in iter_gzip_tsv_rows(cluster_info_path):
        cell_id = row["ID"]
        if cell_id not in count_cell_set:
            continue

        disease, region = parse_region_condition(row["condition"])
        row_by_cell[cell_id] = {
            "cell_id": cell_id,
            "sample": row["sample"],
            "condition": row["condition"],
            "disease": disease,
            "region": region,
            "group": row["group"],
            "cell_type": row["ident"],
            "age": row.get("Age", ""),
        }

    missing = [cell_id for cell_id in count_cells if cell_id not in row_by_cell]
    if missing:
        raise ValueError(
            f"Missing {len(missing)} count-matrix cells from cluster annotations; "
            f"first 5 missing: {missing[:5]}"
        )

    return [row_by_cell[cell_id] for cell_id in count_cells]


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


def build_pseudobulk_groups(covariates, mode: str):
    groups = {}
    ordered_labels = []

    for index, row in enumerate(covariates):
        if mode == "sample_region":
            label = f"{row['sample']}_{row['region']}"
            if label not in groups:
                groups[label] = {
                    "label": label,
                    "sample": label,
                    "donor_sample": row["sample"],
                    "disease": row["disease"],
                    "region": row["region"],
                    "indices": [],
                    "n_cells": 0,
                }
                ordered_labels.append(label)
        elif mode == "sample":
            label = row["sample"]
            if label not in groups:
                groups[label] = {
                    "label": label,
                    "sample": row["sample"],
                    "disease": row["disease"],
                    "regions": [],
                    "indices": [],
                    "n_cells": 0,
                }
                ordered_labels.append(label)
            if groups[label]["disease"] != row["disease"]:
                raise ValueError(
                    f"Sample {label} mixes disease labels "
                    f"{groups[label]['disease']} and {row['disease']}"
                )
            if row["region"] not in groups[label]["regions"]:
                groups[label]["regions"].append(row["region"])
        else:
            raise ValueError(f"Unsupported pseudobulk mode: {mode}")

        groups[label]["indices"].append(index)
        groups[label]["n_cells"] += 1

    ordered_groups = []
    for label in ordered_labels:
        group = groups[label]
        group["indices"] = np.array(group["indices"], dtype=np.int64)
        ordered_groups.append(group)
    return ordered_groups


def build_metadata_rows(groups, mode: str):
    rows = []
    if mode == "sample_region":
        for group in groups:
            rows.append(
                {
                    "sample": group["sample"],
                    "donor_sample": group["donor_sample"],
                    "disease": group["disease"],
                    "region": group["region"],
                    "n_cells": group["n_cells"],
                }
            )
    elif mode == "sample":
        for group in groups:
            rows.append(
                {
                    "sample": group["sample"],
                    "disease": group["disease"],
                    "regions": "+".join(group["regions"]),
                    "n_cells": group["n_cells"],
                }
            )
    else:
        raise ValueError(f"Unsupported metadata mode: {mode}")
    return rows


def filter_groups_by_disease(groups, allowed_diseases):
    allowed = set(allowed_diseases)
    return [group for group in groups if group["disease"] in allowed]


def write_multiple_count_matrices(counts_path: Path, matrix_specs):
    if not matrix_specs:
        return 0

    prepared_specs = []
    handles = []
    try:
        for spec in matrix_specs:
            output_path = spec["output_path"]
            output_path.parent.mkdir(parents=True, exist_ok=True)
            handle = open(output_path, "w", newline="")
            writer = csv.writer(handle)
            writer.writerow(["gene", *[group["label"] for group in spec["groups"]]])
            prepared_specs.append(
                {
                    "writer": writer,
                    "index_arrays": [group["indices"] for group in spec["groups"]],
                }
            )
            handles.append(handle)

        gene_count = 0
        with gzip.open(counts_path, "rt", newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader)
            n_cells = len(header) - 1

            for row in reader:
                gene = row[0]
                values = np.fromiter(
                    (int(value) for value in row[1:]),
                    dtype=np.int64,
                    count=n_cells,
                )
                for spec in prepared_specs:
                    sums = [int(values[index_array].sum()) for index_array in spec["index_arrays"]]
                    spec["writer"].writerow([gene, *sums])
                gene_count += 1
        return gene_count
    finally:
        for handle in handles:
            handle.close()


def main():
    args = parse_args()
    data_dir = args.data_dir.resolve()
    output_dir = data_dir / "output_deseq"
    counts_path = data_dir / "GSE121893_human_heart_sc_umi.csv.gz"
    cluster_info_path = data_dir / "GSE121893_all_heart_cell_cluster_info.txt.gz"
    covariates_path = data_dir / "GSE121893_covariates.csv.gz"

    log_stage("Configuration")
    log(f"dataset dir: {data_dir}")
    log(f"output dir: {output_dir}")
    validate_required_file(counts_path, "counts matrix")
    validate_required_file(cluster_info_path, "cluster annotation table")

    log_stage("Load Cells")
    count_cells = read_count_cells(counts_path)
    log(f"cells in counts matrix: {len(count_cells)}")

    covariates = load_filtered_covariates(cluster_info_path, count_cells)
    disease_counts = Counter(row["disease"] for row in covariates)
    sample_counts = Counter(row["sample"] for row in covariates)
    cell_type_counts = Counter(row["cell_type"] for row in covariates)
    log(f"disease counts: {dict(disease_counts)}")
    log(f"samples in counts matrix: {dict(sample_counts)}")
    log(f"top 10 cell types: {dict(cell_type_counts.most_common(10))}")

    log_stage("Write Covariates")
    covariate_columns = [
        "cell_id",
        "sample",
        "condition",
        "disease",
        "region",
        "group",
        "cell_type",
        "age",
    ]
    write_gzip_csv(covariates_path, covariate_columns, covariates)
    log(f"saved covariates: {covariates_path}")

    log_stage("Build Pseudobulk Metadata")
    sample_region_groups = build_pseudobulk_groups(covariates, mode="sample_region")
    sample_groups = build_pseudobulk_groups(covariates, mode="sample")
    main_diseases = {MAIN_CASE_GROUP, MAIN_CONTROL_GROUP}
    sample_region_groups_main = filter_groups_by_disease(sample_region_groups, main_diseases)
    sample_groups_main = filter_groups_by_disease(sample_groups, main_diseases)

    metadata_specs = [
        {
            "path": output_dir / "metadata.csv",
            "rows": build_metadata_rows(sample_region_groups, mode="sample_region"),
            "fieldnames": ["sample", "donor_sample", "disease", "region", "n_cells"],
        },
        {
            "path": output_dir / "metadata_by_sample.csv",
            "rows": build_metadata_rows(sample_groups, mode="sample"),
            "fieldnames": ["sample", "disease", "regions", "n_cells"],
        },
        {
            "path": output_dir / "metadata_dhf_vs_n.csv",
            "rows": build_metadata_rows(sample_region_groups_main, mode="sample_region"),
            "fieldnames": ["sample", "donor_sample", "disease", "region", "n_cells"],
        },
        {
            "path": output_dir / "metadata_dhf_vs_n_by_sample.csv",
            "rows": build_metadata_rows(sample_groups_main, mode="sample"),
            "fieldnames": ["sample", "disease", "regions", "n_cells"],
        },
    ]
    for spec in metadata_specs:
        write_csv(spec["path"], spec["fieldnames"], spec["rows"])
        log(f"saved metadata: {spec['path']}")

    log_stage("Write Pseudobulk Count Matrices")
    gene_count = write_multiple_count_matrices(
        counts_path=counts_path,
        matrix_specs=[
            {"output_path": output_dir / "count_matrix.csv", "groups": sample_region_groups},
            {
                "output_path": output_dir / "count_matrix_by_sample.csv",
                "groups": sample_groups,
            },
            {
                "output_path": output_dir / "count_matrix_dhf_vs_n.csv",
                "groups": sample_region_groups_main,
            },
            {
                "output_path": output_dir / "count_matrix_dhf_vs_n_by_sample.csv",
                "groups": sample_groups_main,
            },
        ],
    )
    log(f"genes written per matrix: {gene_count}")
    log(f"sample-region pseudo-samples: {len(sample_region_groups)}")
    log(f"sample-level pseudo-samples: {len(sample_groups)}")
    log(f"{MAIN_CASE_GROUP} vs {MAIN_CONTROL_GROUP} sample-region pseudo-samples: {len(sample_region_groups_main)}")
    log(f"{MAIN_CASE_GROUP} vs {MAIN_CONTROL_GROUP} sample-level pseudo-samples: {len(sample_groups_main)}")


if __name__ == "__main__":
    main()
