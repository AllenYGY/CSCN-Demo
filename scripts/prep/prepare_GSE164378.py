from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_SET = "GSE164378"
RAW_TAR_NAME = "GSE164378_RAW.tar"
METADATA_3P_NAME = "GSE164378_sc.meta.data_3P.csv.gz"

RNA_3P_BARCODES = "GSM5008737_RNA_3P-barcodes.tsv.gz"
RNA_3P_FEATURES = "GSM5008737_RNA_3P-features.tsv.gz"
RNA_3P_MATRIX = "GSM5008737_RNA_3P-matrix.mtx.gz"

ADT_3P_BARCODES = "GSM5008738_ADT_3P-barcodes.tsv.gz"
ADT_3P_FEATURES = "GSM5008738_ADT_3P-features.tsv.gz"
ADT_3P_MATRIX = "GSM5008738_ADT_3P-matrix.mtx.gz"

RNA_ONLY_TOP_N = 150
ADT_ONLY_TOP_N = 150
JOINT_RNA_TOP_N = 100
JOINT_ADT_TOP_N = 50
RNA_TARGET_SUM = 1e6
CHUNK_ROWS = 2_000_000
SHARED_TOTAL_CELLS = 2000
SHARED_STRATIFY_KEY = "celltype.l1"
SHARED_RANDOM_SEED = 42


def log(message: str) -> None:
    print(f"[{DATA_SET}] {message}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare GSE164378 3P CSCN inputs for RNA-only, ADT-only, and "
            "RNA+ADT joint runs with fixed feature budgets."
        )
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=REPO_ROOT / "data" / DATA_SET,
        help="Dataset directory containing GSE164378_RAW.tar and 3P metadata.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory for CSCN-ready inputs. Defaults to <data-dir>/cscn_inputs.",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=None,
        help="Scratch directory for extracted 3P files. Defaults to <data-dir>/.gse164378_3p_work.",
    )
    return parser.parse_args()


def validate_required_file(path: Path, description: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {description}: {path}")
    log(f"{description}: {path}")


def read_gzip_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt", newline="") as handle:
        return [line.rstrip("\n") for line in handle if line.strip()]


def read_feature_names(path: Path) -> list[str]:
    names: list[str] = []
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if row:
                names.append(str(row[0]))
    return names


def read_metadata_rows(path: Path) -> tuple[list[str], list[str], list[list[str]]]:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        normalized_header = ["cell_id" if item == "" else str(item) for item in header]
        rows: list[list[str]] = []
        cell_ids: list[str] = []
        for row in reader:
            if not row:
                continue
            row = [str(item) for item in row]
            rows.append(row)
            cell_ids.append(row[0])
    return normalized_header, cell_ids, rows


def write_gzip_csv(path: Path, header: list[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def write_expression_table(
    path: Path,
    cell_ids: list[str],
    feature_names: list[str],
    matrix: np.ndarray,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cell_id", *feature_names])
        for row_index, cell_id in enumerate(cell_ids):
            writer.writerow([cell_id, *matrix[row_index].astype(np.float32).tolist()])


def write_text_lines(path: Path, values: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(values) + "\n", encoding="utf-8")


def sample_shared_cells(
    metadata_header: list[str],
    metadata_rows: list[list[str]],
    *,
    total_cells: int = SHARED_TOTAL_CELLS,
    stratify_key: str = SHARED_STRATIFY_KEY,
    random_seed: int = SHARED_RANDOM_SEED,
) -> tuple[list[str], list[list[str]], dict[str, int]]:
    metadata = pd.DataFrame(metadata_rows, columns=metadata_header)
    if "cell_id" not in metadata.columns:
        raise ValueError("Metadata must contain a `cell_id` column.")
    if stratify_key not in metadata.columns:
        raise ValueError(f"Metadata is missing stratification column: {stratify_key}")

    metadata["cell_id"] = metadata["cell_id"].astype(str)
    grouped = list(metadata.groupby(stratify_key, sort=False))
    if not grouped:
        raise ValueError("No strata available for shared-cell sampling.")
    if total_cells < len(grouped):
        raise ValueError(
            f"Requested total_cells={total_cells} is smaller than the number of strata={len(grouped)}."
        )

    base = total_cells // len(grouped)
    remainder = total_cells % len(grouped)
    rng = np.random.default_rng(random_seed)

    sampled_frames: list[pd.DataFrame] = []
    per_stratum_counts: dict[str, int] = {}
    for idx, (label, frame) in enumerate(grouped):
        target_n = base + (1 if idx < remainder else 0)
        if len(frame) < target_n:
            raise ValueError(
                f"Stratum {label} has only {len(frame)} cells, fewer than requested {target_n}."
            )
        sampled_positions = np.sort(
            rng.choice(len(frame), size=target_n, replace=False)
        )
        sampled_frame = frame.iloc[sampled_positions].copy()
        sampled_frames.append(sampled_frame)
        per_stratum_counts[str(label)] = int(target_n)

    sampled = pd.concat(sampled_frames, axis=0)
    sampled_cell_ids = sampled["cell_id"].astype(str).tolist()
    sampled_rows = sampled.astype(str).values.tolist()
    return sampled_cell_ids, sampled_rows, per_stratum_counts


def extract_member_if_missing(tar_path: Path, member_name: str, output_path: Path) -> None:
    if output_path.is_file():
        return
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r") as archive:
        member = archive.getmember(member_name)
        source = archive.extractfile(member)
        if source is None:
            raise FileNotFoundError(f"Could not extract {member_name} from {tar_path}")
        with output_path.open("wb") as target:
            shutil.copyfileobj(source, target, length=1024 * 1024)


def extract_3p_inputs(tar_path: Path, work_dir: Path) -> dict[str, Path]:
    work_dir.mkdir(parents=True, exist_ok=True)
    targets = {
        RNA_3P_BARCODES: work_dir / RNA_3P_BARCODES,
        RNA_3P_FEATURES: work_dir / RNA_3P_FEATURES,
        RNA_3P_MATRIX: work_dir / RNA_3P_MATRIX,
        ADT_3P_BARCODES: work_dir / ADT_3P_BARCODES,
        ADT_3P_FEATURES: work_dir / ADT_3P_FEATURES,
        ADT_3P_MATRIX: work_dir / ADT_3P_MATRIX,
    }
    for member_name, output_path in targets.items():
        extract_member_if_missing(tar_path, member_name, output_path)
    return targets


def read_matrix_market_shape(matrix_path: Path) -> tuple[int, int, int]:
    with gzip.open(matrix_path, "rt", newline="") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("%"):
                continue
            n_features, n_cells, n_entries = map(int, line.split())
            return n_features, n_cells, n_entries
    raise ValueError(f"Could not parse MatrixMarket header: {matrix_path}")


def iter_matrix_chunks(matrix_path: Path, chunk_rows: int = CHUNK_ROWS):
    reader = pd.read_csv(
        matrix_path,
        sep=r"\s+",
        compression="gzip",
        comment="%",
        header=None,
        names=["feature_idx", "cell_idx", "value"],
        skiprows=1,
        chunksize=chunk_rows,
        dtype={"feature_idx": np.int32, "cell_idx": np.int32, "value": np.float32},
        engine="c",
    )
    for chunk in reader:
        chunk["feature_idx"] = chunk["feature_idx"].astype(np.int64) - 1
        chunk["cell_idx"] = chunk["cell_idx"].astype(np.int64) - 1
        yield chunk


def compute_rna_totals_and_variance(
    matrix_path: Path,
    n_features: int,
    n_cells: int,
) -> tuple[np.ndarray, np.ndarray]:
    log("RNA pass 1/2: computing cell library sizes")
    cell_totals = np.zeros(n_cells, dtype=np.float64)
    for chunk in iter_matrix_chunks(matrix_path):
        np.add.at(cell_totals, chunk["cell_idx"].to_numpy(), chunk["value"].to_numpy(dtype=np.float64))
    cell_totals[cell_totals == 0.0] = 1.0

    log("RNA pass 2/2: computing normalized log1p variance")
    sums = np.zeros(n_features, dtype=np.float64)
    sums_sq = np.zeros(n_features, dtype=np.float64)
    for chunk in iter_matrix_chunks(matrix_path):
        feature_idx = chunk["feature_idx"].to_numpy()
        cell_idx = chunk["cell_idx"].to_numpy()
        values = chunk["value"].to_numpy(dtype=np.float64)
        transformed = np.log1p((values / cell_totals[cell_idx]) * RNA_TARGET_SUM)
        np.add.at(sums, feature_idx, transformed)
        np.add.at(sums_sq, feature_idx, transformed * transformed)
    means = sums / float(n_cells)
    variances = (sums_sq / float(n_cells)) - (means * means)
    return cell_totals, variances


def compute_adt_denominators_and_variance(
    matrix_path: Path,
    n_features: int,
    n_cells: int,
) -> tuple[np.ndarray, np.ndarray]:
    log("ADT pass 1/2: computing CLR denominators")
    log_sums = np.zeros(n_cells, dtype=np.float64)
    for chunk in iter_matrix_chunks(matrix_path):
        np.add.at(
            log_sums,
            chunk["cell_idx"].to_numpy(),
            np.log1p(chunk["value"].to_numpy(dtype=np.float64)),
        )
    denominators = np.exp(log_sums / float(n_features))

    log("ADT pass 2/2: computing CLR-like variance")
    sums = np.zeros(n_features, dtype=np.float64)
    sums_sq = np.zeros(n_features, dtype=np.float64)
    for chunk in iter_matrix_chunks(matrix_path):
        feature_idx = chunk["feature_idx"].to_numpy()
        cell_idx = chunk["cell_idx"].to_numpy()
        values = chunk["value"].to_numpy(dtype=np.float64)
        transformed = np.log1p(values / denominators[cell_idx])
        np.add.at(sums, feature_idx, transformed)
        np.add.at(sums_sq, feature_idx, transformed * transformed)
    means = sums / float(n_cells)
    variances = (sums_sq / float(n_cells)) - (means * means)
    return denominators, variances


def top_indices_by_variance(variances: np.ndarray, top_n: int) -> np.ndarray:
    top_n = min(top_n, variances.shape[0])
    return np.argsort(-variances, kind="stable")[:top_n].astype(np.int64)


def build_dense_subset(
    matrix_path: Path,
    n_cells: int,
    selected_indices: np.ndarray,
) -> np.ndarray:
    selected_map = {int(idx): col for col, idx in enumerate(selected_indices.tolist())}
    dense = np.zeros((n_cells, len(selected_indices)), dtype=np.float32)
    for chunk in iter_matrix_chunks(matrix_path):
        feature_idx = chunk["feature_idx"].to_numpy()
        cell_idx = chunk["cell_idx"].to_numpy()
        values = chunk["value"].to_numpy(dtype=np.float32)
        mask = np.isin(feature_idx, selected_indices, assume_unique=False)
        if not mask.any():
            continue
        for feat, cell, value in zip(feature_idx[mask], cell_idx[mask], values[mask], strict=False):
            dense[cell, selected_map[int(feat)]] = value
    return dense


def transform_rna(raw_counts: np.ndarray, cell_totals: np.ndarray) -> np.ndarray:
    matrix = raw_counts.astype(np.float64, copy=True)
    matrix /= cell_totals[:, None]
    matrix *= RNA_TARGET_SUM
    np.log1p(matrix, out=matrix)
    return matrix.astype(np.float32)


def transform_adt(raw_counts: np.ndarray, denominators: np.ndarray) -> np.ndarray:
    matrix = raw_counts.astype(np.float64, copy=True)
    matrix /= denominators[:, None]
    np.log1p(matrix, out=matrix)
    return matrix.astype(np.float32)


def main() -> None:
    args = parse_args()
    data_dir = args.data_dir.resolve()
    output_dir = (args.output_dir or (data_dir / "cscn_inputs")).resolve()
    work_dir = (args.work_dir or (data_dir / ".gse164378_3p_work")).resolve()
    raw_tar_path = data_dir / RAW_TAR_NAME
    metadata_path = data_dir / METADATA_3P_NAME

    validate_required_file(raw_tar_path, "raw GEO tar")
    validate_required_file(metadata_path, "3P metadata")

    log(f"extracting 3P members to {work_dir}")
    extracted = extract_3p_inputs(raw_tar_path, work_dir)

    log("loading barcodes, features, and metadata")
    rna_barcodes = read_gzip_lines(extracted[RNA_3P_BARCODES])
    adt_barcodes = read_gzip_lines(extracted[ADT_3P_BARCODES])
    rna_features = read_feature_names(extracted[RNA_3P_FEATURES])
    adt_features = read_feature_names(extracted[ADT_3P_FEATURES])
    metadata_header, metadata_cell_ids, metadata_rows = read_metadata_rows(metadata_path)

    if rna_barcodes != adt_barcodes:
        raise ValueError("RNA and ADT barcode orders do not match for GSE164378 3P.")
    if rna_barcodes != metadata_cell_ids:
        raise ValueError("3P metadata order does not match 3P barcodes.")

    rna_shape = read_matrix_market_shape(extracted[RNA_3P_MATRIX])
    adt_shape = read_matrix_market_shape(extracted[ADT_3P_MATRIX])
    if rna_shape[0] != len(rna_features) or rna_shape[1] != len(rna_barcodes):
        raise ValueError("RNA matrix shape does not match RNA features/barcodes.")
    if adt_shape[0] != len(adt_features) or adt_shape[1] != len(adt_barcodes):
        raise ValueError("ADT matrix shape does not match ADT features/barcodes.")

    rna_cell_totals, rna_variances = compute_rna_totals_and_variance(
        extracted[RNA_3P_MATRIX],
        n_features=rna_shape[0],
        n_cells=rna_shape[1],
    )
    adt_denominators, adt_variances = compute_adt_denominators_and_variance(
        extracted[ADT_3P_MATRIX],
        n_features=adt_shape[0],
        n_cells=adt_shape[1],
    )

    rna_top150 = top_indices_by_variance(rna_variances, RNA_ONLY_TOP_N)
    adt_top150 = top_indices_by_variance(adt_variances, ADT_ONLY_TOP_N)
    joint_rna_indices = rna_top150[: min(JOINT_RNA_TOP_N, len(rna_top150))]
    joint_adt_indices = adt_top150[: min(JOINT_ADT_TOP_N, len(adt_top150))]

    log("materializing selected dense subsets")
    rna_top150_counts = build_dense_subset(
        extracted[RNA_3P_MATRIX],
        n_cells=rna_shape[1],
        selected_indices=rna_top150,
    )
    adt_top150_counts = build_dense_subset(
        extracted[ADT_3P_MATRIX],
        n_cells=adt_shape[1],
        selected_indices=adt_top150,
    )

    log("applying modality-specific transforms")
    rna_top150_expr = transform_rna(rna_top150_counts, rna_cell_totals)
    adt_top150_expr = transform_adt(adt_top150_counts, adt_denominators)

    rna_top150_names = [rna_features[index] for index in rna_top150.tolist()]
    adt_top150_names = [adt_features[index] for index in adt_top150.tolist()]

    joint_rna_count = len(joint_rna_indices)
    joint_adt_count = len(joint_adt_indices)
    joint_expr = np.concatenate(
        [
            rna_top150_expr[:, :joint_rna_count],
            adt_top150_expr[:, :joint_adt_count],
        ],
        axis=1,
    ).astype(np.float32)
    joint_feature_names = (
        rna_top150_names[:joint_rna_count]
        + [f"ADT_{name}" for name in adt_top150_names[:joint_adt_count]]
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    metadata_output = output_dir / "gse164378_3p_metadata.csv.gz"
    rna_output = output_dir / "gse164378_3p_rna_only_expression.csv.gz"
    adt_output = output_dir / "gse164378_3p_adt_only_expression.csv.gz"
    joint_output = output_dir / "gse164378_3p_rna_adt_joint_expression.csv.gz"

    log(f"writing outputs to {output_dir}")
    write_gzip_csv(metadata_output, metadata_header, metadata_rows)
    write_expression_table(rna_output, rna_barcodes, rna_top150_names, rna_top150_expr)
    write_expression_table(adt_output, adt_barcodes, adt_top150_names, adt_top150_expr)
    write_expression_table(joint_output, rna_barcodes, joint_feature_names, joint_expr)

    write_text_lines(output_dir / "gse164378_3p_rna_only_features.txt", rna_top150_names)
    write_text_lines(output_dir / "gse164378_3p_adt_only_features.txt", adt_top150_names)
    write_text_lines(
        output_dir / "gse164378_3p_joint_rna_features.txt",
        rna_top150_names[:joint_rna_count],
    )
    write_text_lines(
        output_dir / "gse164378_3p_joint_adt_features.txt",
        adt_top150_names[:joint_adt_count],
    )

    log("building fixed shared-cell subset for fair modality comparison")
    shared_cell_ids, shared_metadata_rows, shared_strata_counts = sample_shared_cells(
        metadata_header,
        metadata_rows,
    )
    shared_index = {cell_id: idx for idx, cell_id in enumerate(rna_barcodes)}
    shared_positions = np.asarray(
        [shared_index[cell_id] for cell_id in shared_cell_ids],
        dtype=np.int64,
    )
    shared_prefix = f"gse164378_3p_shared{SHARED_TOTAL_CELLS}"

    shared_metadata_output = output_dir / f"{shared_prefix}_metadata.csv.gz"
    shared_cells_output = output_dir / f"{shared_prefix}_cells.txt"
    shared_rna_output = output_dir / f"{shared_prefix}_rna_only_expression.csv.gz"
    shared_adt_output = output_dir / f"{shared_prefix}_adt_only_expression.csv.gz"
    shared_joint_output = output_dir / f"{shared_prefix}_rna_adt_joint_expression.csv.gz"

    write_gzip_csv(shared_metadata_output, metadata_header, shared_metadata_rows)
    write_text_lines(shared_cells_output, shared_cell_ids)
    write_expression_table(
        shared_rna_output,
        shared_cell_ids,
        rna_top150_names,
        rna_top150_expr[shared_positions],
    )
    write_expression_table(
        shared_adt_output,
        shared_cell_ids,
        adt_top150_names,
        adt_top150_expr[shared_positions],
    )
    write_expression_table(
        shared_joint_output,
        shared_cell_ids,
        joint_feature_names,
        joint_expr[shared_positions],
    )
    write_text_lines(
        output_dir / f"{shared_prefix}_rna_only_features.txt",
        rna_top150_names,
    )
    write_text_lines(
        output_dir / f"{shared_prefix}_adt_only_features.txt",
        adt_top150_names,
    )
    write_text_lines(
        output_dir / f"{shared_prefix}_joint_rna_features.txt",
        rna_top150_names[:joint_rna_count],
    )
    write_text_lines(
        output_dir / f"{shared_prefix}_joint_adt_features.txt",
        adt_top150_names[:joint_adt_count],
    )

    summary = {
        "dataset": DATA_SET,
        "mode": "3P",
        "n_cells": len(rna_barcodes),
        "raw_rna_features": len(rna_features),
        "raw_adt_features": len(adt_features),
        "rna_only_features": len(rna_top150_names),
        "adt_only_features": len(adt_top150_names),
        "joint_rna_features": joint_rna_count,
        "joint_adt_features": joint_adt_count,
        "joint_total_features": int(joint_expr.shape[1]),
        "work_dir": str(work_dir),
        "outputs": {
            "metadata": str(metadata_output),
            "rna_only_expression": str(rna_output),
            "adt_only_expression": str(adt_output),
            "joint_expression": str(joint_output),
        },
        "shared_subset": {
            "total_cells": len(shared_cell_ids),
            "stratify_key": SHARED_STRATIFY_KEY,
            "random_seed": SHARED_RANDOM_SEED,
            "per_stratum_counts": shared_strata_counts,
            "outputs": {
                "metadata": str(shared_metadata_output),
                "cells": str(shared_cells_output),
                "rna_only_expression": str(shared_rna_output),
                "adt_only_expression": str(shared_adt_output),
                "joint_expression": str(shared_joint_output),
            },
        },
    }
    (output_dir / "gse164378_3p_cscn_input_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    log("done")


if __name__ == "__main__":
    main()
