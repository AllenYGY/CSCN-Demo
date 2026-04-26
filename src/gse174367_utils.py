from __future__ import annotations

from pathlib import Path
import re

import h5py
import numpy as np
import pandas as pd


DATA_SET = "GSE174367"
DEFAULT_CASE_GROUP = "AD"
DEFAULT_CONTROL_GROUP = "Control"
REQUIRED_METADATA_COLUMNS = (
    "Barcode",
    "SampleID",
    "Diagnosis",
    "Batch",
    "Cell.Type",
    "cluster",
    "Age",
    "Sex",
    "PMI",
    "Tangle.Stage",
    "Plaque.Stage",
    "RIN",
)


def normalize_slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip()).strip("_").lower()


def sample_sort_key(sample_id: str) -> tuple[int, str]:
    match = re.search(r"(\d+)$", str(sample_id))
    if match:
        return int(match.group(1)), str(sample_id)
    return 10**9, str(sample_id)


def build_run_slug(
    cell_type: str,
    case_group: str = DEFAULT_CASE_GROUP,
    control_group: str = DEFAULT_CONTROL_GROUP,
) -> str:
    return (
        f"{normalize_slug(cell_type)}_"
        f"{normalize_slug(case_group)}_vs_{normalize_slug(control_group)}"
    )


def deduplicate_gene_names(names: list[str]) -> list[str]:
    counts: dict[str, int] = {}
    deduped: list[str] = []
    for name in names:
        base = str(name)
        seen = counts.get(base, 0)
        counts[base] = seen + 1
        if seen == 0:
            deduped.append(base)
        else:
            deduped.append(f"{base}__{seen}")
    return deduped


def _decode_string_dataset(dataset) -> np.ndarray:
    if hasattr(dataset, "asstr"):
        return np.asarray(dataset.asstr()[:], dtype=object)
    values = dataset[:]
    decoded = []
    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8"))
        else:
            decoded.append(str(value))
    return np.asarray(decoded, dtype=object)


def load_gse174367_metadata(metadata_path: Path) -> pd.DataFrame:
    metadata_df = pd.read_csv(metadata_path)
    missing = [column for column in REQUIRED_METADATA_COLUMNS if column not in metadata_df.columns]
    if missing:
        raise ValueError(
            f"Missing required metadata columns in {metadata_path}: {missing}"
        )
    metadata_df = metadata_df.copy()
    for column in ("Barcode", "SampleID", "Diagnosis", "Batch", "Cell.Type", "cluster", "Sex"):
        metadata_df[column] = metadata_df[column].astype(str)
    for column in ("Age", "PMI", "RIN"):
        metadata_df[column] = pd.to_numeric(metadata_df[column], errors="coerce")
    return metadata_df


def filter_metadata_for_cell_type(
    metadata_df: pd.DataFrame,
    *,
    cell_type: str,
    case_group: str = DEFAULT_CASE_GROUP,
    control_group: str = DEFAULT_CONTROL_GROUP,
    min_cells_per_sample: int | None = None,
) -> pd.DataFrame:
    selected = metadata_df.loc[
        (metadata_df["Cell.Type"] == str(cell_type))
        & (metadata_df["Diagnosis"].isin([control_group, case_group]))
    ].copy()
    if selected.empty:
        raise ValueError(
            f"No metadata rows found for cell type {cell_type!r} and groups "
            f"{control_group!r}/{case_group!r}."
        )
    if min_cells_per_sample is not None:
        sample_counts = selected["SampleID"].value_counts()
        keep_samples = sample_counts[sample_counts >= min_cells_per_sample].index.astype(str)
        selected = selected.loc[selected["SampleID"].isin(keep_samples)].copy()
        if selected.empty:
            raise ValueError(
                f"No cells remain for {cell_type!r} after applying "
                f"min_cells_per_sample={min_cells_per_sample}."
            )
    return selected


def build_sample_metadata(selected_metadata_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample_id, frame in selected_metadata_df.groupby("SampleID", sort=False):
        row = {"sample": str(sample_id), "n_cells": int(frame.shape[0])}
        for source, target in (
            ("Diagnosis", "diagnosis"),
            ("Batch", "batch"),
            ("Age", "age"),
            ("Sex", "sex"),
            ("PMI", "pmi"),
            ("Tangle.Stage", "tangle_stage"),
            ("Plaque.Stage", "plaque_stage"),
            ("RIN", "rin"),
            ("Cell.Type", "cell_type"),
        ):
            unique_values = pd.unique(frame[source])
            non_missing = [
                value
                for value in unique_values
                if not (pd.isna(value) or str(value) == "nan")
            ]
            if len(non_missing) > 1:
                raise ValueError(
                    f"Sample {sample_id} has inconsistent values for {source}: {non_missing}"
                )
            row[target] = non_missing[0] if non_missing else np.nan
        rows.append(row)
    sample_metadata_df = pd.DataFrame(rows)
    if sample_metadata_df.empty:
        raise ValueError("No sample-level metadata rows were produced.")
    sample_metadata_df["sample_sort"] = sample_metadata_df["sample"].map(sample_sort_key)
    sample_metadata_df = sample_metadata_df.sort_values(
        by=["diagnosis", "sample_sort", "sample"],
        ascending=[True, True, True],
    ).drop(columns=["sample_sort"])
    sample_metadata_df = sample_metadata_df.reset_index(drop=True)
    return sample_metadata_df


def summarize_group_cells(
    selected_metadata_df: pd.DataFrame,
    *,
    case_group: str = DEFAULT_CASE_GROUP,
    control_group: str = DEFAULT_CONTROL_GROUP,
) -> dict[str, list[str]]:
    grouped = {
        control_group: [],
        case_group: [],
    }
    for diagnosis, frame in selected_metadata_df.groupby("Diagnosis", sort=False):
        if diagnosis in grouped:
            grouped[diagnosis] = frame["Barcode"].astype(str).tolist()
    return grouped


def resolve_default_sample_size(
    cells_by_group: dict[str, list[str]],
    *,
    cap: int = 3000,
) -> int:
    available = [len(cells) for cells in cells_by_group.values()]
    if not available or min(available) <= 0:
        raise ValueError("At least one group has no available cells.")
    return min(cap, min(available))


def _load_h5_layout(h5_path: Path) -> tuple[list[str], list[str], np.ndarray]:
    with h5py.File(h5_path, "r") as handle:
        matrix_group = handle["matrix"]
        barcodes = _decode_string_dataset(matrix_group["barcodes"]).astype(str).tolist()
        feature_names = _decode_string_dataset(matrix_group["features"]["name"]).astype(str)
        if "feature_type" in matrix_group["features"]:
            feature_types = _decode_string_dataset(matrix_group["features"]["feature_type"]).astype(str)
            keep_mask = feature_types == "Gene Expression"
        else:
            keep_mask = np.ones(feature_names.shape[0], dtype=bool)
        kept_rows = np.flatnonzero(keep_mask)
        gene_names = deduplicate_gene_names(feature_names[keep_mask].tolist())
        row_to_gene_pos = np.full(feature_names.shape[0], -1, dtype=np.int64)
        row_to_gene_pos[kept_rows] = np.arange(len(gene_names), dtype=np.int64)
    return barcodes, gene_names, row_to_gene_pos


def build_pseudobulk_count_matrix(
    h5_path: Path,
    selected_metadata_df: pd.DataFrame,
    sample_metadata_df: pd.DataFrame,
) -> pd.DataFrame:
    barcodes, gene_names, row_to_gene_pos = _load_h5_layout(h5_path)
    barcode_to_col = {barcode: idx for idx, barcode in enumerate(barcodes)}
    sample_ids = sample_metadata_df["sample"].astype(str).tolist()
    sample_to_index = {sample_id: idx for idx, sample_id in enumerate(sample_ids)}

    selected_columns: list[tuple[int, int]] = []
    missing_barcodes = []
    for row in selected_metadata_df.itertuples(index=False):
        barcode = str(getattr(row, "Barcode"))
        sample_id = str(getattr(row, "SampleID"))
        col_idx = barcode_to_col.get(barcode)
        if col_idx is None:
            missing_barcodes.append(barcode)
            continue
        selected_columns.append((col_idx, sample_to_index[sample_id]))
    if missing_barcodes:
        raise ValueError(
            f"Missing {len(missing_barcodes)} selected barcodes from the 10x matrix; "
            f"first 5: {missing_barcodes[:5]}"
        )

    selected_columns.sort(key=lambda item: item[0])
    count_matrix = np.zeros((len(sample_ids), len(gene_names)), dtype=np.int64)

    with h5py.File(h5_path, "r") as handle:
        matrix_group = handle["matrix"]
        data_ds = matrix_group["data"]
        indices_ds = matrix_group["indices"]
        indptr = matrix_group["indptr"][:]

        for col_idx, sample_idx in selected_columns:
            start = int(indptr[col_idx])
            end = int(indptr[col_idx + 1])
            row_ids = indices_ds[start:end]
            values = data_ds[start:end].astype(np.int64, copy=False)
            mapped = row_to_gene_pos[row_ids]
            valid = mapped >= 0
            if np.any(valid):
                np.add.at(count_matrix[sample_idx], mapped[valid], values[valid])

    count_df = pd.DataFrame(count_matrix.T, columns=sample_ids)
    count_df.insert(0, "gene", gene_names)
    return count_df


def extract_sampled_expression_from_h5(
    h5_path: Path,
    sampled_cells: dict[str, list[str]],
    top_genes: list[str],
) -> tuple[dict[str, np.ndarray], list[str]]:
    barcodes, all_gene_names, row_to_gene_pos = _load_h5_layout(h5_path)
    barcode_to_col = {barcode: idx for idx, barcode in enumerate(barcodes)}
    gene_name_to_pos = {gene: idx for idx, gene in enumerate(all_gene_names)}
    used_genes = [gene for gene in top_genes if gene in gene_name_to_pos]
    if not used_genes:
        raise ValueError("None of the requested top genes were found in the 10x matrix.")
    gene_pos_to_target = np.full(len(all_gene_names), -1, dtype=np.int64)
    for target_pos, gene_name in enumerate(used_genes):
        gene_pos_to_target[gene_name_to_pos[gene_name]] = target_pos
    row_to_target_pos = np.full(row_to_gene_pos.shape[0], -1, dtype=np.int64)
    valid_rows = row_to_gene_pos >= 0
    row_to_target_pos[valid_rows] = gene_pos_to_target[row_to_gene_pos[valid_rows]]

    matrices = {
        group: np.zeros((len(cell_ids), len(used_genes)), dtype=np.float64)
        for group, cell_ids in sampled_cells.items()
    }
    selected_columns: list[tuple[int, str, int]] = []
    for group, cell_ids in sampled_cells.items():
        missing = []
        for row_idx, barcode in enumerate(cell_ids):
            col_idx = barcode_to_col.get(str(barcode))
            if col_idx is None:
                missing.append(str(barcode))
                continue
            selected_columns.append((col_idx, group, row_idx))
        if missing:
            raise ValueError(
                f"Missing sampled barcodes from the 10x matrix for group {group}: {missing[:5]}"
            )

    selected_columns.sort(key=lambda item: item[0])
    with h5py.File(h5_path, "r") as handle:
        matrix_group = handle["matrix"]
        data_ds = matrix_group["data"]
        indices_ds = matrix_group["indices"]
        indptr = matrix_group["indptr"][:]

        for col_idx, group, row_idx in selected_columns:
            start = int(indptr[col_idx])
            end = int(indptr[col_idx + 1])
            row_ids = indices_ds[start:end]
            values = data_ds[start:end].astype(np.float64, copy=False)
            mapped = row_to_target_pos[row_ids]
            valid = mapped >= 0
            if np.any(valid):
                matrices[group][row_idx, mapped[valid]] = values[valid]

    return matrices, used_genes
