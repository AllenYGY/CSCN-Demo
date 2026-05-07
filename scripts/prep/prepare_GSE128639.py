from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shutil
import tarfile
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_SET = "GSE128639"
RAW_TAR_NAME = "GSE128639_RAW.tar"
ADT_BARCODES_NAME = "GSE128639_MNC_ADT_Barcodes.csv.gz"
HTO_BARCODES_NAME = "GSE128639_MNC_HTO_Barcodes.csv.gz"

RNA_MEMBER = "GSM3681518_MNC_RNA_counts.tsv.gz"
ADT_MEMBER = "GSM3681519_MNC_ADT_counts.tsv.gz"
HTO_MEMBER = "GSM3681520_MNC_HTO_counts.tsv.gz"

RNA_ONLY_TOP_N = 150
ADT_ONLY_TOP_N = 150
JOINT_RNA_TOP_N = 100
JOINT_ADT_TOP_N = 50
RNA_TARGET_SUM = 1e6
SHARED_TOTAL_CELLS = 3000
SHARED_STRATIFY_KEY = "hto_best_label"
SHARED_RANDOM_SEED = 42


def log(message: str) -> None:
    print(f"[{DATA_SET}] {message}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare GSE128639 MNC CSCN inputs for RNA-only, ADT-only, and "
            "RNA+ADT joint runs."
        )
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=REPO_ROOT / "data" / DATA_SET,
        help="Dataset directory containing GSE128639_RAW.tar.",
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
        help="Scratch directory for extracted members. Defaults to <data-dir>/.gse128639_work.",
    )
    return parser.parse_args(argv)


def validate_required_file(path: Path, description: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {description}: {path}")
    log(f"{description}: {path}")


def normalize_cell_id(cell_id: str) -> str:
    return re.sub(r"\.(\d+)$", r"-\1", str(cell_id))


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


def extract_inputs(tar_path: Path, work_dir: Path) -> dict[str, Path]:
    work_dir.mkdir(parents=True, exist_ok=True)
    targets = {
        RNA_MEMBER: work_dir / RNA_MEMBER,
        ADT_MEMBER: work_dir / ADT_MEMBER,
        HTO_MEMBER: work_dir / HTO_MEMBER,
    }
    for member_name, output_path in targets.items():
        extract_member_if_missing(tar_path, member_name, output_path)
    return targets


def read_header_cell_ids(path: Path) -> list[str]:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
    return [normalize_cell_id(cell_id) for cell_id in header]


def iter_dense_rows(path: Path):
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        normalized_header = [normalize_cell_id(cell_id) for cell_id in header]
        yield normalized_header
        for row in reader:
            if not row:
                continue
            feature_name = str(row[0])
            values = np.asarray(row[1:], dtype=np.float64)
            yield feature_name, values


def compute_rna_totals_and_variance(path: Path) -> tuple[list[str], np.ndarray, list[str], np.ndarray]:
    iterator = iter_dense_rows(path)
    cell_ids = next(iterator)

    log("RNA pass 1/3: computing cell library sizes")
    cell_totals = np.zeros(len(cell_ids), dtype=np.float64)
    for _, values in iterator:
        cell_totals += values
    cell_totals[cell_totals == 0.0] = 1.0

    log("RNA pass 2/3: computing normalized log1p variance")
    feature_names: list[str] = []
    variances: list[float] = []
    iterator = iter_dense_rows(path)
    _ = next(iterator)
    for feature_name, values in iterator:
        transformed = np.log1p((values / cell_totals) * RNA_TARGET_SUM)
        feature_names.append(feature_name)
        variances.append(float(np.var(transformed, dtype=np.float64)))
    return cell_ids, cell_totals, feature_names, np.asarray(variances, dtype=np.float64)


def build_selected_rna_matrix(
    path: Path,
    cell_totals: np.ndarray,
    selected_indices: np.ndarray,
) -> tuple[list[str], np.ndarray]:
    selected_lookup = {int(idx): pos for pos, idx in enumerate(selected_indices.tolist())}
    selected_names = [""] * len(selected_indices)
    dense = np.zeros((len(selected_indices), len(cell_totals)), dtype=np.float32)

    log("RNA pass 3/3: materializing selected normalized rows")
    iterator = iter_dense_rows(path)
    _ = next(iterator)
    for row_index, payload in enumerate(iterator):
        feature_name, values = payload
        if row_index not in selected_lookup:
            continue
        target_pos = selected_lookup[row_index]
        transformed = np.log1p((values / cell_totals) * RNA_TARGET_SUM)
        selected_names[target_pos] = feature_name
        dense[target_pos, :] = transformed.astype(np.float32)
    return selected_names, dense


def load_dense_matrix(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    iterator = iter_dense_rows(path)
    cell_ids = next(iterator)
    feature_names: list[str] = []
    rows: list[np.ndarray] = []
    for feature_name, values in iterator:
        feature_names.append(feature_name)
        rows.append(values.astype(np.float32))
    if not rows:
        return cell_ids, feature_names, np.zeros((0, len(cell_ids)), dtype=np.float32)
    return cell_ids, feature_names, np.vstack(rows).astype(np.float32)


def compute_adt_transform(raw_counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if raw_counts.shape[0] == 0:
        return np.ones(raw_counts.shape[1], dtype=np.float64), raw_counts.astype(np.float32)
    denominators = np.exp(np.log1p(raw_counts.astype(np.float64)).sum(axis=0) / float(raw_counts.shape[0]))
    transformed = np.log1p(raw_counts.astype(np.float64) / denominators[None, :])
    return denominators, transformed.astype(np.float32)


def top_indices_by_variance(variances: np.ndarray, top_n: int) -> np.ndarray:
    top_n = min(top_n, variances.shape[0])
    if top_n <= 0:
        return np.zeros(0, dtype=np.int64)
    return np.argsort(-variances, kind="stable")[:top_n].astype(np.int64)


def read_hto_label_map(path: Path | None) -> dict[str, str]:
    if path is None or not path.is_file():
        return {}
    labels: dict[str, str] = {}
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            raw_name = str(row.get("Name") or "").strip()
            match = re.search(r"(\d+)", raw_name)
            if not match:
                continue
            labels[f"HumanHTO{match.group(1)}"] = raw_name
    return labels


def build_metadata_rows(
    cell_ids: list[str],
    hto_feature_names: list[str],
    hto_counts: np.ndarray,
    hto_label_map: dict[str, str],
) -> tuple[list[str], list[list[str]]]:
    header = [
        "cell_id",
        "hto_best_feature",
        "hto_best_label",
        "hto_best_count",
        "hto_second_feature",
        "hto_second_label",
        "hto_second_count",
        "hto_total_count",
        "hto_detected_features",
    ]
    if hto_counts.shape[0] == 0:
        return header, [[cell_id, "", "", "0", "", "", "0", "0", "0"] for cell_id in cell_ids]

    rows: list[list[str]] = []
    for col_idx, cell_id in enumerate(cell_ids):
        counts = hto_counts[:, col_idx].astype(np.float64, copy=False)
        ranked = np.argsort(-counts, kind="stable")
        best_idx = int(ranked[0])
        second_idx = int(ranked[1]) if ranked.shape[0] > 1 else best_idx
        best_feature = hto_feature_names[best_idx]
        second_feature = hto_feature_names[second_idx]
        rows.append(
            [
                cell_id,
                best_feature,
                hto_label_map.get(best_feature, best_feature),
                str(int(counts[best_idx])),
                second_feature,
                hto_label_map.get(second_feature, second_feature),
                str(int(counts[second_idx])),
                str(int(counts.sum())),
                str(int(np.count_nonzero(counts))),
            ]
        )
    return header, rows


def write_gzip_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def write_feature_matrix(
    path: Path,
    cell_ids: list[str],
    feature_names: list[str],
    matrix: np.ndarray,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["feature_name", *cell_ids])
        for row_idx, feature_name in enumerate(feature_names):
            writer.writerow([feature_name, *matrix[row_idx, :].astype(np.float32).tolist()])


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
    column_to_idx = {name: idx for idx, name in enumerate(metadata_header)}
    if "cell_id" not in column_to_idx:
        raise ValueError("Metadata must contain a `cell_id` column.")
    if stratify_key not in column_to_idx:
        raise ValueError(f"Metadata is missing stratification column: {stratify_key}")

    grouped: dict[str, list[list[str]]] = {}
    for row in metadata_rows:
        label = str(row[column_to_idx[stratify_key]])
        grouped.setdefault(label, []).append(row)

    groups = list(grouped.items())
    if not groups:
        raise ValueError("No strata available for shared-cell sampling.")
    if total_cells < len(groups):
        raise ValueError(
            f"Requested total_cells={total_cells} is smaller than the number of strata={len(groups)}."
        )

    total_available = sum(len(rows) for _, rows in groups)
    target_total = min(total_cells, total_available)
    base = target_total // len(groups)
    remainder = target_total % len(groups)
    rng = np.random.default_rng(random_seed)

    sampled_rows: list[list[str]] = []
    per_stratum_counts: dict[str, int] = {}
    leftovers: list[tuple[str, list[list[str]]]] = []
    for idx, (label, rows) in enumerate(groups):
        target_n = min(len(rows), base + (1 if idx < remainder else 0))
        picks = np.sort(rng.choice(len(rows), size=target_n, replace=False)) if target_n > 0 else np.zeros(0, dtype=int)
        chosen = [rows[pos] for pos in picks.tolist()]
        sampled_rows.extend(chosen)
        per_stratum_counts[str(label)] = int(target_n)
        picked_positions = set(picks.tolist())
        remaining_rows = [row for pos, row in enumerate(rows) if pos not in picked_positions]
        if remaining_rows:
            leftovers.append((str(label), remaining_rows))

    remaining_needed = target_total - len(sampled_rows)
    if remaining_needed > 0:
        flat_leftovers: list[tuple[str, list[str]]] = []
        for label, rows in leftovers:
            for row in rows:
                flat_leftovers.append((label, row))
        if remaining_needed > len(flat_leftovers):
            raise ValueError(
                f"Could not reach requested sample size {target_total}; only {len(sampled_rows)} cells sampled."
            )
        extra_positions = np.sort(rng.choice(len(flat_leftovers), size=remaining_needed, replace=False))
        for pos in extra_positions.tolist():
            label, row = flat_leftovers[pos]
            sampled_rows.append(row)
            per_stratum_counts[label] = int(per_stratum_counts.get(label, 0) + 1)

    cell_idx = column_to_idx["cell_id"]
    sampled_cell_ids = [str(row[cell_idx]) for row in sampled_rows]
    return sampled_cell_ids, sampled_rows, per_stratum_counts


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    data_dir = args.data_dir.resolve()
    output_dir = (args.output_dir or (data_dir / "cscn_inputs")).resolve()
    work_dir = (args.work_dir or (data_dir / ".gse128639_work")).resolve()

    raw_tar_path = data_dir / RAW_TAR_NAME
    adt_barcodes_path = data_dir / ADT_BARCODES_NAME
    hto_barcodes_path = data_dir / HTO_BARCODES_NAME

    validate_required_file(raw_tar_path, "raw GEO tar")
    extracted = extract_inputs(raw_tar_path, work_dir)

    log("reading RNA header and computing RNA feature ranking")
    rna_cell_ids, rna_totals, rna_feature_names, rna_variances = compute_rna_totals_and_variance(
        extracted[RNA_MEMBER]
    )

    log("loading ADT and HTO dense matrices")
    adt_cell_ids, adt_feature_names, adt_raw = load_dense_matrix(extracted[ADT_MEMBER])
    hto_cell_ids, hto_feature_names, hto_raw = load_dense_matrix(extracted[HTO_MEMBER])

    if rna_cell_ids != adt_cell_ids:
        raise ValueError("RNA and ADT cell ids do not match after normalization.")
    if rna_cell_ids != hto_cell_ids:
        raise ValueError("RNA and HTO cell ids do not match after normalization.")

    _, adt_expr = compute_adt_transform(adt_raw)
    adt_variances = np.var(adt_expr.astype(np.float64), axis=1, dtype=np.float64)

    rna_top = top_indices_by_variance(rna_variances, RNA_ONLY_TOP_N)
    adt_top = top_indices_by_variance(adt_variances, ADT_ONLY_TOP_N)
    joint_rna_indices = rna_top[: min(JOINT_RNA_TOP_N, len(rna_top))]
    joint_adt_indices = adt_top[: min(JOINT_ADT_TOP_N, len(adt_top))]

    rna_only_names, rna_only_expr = build_selected_rna_matrix(
        extracted[RNA_MEMBER],
        rna_totals,
        rna_top,
    )
    adt_only_names = [adt_feature_names[idx] for idx in adt_top.tolist()]
    adt_only_expr = adt_expr[adt_top, :]

    joint_rna_names, joint_rna_expr = build_selected_rna_matrix(
        extracted[RNA_MEMBER],
        rna_totals,
        joint_rna_indices,
    )
    joint_adt_names = [adt_feature_names[idx] for idx in joint_adt_indices.tolist()]
    joint_adt_expr = adt_expr[joint_adt_indices, :]
    joint_feature_names = joint_rna_names + [f"ADT_{name}" for name in joint_adt_names]
    joint_expr = np.vstack([joint_rna_expr, joint_adt_expr]).astype(np.float32)

    hto_label_map = read_hto_label_map(hto_barcodes_path if hto_barcodes_path.is_file() else None)
    metadata_header, metadata_rows = build_metadata_rows(
        rna_cell_ids,
        hto_feature_names,
        hto_raw,
        hto_label_map,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    metadata_output = output_dir / "gse128639_mnc_metadata.csv.gz"
    rna_output = output_dir / "gse128639_mnc_rna_only_expression.tsv.gz"
    adt_output = output_dir / "gse128639_mnc_adt_only_expression.tsv.gz"
    joint_output = output_dir / "gse128639_mnc_rna_adt_joint_expression.tsv.gz"

    log(f"writing outputs to {output_dir}")
    write_gzip_csv(metadata_output, metadata_header, metadata_rows)
    write_feature_matrix(rna_output, rna_cell_ids, rna_only_names, rna_only_expr)
    write_feature_matrix(adt_output, rna_cell_ids, adt_only_names, adt_only_expr)
    write_feature_matrix(joint_output, rna_cell_ids, joint_feature_names, joint_expr)

    write_text_lines(output_dir / "gse128639_mnc_rna_only_features.txt", rna_only_names)
    write_text_lines(output_dir / "gse128639_mnc_adt_only_features.txt", adt_only_names)
    write_text_lines(output_dir / "gse128639_mnc_joint_rna_features.txt", joint_rna_names)
    write_text_lines(output_dir / "gse128639_mnc_joint_adt_features.txt", joint_adt_names)

    log("building fixed shared-cell subset for fairer medium-scale runs")
    shared_cell_ids, shared_metadata_rows, shared_strata_counts = sample_shared_cells(
        metadata_header,
        metadata_rows,
    )
    shared_index = {cell_id: idx for idx, cell_id in enumerate(rna_cell_ids)}
    shared_positions = np.asarray([shared_index[cell_id] for cell_id in shared_cell_ids], dtype=np.int64)
    shared_prefix = f"gse128639_mnc_shared{SHARED_TOTAL_CELLS}"

    shared_metadata_output = output_dir / f"{shared_prefix}_metadata.csv.gz"
    shared_cells_output = output_dir / f"{shared_prefix}_cells.txt"
    shared_rna_output = output_dir / f"{shared_prefix}_rna_only_expression.tsv.gz"
    shared_adt_output = output_dir / f"{shared_prefix}_adt_only_expression.tsv.gz"
    shared_joint_output = output_dir / f"{shared_prefix}_rna_adt_joint_expression.tsv.gz"

    write_gzip_csv(shared_metadata_output, metadata_header, shared_metadata_rows)
    write_text_lines(shared_cells_output, shared_cell_ids)
    write_feature_matrix(shared_rna_output, shared_cell_ids, rna_only_names, rna_only_expr[:, shared_positions])
    write_feature_matrix(shared_adt_output, shared_cell_ids, adt_only_names, adt_only_expr[:, shared_positions])
    write_feature_matrix(
        shared_joint_output,
        shared_cell_ids,
        joint_feature_names,
        joint_expr[:, shared_positions],
    )
    write_text_lines(output_dir / f"{shared_prefix}_rna_only_features.txt", rna_only_names)
    write_text_lines(output_dir / f"{shared_prefix}_adt_only_features.txt", adt_only_names)
    write_text_lines(output_dir / f"{shared_prefix}_joint_rna_features.txt", joint_rna_names)
    write_text_lines(output_dir / f"{shared_prefix}_joint_adt_features.txt", joint_adt_names)

    summary = {
        "dataset": DATA_SET,
        "mode": "MNC",
        "n_cells": len(rna_cell_ids),
        "raw_rna_features": len(rna_feature_names),
        "raw_adt_features": len(adt_feature_names),
        "raw_hto_features": len(hto_feature_names),
        "rna_only_features": len(rna_only_names),
        "adt_only_features": len(adt_only_names),
        "joint_rna_features": len(joint_rna_names),
        "joint_adt_features": len(joint_adt_names),
        "joint_total_features": int(joint_expr.shape[0]),
        "cell_id_normalization": ".<lane> suffixes normalized to -<lane>",
        "work_dir": str(work_dir),
        "source_files": {
            "raw_tar": str(raw_tar_path),
            "adt_barcodes": str(adt_barcodes_path) if adt_barcodes_path.is_file() else None,
            "hto_barcodes": str(hto_barcodes_path) if hto_barcodes_path.is_file() else None,
        },
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
    (output_dir / "gse128639_mnc_cscn_input_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    log("done")


if __name__ == "__main__":
    main()
