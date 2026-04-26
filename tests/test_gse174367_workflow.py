from __future__ import annotations

from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from gse174367_utils import (
    build_pseudobulk_count_matrix,
    build_sample_metadata,
    extract_sampled_expression_from_h5,
    filter_metadata_for_cell_type,
)


def write_fake_10x_h5(path: Path) -> None:
    dense = np.asarray(
        [
            [1, 0, 2, 0],
            [0, 3, 0, 4],
            [5, 0, 6, 0],
        ],
        dtype=np.int32,
    )
    indptr = [0]
    indices = []
    data = []
    for col_idx in range(dense.shape[1]):
        rows = np.flatnonzero(dense[:, col_idx])
        indices.extend(rows.tolist())
        data.extend(dense[rows, col_idx].tolist())
        indptr.append(len(indices))

    with h5py.File(path, "w") as handle:
        matrix = handle.create_group("matrix")
        matrix.create_dataset("data", data=np.asarray(data, dtype=np.int32))
        matrix.create_dataset("indices", data=np.asarray(indices, dtype=np.int64))
        matrix.create_dataset("indptr", data=np.asarray(indptr, dtype=np.int64))
        matrix.create_dataset("shape", data=np.asarray(dense.shape, dtype=np.int32))
        matrix.create_dataset(
            "barcodes",
            data=np.asarray([b"cellA", b"cellB", b"cellC", b"cellD"]),
        )
        features = matrix.create_group("features")
        features.create_dataset(
            "name",
            data=np.asarray([b"G1", b"G1", b"G2"]),
        )
        features.create_dataset(
            "feature_type",
            data=np.asarray([b"Gene Expression", b"Gene Expression", b"Gene Expression"]),
        )
        features.create_dataset(
            "id",
            data=np.asarray([b"id1", b"id2", b"id3"]),
        )


def build_fake_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Barcode": "cellA",
                "SampleID": "Sample-1",
                "Diagnosis": "AD",
                "Batch": "1",
                "Cell.Type": "ODC",
                "cluster": "ODC1",
                "Age": 80,
                "Sex": "F",
                "PMI": 4.0,
                "Tangle.Stage": "Stage 6",
                "Plaque.Stage": "Stage C",
                "RIN": 8.0,
            },
            {
                "Barcode": "cellB",
                "SampleID": "Sample-2",
                "Diagnosis": "Control",
                "Batch": "1",
                "Cell.Type": "ODC",
                "cluster": "ODC1",
                "Age": 79,
                "Sex": "M",
                "PMI": 5.0,
                "Tangle.Stage": "Stage 1",
                "Plaque.Stage": "Stage A",
                "RIN": 7.5,
            },
            {
                "Barcode": "cellC",
                "SampleID": "Sample-1",
                "Diagnosis": "AD",
                "Batch": "1",
                "Cell.Type": "ODC",
                "cluster": "ODC2",
                "Age": 80,
                "Sex": "F",
                "PMI": 4.0,
                "Tangle.Stage": "Stage 6",
                "Plaque.Stage": "Stage C",
                "RIN": 8.0,
            },
            {
                "Barcode": "cellD",
                "SampleID": "Sample-2",
                "Diagnosis": "Control",
                "Batch": "1",
                "Cell.Type": "ODC",
                "cluster": "ODC2",
                "Age": 79,
                "Sex": "M",
                "PMI": 5.0,
                "Tangle.Stage": "Stage 1",
                "Plaque.Stage": "Stage A",
                "RIN": 7.5,
            },
        ]
    )


def test_build_pseudobulk_count_matrix_and_extract_expression(tmp_path):
    h5_path = tmp_path / "fake_10x.h5"
    write_fake_10x_h5(h5_path)
    metadata_df = build_fake_metadata()

    selected_df = filter_metadata_for_cell_type(
        metadata_df,
        cell_type="ODC",
        case_group="AD",
        control_group="Control",
        min_cells_per_sample=1,
    )
    sample_metadata_df = build_sample_metadata(selected_df)
    count_df = build_pseudobulk_count_matrix(h5_path, selected_df, sample_metadata_df)

    assert count_df.columns.tolist() == ["gene", "Sample-1", "Sample-2"]
    assert count_df["gene"].tolist() == ["G1", "G1__1", "G2"]
    assert count_df["Sample-1"].tolist() == [3, 0, 11]
    assert count_df["Sample-2"].tolist() == [0, 7, 0]

    sampled_cells = {
        "Control": ["cellB", "cellD"],
        "AD": ["cellA", "cellC"],
    }
    matrices, used_genes = extract_sampled_expression_from_h5(
        h5_path=h5_path,
        sampled_cells=sampled_cells,
        top_genes=["G2", "G1__1", "missing"],
    )

    assert used_genes == ["G2", "G1__1"]
    assert matrices["AD"].shape == (2, 2)
    assert matrices["Control"].shape == (2, 2)
    assert matrices["AD"].tolist() == [[5.0, 0.0], [6.0, 0.0]]
    assert matrices["Control"].tolist() == [[0.0, 3.0], [0.0, 4.0]]
