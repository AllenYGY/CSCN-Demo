from __future__ import annotations

import csv
import gzip
import json
import sys
import tarfile
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from cscn.config import load_config
from cscn.workflow import prepare_run
from scripts.prep.prepare_GSE128639 import main, normalize_cell_id


def _write_tsv_gz(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def _write_csv_gz(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def _build_synthetic_dataset(data_dir: Path) -> None:
    data_dir.mkdir(parents=True, exist_ok=True)

    scratch = data_dir / "scratch"
    scratch.mkdir(exist_ok=True)

    rna_path = scratch / "GSM3681518_MNC_RNA_counts.tsv.gz"
    adt_path = scratch / "GSM3681519_MNC_ADT_counts.tsv.gz"
    hto_path = scratch / "GSM3681520_MNC_HTO_counts.tsv.gz"

    rna_header = ["a_cell1.1", "a_cell2.1", "b_cell3.1"]
    adt_header = ["a_cell1.1", "a_cell2.1", "b_cell3.1"]
    hto_header = ["a_cell1-1", "a_cell2-1", "b_cell3-1"]

    _write_tsv_gz(
        rna_path,
        rna_header,
        [
            ["GeneA", 100, 0, 0],
            ["GeneB", 0, 20, 20],
            ["GeneC", 5, 5, 5],
        ],
    )
    _write_tsv_gz(
        adt_path,
        adt_header,
        [
            ["CD3", 100, 0, 5],
            ["CD19", 0, 200, 5],
        ],
    )
    _write_tsv_gz(
        hto_path,
        hto_header,
        [
            ["HumanHTO1", 50, 0, 2],
            ["HumanHTO2", 1, 60, 3],
        ],
    )

    raw_tar = data_dir / "GSE128639_RAW.tar"
    with tarfile.open(raw_tar, "w") as archive:
        archive.add(rna_path, arcname=rna_path.name)
        archive.add(adt_path, arcname=adt_path.name)
        archive.add(hto_path, arcname=hto_path.name)

    _write_csv_gz(
        data_dir / "GSE128639_MNC_HTO_Barcodes.csv.gz",
        ["Name", "Catalogue #", "Clone ", "Barcode sequence"],
        [
            ["Hashtag 1", "x", "x", "AAAA"],
            ["Hashtag 2", "x", "x", "BBBB"],
        ],
    )
    _write_csv_gz(
        data_dir / "GSE128639_MNC_ADT_Barcodes.csv.gz",
        ["Name", "Catalogue #", "Clone", "Barcode Sequence"],
        [
            ["CD3", "x", "x", "CCCC"],
            ["CD19", "x", "x", "DDDD"],
        ],
    )


def test_normalize_cell_id_converts_dot_suffix_to_dash_suffix():
    assert normalize_cell_id("a_AAAC.1") == "a_AAAC-1"
    assert normalize_cell_id("plain-cell") == "plain-cell"


def test_prepare_gse128639_builds_joint_inputs_and_is_cscn_compatible(tmp_path):
    data_dir = tmp_path / "data" / "GSE128639"
    _build_synthetic_dataset(data_dir)

    output_dir = data_dir / "cscn_inputs"
    work_dir = data_dir / "work"
    main(["--data-dir", str(data_dir), "--output-dir", str(output_dir), "--work-dir", str(work_dir)])

    metadata_path = output_dir / "gse128639_mnc_metadata.csv.gz"
    joint_path = output_dir / "gse128639_mnc_rna_adt_joint_expression.tsv.gz"
    adt_path = output_dir / "gse128639_mnc_adt_only_expression.tsv.gz"
    summary_path = output_dir / "gse128639_mnc_cscn_input_summary.json"

    assert metadata_path.is_file()
    assert joint_path.is_file()
    assert adt_path.is_file()
    assert summary_path.is_file()

    metadata = pd.read_csv(metadata_path)
    assert metadata["cell_id"].tolist() == ["a_cell1-1", "a_cell2-1", "b_cell3-1"]
    assert metadata["hto_best_label"].tolist() == ["Hashtag 1", "Hashtag 2", "Hashtag 2"]

    joint = pd.read_csv(joint_path, sep="\t")
    assert joint.columns[0] == "feature_name"
    assert joint.columns[1:].tolist() == ["a_cell1-1", "a_cell2-1", "b_cell3-1"]
    assert "ADT_CD3" in joint["feature_name"].tolist()
    assert "ADT_CD19" in joint["feature_name"].tolist()

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["n_cells"] == 3
    assert summary["raw_rna_features"] == 3
    assert summary["raw_adt_features"] == 2
    assert summary["raw_hto_features"] == 2
    assert summary["shared_subset"]["total_cells"] == 3
    assert summary["shared_subset"]["stratify_key"] == "hto_best_label"

    shared_metadata_path = output_dir / "gse128639_mnc_shared3000_metadata.csv.gz"
    shared_joint_path = output_dir / "gse128639_mnc_shared3000_rna_adt_joint_expression.tsv.gz"
    assert shared_metadata_path.is_file()
    assert shared_joint_path.is_file()

    shared_metadata = pd.read_csv(shared_metadata_path)
    assert shared_metadata["cell_id"].tolist() == ["a_cell1-1", "a_cell2-1", "b_cell3-1"]

    config_path = tmp_path / "gse128639_joint.yaml"
    config_path.write_text(
        "\n".join(
            [
                "run_name: gse128639_test_joint",
                "input:",
                "  format: tables",
                f"  expr_path: {joint_path}",
                f"  metadata_path: {metadata_path}",
                "  expr_orientation: genes_by_cells",
                "  gene_key: feature_name",
                "  metadata_cell_id_column: cell_id",
                "preprocess:",
                "  normalize: false",
                "  log1p: false",
                "  gene_selection:",
                "    top_n: 10",
                "run:",
                f"  output_dir: {tmp_path / 'runs'}",
                "aggregate:",
                "  consensus: false",
                "biomarker:",
                "  enabled: false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    config = load_config(config_path)
    summary = prepare_run(config)
    assert summary.groups == {"all": 3}
