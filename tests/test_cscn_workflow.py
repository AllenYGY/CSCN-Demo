from __future__ import annotations

import csv
import pickle
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from cscn.config import load_config
from cscn.layout import RunLayout
from cscn.viewer import discover_run_roots, load_run_graph, scan_run_root
from cscn.workflow import aggregate_run, prepare_run, run_biomarker_workflow
from tests.support.fake_graph import FakeGraph


def test_prepare_run_writes_standard_layout_from_table_inputs(tmp_path):
    config_path = build_table_config(tmp_path, sample_per_group=None, top_n=2)

    config = load_config(config_path)
    summary = prepare_run(config)
    layout = RunLayout.from_config(config)

    assert summary.run_dir == layout.run_dir
    assert layout.genes_path.is_file()
    assert layout.group_manifest_path.is_file()
    assert layout.cell_metadata_path.is_file()
    assert layout.matrix_path("A").is_file()
    assert layout.matrix_path("B").is_file()

    genes = pd.read_csv(layout.genes_path)["gene_name"].tolist()
    assert genes == ["G2", "G3"]

    groups = pd.read_csv(layout.group_manifest_path)
    assert groups["group_key"].tolist() == ["A", "B"]
    assert groups["n_cells"].tolist() == [2, 2]


def test_aggregate_run_writes_consensus_csv_and_viewer_can_read_run_layout(tmp_path):
    config_path = build_table_config(tmp_path, sample_per_group=None, top_n=2)
    config = load_config(config_path)
    prepare_run(config)
    layout = RunLayout.from_config(config)
    create_fake_dags(layout, "A", [FakeGraph(nodes=[0, 1], edges=[(0, 1)]), FakeGraph(nodes=[0, 1], edges=[(0, 1)])])
    create_fake_dags(layout, "B", [FakeGraph(nodes=[0, 1], edges=[(1, 0)]), FakeGraph(nodes=[0, 1], edges=[(1, 0)])])

    aggregate_run(config)

    consensus_path = layout.consensus_csv_path("A")
    assert consensus_path.is_file()
    rows = list(csv.DictReader(consensus_path.open("r", encoding="utf-8", newline="")))
    assert rows[0]["from"] == "G2"
    assert rows[0]["to"] == "G3"
    assert rows[0]["threshold"] == "2"

    roots = discover_run_roots(layout.base_output_dir, default_root=layout.run_dir)
    assert roots["defaultRootPath"] == str(layout.run_dir)

    scan_payload = scan_run_root(layout.run_dir)
    assert scan_payload["rootKind"] == "run"
    assert scan_payload["collectionKind"] == "single"
    assert scan_payload["groups"][0]["method"] == layout.run_name

    single_payload = load_run_graph(
        root_path=layout.run_dir,
        method=layout.run_name,
        group_key="A",
        mode="single",
        result_id=0,
    )
    assert single_payload["metadata"]["cellRunId"] == "c1"
    assert {node["label"] for node in single_payload["nodes"]} == {"G2", "G3"}

    consensus_payload = load_run_graph(
        root_path=layout.run_dir,
        method=layout.run_name,
        group_key="A",
        mode="consensus",
    )
    assert consensus_payload["metadata"]["threshold"] == 2
    assert consensus_payload["edges"][0]["count"] == 2


def test_biomarker_requires_case_and_control_groups(tmp_path):
    config_path = build_table_config(tmp_path, sample_per_group=None, top_n=2)
    config = load_config(config_path)
    prepare_run(config)

    try:
        run_biomarker_workflow(config)
    except ValueError as exc:
        assert "case_group" in str(exc)
    else:
        raise AssertionError("Expected biomarker workflow to reject missing case/control groups")


def test_prepare_run_accepts_blank_metadata_cell_id_column(tmp_path):
    expr_path = tmp_path / "expression.csv"
    metadata_path = tmp_path / "metadata.csv"
    expr_path.write_text(
        "\n".join(
            [
                "cell_id,G1,G2",
                "c1,1,10",
                "c2,2,20",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    metadata_path.write_text(
        "\n".join(
            [
                ",condition",
                "c1,A",
                "c2,A",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "\n".join(
            [
                "run_name: blank_metadata_id",
                "input:",
                "  format: tables",
                "  expr_path: expression.csv",
                "  metadata_path: metadata.csv",
                "  expr_orientation: cells_by_genes",
                "  expr_cell_id_column: cell_id",
                '  metadata_cell_id_column: ""',
                "  obs_group_key: condition",
                "preprocess:",
                "  gene_selection:",
                "    top_n: 2",
                "run:",
                "  output_dir: runs",
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

    assert config.input.metadata_cell_id_column == ""
    assert summary.groups == {"A": 2}


def build_table_config(tmp_path: Path, *, sample_per_group: int | None, top_n: int) -> Path:
    expr_path = tmp_path / "expression.csv"
    metadata_path = tmp_path / "metadata.csv"
    expr_path.write_text(
        "\n".join(
            [
                "cell_id,G1,G2,G3",
                "c1,1,10,1",
                "c2,2,20,1",
                "c3,3,30,10",
                "c4,4,40,10",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    metadata_path.write_text(
        "\n".join(
            [
                "cell_id,condition",
                "c1,A",
                "c2,A",
                "c3,B",
                "c4,B",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    sample_line = "null" if sample_per_group is None else str(sample_per_group)
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "\n".join(
            [
                "run_name: test_run",
                "input:",
                "  format: tables",
                "  expr_path: expression.csv",
                "  metadata_path: metadata.csv",
                "  expr_orientation: cells_by_genes",
                "  expr_cell_id_column: cell_id",
                "  metadata_cell_id_column: cell_id",
                "  obs_group_key: condition",
                "preprocess:",
                f"  sample_per_group: {sample_line}",
                "  random_seed: 7",
                "  gene_selection:",
                f"    top_n: {top_n}",
                "run:",
                "  output_dir: runs",
                "aggregate:",
                "  consensus: true",
                "  consensus_threshold_mode: auto",
                "biomarker:",
                "  enabled: false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return config_path


def create_fake_dags(layout: RunLayout, group_key: str, graphs: list[FakeGraph]) -> None:
    dag_dir = layout.dag_group_dir(group_key)
    dag_dir.mkdir(parents=True, exist_ok=True)
    for index, graph in enumerate(graphs):
        with (dag_dir / f"result_{index}.pkl").open("wb") as handle:
            pickle.dump(graph, handle)
