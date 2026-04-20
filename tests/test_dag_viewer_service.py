from __future__ import annotations

import csv
import pickle
import sys
from pathlib import Path
from unittest import mock

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from biomarker.dag_viewer.service import (
    DagViewerError,
    discover_root_options,
    load_graph,
    resolve_viewer_root,
    scan_root,
)
from tests.support.fake_graph import FakeGraph


def test_resolve_viewer_root_accepts_dataset_and_dag_roots(tmp_path):
    dataset_root = build_dataset_root(tmp_path / "dataset")
    dataset_spec = resolve_viewer_root(dataset_root)
    dag_spec = resolve_viewer_root(dataset_root / "DAG")

    assert dataset_spec.root_kind == "dataset"
    assert dataset_spec.dataset_root == dataset_root
    assert dag_spec.root_kind == "dag"
    assert dag_spec.dataset_root == dataset_root
    assert dag_spec.dag_root == dataset_root / "DAG"


def test_scan_root_reports_group_metadata_and_result_ids(tmp_path):
    dataset_root = build_dataset_root(tmp_path / "dataset")

    payload = scan_root(dataset_root)

    assert payload["rootKind"] == "dataset"
    assert payload["methods"] == ["kTotal"]
    assert payload["collectionKind"] == "method"
    assert payload["collectionLabel"] == "方法"
    assert payload["showCollectionSelector"] is False
    assert len(payload["groups"]) == 1
    group = payload["groups"][0]
    assert group["method"] == "kTotal"
    assert group["groupKey"] == "Day54_cortical_interneuron"
    assert group["time"] == "54"
    assert group["cellType"] == "cortical interneuron"
    assert group["numDags"] == 2
    assert group["hasNodeMap"] is True
    assert group["hasConsensusCsv"] is True
    assert group["resultIds"] == [0, 1]


def test_scan_root_uses_single_collection_label_for_run_based_dataset(tmp_path):
    dataset_root = build_dataset_root(
        tmp_path / "dataset",
        collection_key="GSE121893_dhf_vs_n_all_all",
        group_key="dHF",
        with_used_genes=True,
        with_metadata=False,
        with_sampled_cells=True,
    )

    payload = scan_root(dataset_root)

    assert payload["methods"] == ["GSE121893_dhf_vs_n_all_all"]
    assert payload["collectionKind"] == "single"
    assert payload["collectionLabel"] == "当前视图"
    assert payload["showCollectionSelector"] is False
    assert payload["collections"][0]["displayLabel"] == "默认视图"


def test_load_graph_single_maps_nodes_and_resolves_cell_run_id(tmp_path):
    dataset_root = build_dataset_root(tmp_path / "dataset")

    payload = load_graph(
        root_path=dataset_root,
        method="kTotal",
        group_key="Day54_cortical_interneuron",
        mode="single",
        result_id=0,
    )

    node_by_raw = {node["rawId"]: node for node in payload["nodes"]}
    assert payload["mode"] == "single"
    assert payload["metadata"]["cellRunId"] == "SRR0001"
    assert node_by_raw["0"]["label"] == "CD99"
    assert node_by_raw["0"]["mapped"] is True
    assert node_by_raw["1"]["label"] == "1"
    assert node_by_raw["1"]["mapped"] is False


def test_load_graph_single_uses_used_genes_and_sampled_cells_fallbacks(tmp_path):
    dataset_root = build_dataset_root(
        tmp_path / "dataset",
        collection_key="GSE121893_dhf_vs_n_all_all",
        group_key="dHF",
        with_metadata=False,
        with_consensus_csv=False,
        with_used_genes=True,
        with_sampled_cells=True,
        result_graphs=[FakeGraph(nodes=[0, 1], edges=[(0, 1)])],
    )

    payload = load_graph(
        root_path=dataset_root,
        method="GSE121893_dhf_vs_n_all_all",
        group_key="dHF",
        mode="single",
        result_id=0,
    )

    node_by_raw = {node["rawId"]: node for node in payload["nodes"]}
    assert payload["metadata"]["cellRunId"] == "CELL_A"
    assert node_by_raw["0"]["label"] == "GENE_A"
    assert node_by_raw["0"]["mapped"] is True
    assert node_by_raw["1"]["label"] == "GENE_B"
    assert node_by_raw["1"]["mapped"] is True


def test_load_graph_consensus_prefers_csv_fast_path(tmp_path):
    dataset_root = build_dataset_root(tmp_path / "dataset")

    payload = load_graph(
        root_path=dataset_root,
        method="kTotal",
        group_key="Day54_cortical_interneuron",
        mode="consensus",
    )

    assert payload["mode"] == "consensus"
    assert payload["metadata"]["threshold"] == 2
    assert payload["metadata"]["numDags"] == 2
    assert payload["edges"][0]["count"] == 2


def test_load_graph_consensus_falls_back_to_pickle_aggregation(tmp_path):
    dataset_root = build_dataset_root(
        tmp_path / "dataset",
        with_consensus_csv=False,
        result_graphs=[
            FakeGraph(nodes=[0, 1, 2], edges=[(0, 1), (1, 1)]),
            FakeGraph(nodes=[0, 1], edges=[(0, 1)]),
        ],
    )

    payload = load_graph(
        root_path=dataset_root / "DAG",
        method="kTotal",
        group_key="Day54_cortical_interneuron",
        mode="consensus",
    )

    assert payload["metadata"]["threshold"] == 2
    assert payload["metadata"]["numDags"] == 2
    assert payload["edges"] == [
        {
            "id": "0->1",
            "source": "0",
            "target": "1",
            "count": 2,
            "weight": 2,
        }
    ]


def test_load_graph_reports_missing_pgmpy_dependency(tmp_path):
    dataset_root = build_dataset_root(tmp_path / "dataset")
    error = ModuleNotFoundError("No module named 'pgmpy'")
    error.name = "pgmpy"

    with mock.patch("biomarker.dag_viewer.service.pickle.load", side_effect=error):
        try:
            load_graph(
                root_path=dataset_root,
                method="kTotal",
                group_key="Day54_cortical_interneuron",
                mode="single",
                result_id=0,
            )
        except DagViewerError as exc:
            assert "pgmpy" in str(exc)
        else:
            raise AssertionError("DagViewerError was not raised for missing pgmpy")


def test_load_graph_degrades_gracefully_without_metadata_mapping(tmp_path):
    dataset_root = build_dataset_root(tmp_path / "dataset", with_metadata=False)

    payload = load_graph(
        root_path=dataset_root,
        method="kTotal",
        group_key="Day54_cortical_interneuron",
        mode="single",
        result_id=1,
    )

    assert "cellRunId" not in payload["metadata"]


def test_discover_root_options_reports_compatible_and_incompatible_roots(tmp_path):
    build_dataset_root(tmp_path / "data" / "E-GEOD-93593")
    (tmp_path / "data" / "GSE138852" / "DAG").mkdir(parents=True)

    payload = discover_root_options(tmp_path / "data")

    compatible = [root for root in payload["roots"] if root["compatible"]]
    incompatible = [root for root in payload["roots"] if not root["compatible"]]
    assert compatible[0]["label"] == "E-GEOD-93593"
    assert payload["defaultRootPath"] == compatible[0]["path"]
    assert incompatible[0]["label"] == "GSE138852"


def build_dataset_root(
    root: Path,
    *,
    collection_key: str = "kTotal",
    group_key: str = "Day54_cortical_interneuron",
    with_consensus_csv: bool = True,
    with_metadata: bool = True,
    with_sampled_cells: bool = False,
    with_used_genes: bool = False,
    result_graphs: list[FakeGraph] | None = None,
) -> Path:
    dag_group = root / "DAG" / collection_key / group_key
    dag_group.mkdir(parents=True, exist_ok=True)
    (root / "gene").mkdir(parents=True, exist_ok=True)

    if not with_used_genes:
        write_csv(
            root / "gene" / f"{collection_key}_node_map.csv",
            ["node_index", "gene_id", "gene_name"],
            [
                {"node_index": "0", "gene_id": "ENSG0000", "gene_name": "CD99"},
                {"node_index": "1", "gene_id": "ENSG0001", "gene_name": ""},
            ],
        )

    graphs = result_graphs or [
        FakeGraph(nodes=[0, 1], edges=[(0, 1)]),
        FakeGraph(nodes=[0, 1], edges=[(0, 1)]),
    ]
    for index, graph in enumerate(graphs):
        with (dag_group / f"result_{index}.pkl").open("wb") as handle:
            pickle.dump(graph, handle)

    if with_metadata:
        write_expression_metadata(root, group_key=group_key)

    if with_sampled_cells:
        write_sampled_cells(root, collection_key=collection_key, group_key=group_key)

    if with_used_genes:
        write_used_genes(root, collection_key=collection_key)

    if with_consensus_csv:
        consensus_dir = root / "visualizations" / collection_key / "consensus_edges"
        consensus_dir.mkdir(parents=True, exist_ok=True)
        write_csv(
            consensus_dir / f"{group_key}.csv",
            ["time", "cell_type", "from", "to", "edge", "count", "threshold", "num_dags"],
            [
                {
                    "time": "54",
                    "cell_type": "cortical interneuron",
                    "from": "CD99",
                    "to": "1",
                    "edge": "CD99→1",
                    "count": "2",
                    "threshold": "2",
                    "num_dags": "2",
                }
            ],
        )

    return root


def write_expression_metadata(root: Path, group_key: str):
    with (root / "expression_matrix.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["", "geneA", "geneB"])
        writer.writerow(["SRR0001", "1", "2"])
        writer.writerow(["SRR0002", "3", "4"])

    with (root / "Example.sdrf.txt").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Comment [ENA_RUN]",
                "Factor Value[time]",
                "Factor Value[inferred cell type - ontology labels]",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "Comment [ENA_RUN]": "SRR0001",
                "Factor Value[time]": "54",
                "Factor Value[inferred cell type - ontology labels]": group_key.replace("Day54_", "").replace("_", " "),
            }
        )
        writer.writerow(
            {
                "Comment [ENA_RUN]": "SRR0002",
                "Factor Value[time]": "54",
                "Factor Value[inferred cell type - ontology labels]": group_key.replace("Day54_", "").replace("_", " "),
            }
        )


def write_sampled_cells(root: Path, collection_key: str, group_key: str):
    write_csv(
        root / "output_deseq" / f"{collection_key}_{group_key}_sampled_cells.csv",
        ["cell_id"],
        [{"cell_id": "CELL_A"}, {"cell_id": "CELL_B"}],
    )


def write_used_genes(root: Path, collection_key: str):
    write_csv(
        root / "output_deseq" / f"{collection_key}_top150_genes_used.csv",
        ["gene"],
        [{"gene": "GENE_A"}, {"gene": "GENE_B"}],
    )


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
