from __future__ import annotations

import sys
from pathlib import Path

import networkx as nx

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from biomarker.graph_utils import (
    build_biomarker_dag,
    filter_graph_nodes,
    map_node_id_to_gene_directed,
)


def test_map_node_id_to_gene_directed_preserves_edge_direction():
    dag = nx.DiGraph()
    dag.add_edge(0, 1)
    mapped = map_node_id_to_gene_directed([(0, dag)], {0: "G0", 1: "G1"})

    _, mapped_graph = mapped[0]
    assert isinstance(mapped_graph, nx.DiGraph)
    assert mapped_graph.has_edge("G0", "G1")
    assert not mapped_graph.has_edge("G1", "G0")


def test_filter_graph_nodes_returns_same_graph_type():
    graph = nx.DiGraph()
    graph.add_edge("A", "B")
    graph.add_edge("B", "C")

    filtered = filter_graph_nodes(graph, ["A", "B"])

    assert isinstance(filtered, nx.DiGraph)
    assert set(filtered.nodes()) == {"A", "B"}
    assert list(filtered.edges()) == [("A", "B")]


def test_build_biomarker_dag_adds_outcome_and_confounders():
    graph = nx.DiGraph()
    graph.add_edge("C", "T")
    graph.add_edge("C", "DISEASE")
    graph.add_edge("T", "X")

    dag, confounders_by_treatment = build_biomarker_dag(
        global_graph=graph,
        highlighted_treatments=["T"],
        candidate_treatments=["T"],
        outcome_node="DISEASE",
    )

    assert set(dag.nodes()) == {"C", "T", "DISEASE"}
    assert dag.has_edge("C", "T")
    assert dag.has_edge("T", "DISEASE")
    assert confounders_by_treatment == {"T": ["C"]}
