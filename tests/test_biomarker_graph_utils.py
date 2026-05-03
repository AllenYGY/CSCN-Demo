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
    _resolve_node_size_map,
    _resolve_layout,
    _resolve_legend_path,
    build_biomarker_dag,
    filter_graph_nodes,
    get_global_graph,
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


def test_build_biomarker_dag_only_connects_selected_outcome_parents():
    graph = nx.DiGraph()
    graph.add_edge("C", "T")
    graph.add_edge("C", "G2")
    graph.add_edge("C", "DISEASE")

    dag, _ = build_biomarker_dag(
        global_graph=graph,
        highlighted_treatments=["T"],
        candidate_treatments=["T", "G2"],
        outcome_parent_nodes=["T"],
        outcome_node="DISEASE",
    )

    assert dag.has_edge("T", "DISEASE")
    assert not dag.has_edge("G2", "DISEASE")


def test_resolve_node_size_map_scales_degree():
    graph = nx.DiGraph()
    graph.add_edge("A", "B")
    graph.add_edge("A", "C")
    graph.add_edge("B", "C")

    size_map = _resolve_node_size_map(
        graph,
        size_by="degree",
        min_node_size=100,
        max_node_size=300,
    )

    assert size_map["A"] == 300
    assert size_map["B"] == 100
    assert size_map["C"] == 300


def test_resolve_legend_path_appends_legend_suffix():
    legend_path = _resolve_legend_path(Path("/tmp/example.png"))
    assert legend_path == Path("/tmp/example_legend.png")


def test_concentric_layout_places_outcome_center_and_treatments_inside_other_nodes():
    graph = nx.DiGraph()
    graph.add_nodes_from(["DISEASE", "B1", "B2", "O1", "O2"])

    pos = _resolve_layout(
        graph,
        layout="concentric",
        treatment_nodes=["B1", "B2"],
        outcome_nodes=["DISEASE"],
    )

    assert pos["DISEASE"] == (0.0, 0.0)
    assert abs((pos["B1"][0] ** 2 + pos["B1"][1] ** 2) - 1.4**2) < 1e-6
    assert abs((pos["O1"][0] ** 2 + pos["O1"][1] ** 2) - 3.0**2) < 1e-6


def test_concentric_layout_separates_first_and_second_hop_neighbors():
    graph = nx.DiGraph()
    graph.add_edge("B1", "N1")
    graph.add_edge("N1", "N2")
    graph.add_node("DISEASE")

    pos = _resolve_layout(
        graph,
        layout="concentric",
        treatment_nodes=["B1"],
        outcome_nodes=["DISEASE"],
    )

    assert abs((pos["N1"][0] ** 2 + pos["N1"][1] ** 2) - 3.0**2) < 1e-6
    assert abs((pos["N2"][0] ** 2 + pos["N2"][1] ** 2) - 4.6**2) < 1e-6


def test_concentric_layout_opens_new_ring_when_capacity_is_exceeded():
    graph = nx.DiGraph()
    graph.add_node("DISEASE")
    graph.add_nodes_from(["B1", "B2", "B3"])

    pos = _resolve_layout(
        graph,
        layout="concentric",
        treatment_nodes=["B1", "B2", "B3"],
        outcome_nodes=["DISEASE"],
        max_nodes_per_ring=2,
        ring_growth_factor=1.0,
    )

    radii = sorted(
        round((pos[node][0] ** 2 + pos[node][1] ** 2) ** 0.5, 3)
        for node in ["B1", "B2", "B3"]
    )
    assert radii.count(1.4) == 2
    assert radii.count(3.0) == 1


def test_global_scope_style_targets_can_be_selected_without_labeling_all_nodes():
    graph = nx.DiGraph()
    graph.add_nodes_from(["DISEASE", "B1", "O1"])

    pos = _resolve_layout(
        graph,
        layout="concentric",
        treatment_nodes=["B1"],
        outcome_nodes=["DISEASE"],
    )

    assert "DISEASE" in pos
    assert "B1" in pos
    assert "O1" in pos


def test_get_global_graph_preserves_explicit_isolated_nodes():
    dag = nx.DiGraph()
    dag.add_edge("G1", "G2")

    global_graph = get_global_graph([(0, dag)], all_nodes=["G1", "G2", "G3"])

    assert set(global_graph.nodes()) == {"G1", "G2", "G3"}
    assert ("G1", "G2") in global_graph.edges()
