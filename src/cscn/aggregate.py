from __future__ import annotations

import csv
import math
import pickle
from collections import defaultdict
from pathlib import Path
from typing import Any

import networkx as nx
import pandas as pd

from biomarker.dag_viewer.service import _load_pickle_graph, _normalize_graph


def load_node_map(run_dir: Path) -> dict[int, str]:
    genes_path = run_dir / "inputs" / "genes.csv"
    if not genes_path.is_file():
        return {}
    frame = pd.read_csv(genes_path)
    mapping: dict[int, str] = {}
    for row in frame.to_dict(orient="records"):
        try:
            node_index = int(row["node_index"])
        except (KeyError, TypeError, ValueError):
            continue
        mapping[node_index] = str(row.get("gene_name") or row.get("gene") or node_index)
    return mapping


def load_group_dags(group_dir: Path) -> list[tuple[int, Any]]:
    dags: list[tuple[int, Any]] = []
    for path in sorted(group_dir.glob("result_*.pkl")):
        result_id = int(path.stem.split("_", 1)[1])
        dags.append((result_id, _load_pickle_graph(path)))
    return dags


def map_group_dags_to_genes(
    dags: list[tuple[int, Any]],
    id2gene: dict[int, str],
) -> list[tuple[int, nx.DiGraph]]:
    mapped_dags: list[tuple[int, nx.DiGraph]] = []
    for dag_id, dag in dags:
        mapped = nx.DiGraph()
        for node in dag.nodes():
            mapped.add_node(id2gene.get(node, str(node)))
        for source, target in dag.edges():
            mapped.add_edge(id2gene.get(source, str(source)), id2gene.get(target, str(target)))
        mapped_dags.append((dag_id, mapped))
    return mapped_dags


def build_group_graph(dags: list[tuple[int, nx.DiGraph]]) -> nx.DiGraph:
    global_graph = nx.DiGraph()
    for _, dag in dags:
        global_graph.add_edges_from(dag.edges())
    return global_graph


def resolve_consensus_threshold(mode: str | int, num_dags: int) -> int:
    if isinstance(mode, int):
        return max(1, mode)
    if str(mode).lower() == "auto":
        return max(2, int(math.ceil(0.05 * num_dags)))
    try:
        return max(1, int(mode))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Unsupported consensus threshold mode: {mode}") from exc


def build_consensus_edges(
    group_dir: Path,
    node_map: dict[int, str],
    threshold_mode: str | int = "auto",
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    dags = load_group_dags(group_dir)
    if not dags:
        raise FileNotFoundError(f"No DAG files found in {group_dir}")

    edge_counts: defaultdict[tuple[str, str], int] = defaultdict(int)
    num_dags = 0
    for _, dag in dags:
        normalized = _normalize_graph(dag, node_map)
        node_labels = {
            str(node["id"]): str(node.get("label") or node["id"])
            for node in normalized["nodes"]
        }
        num_dags += 1
        for edge in normalized["edges"]:
            source = node_labels.get(str(edge["source"]), str(edge["source"]))
            target = node_labels.get(str(edge["target"]), str(edge["target"]))
            if source == target:
                continue
            edge_counts[(source, target)] += 1

    threshold = resolve_consensus_threshold(threshold_mode, num_dags)
    edges = []
    for (source, target), count in sorted(edge_counts.items()):
        if count < threshold:
            continue
        edges.append(
            {
                "from": source,
                "to": target,
                "count": count,
                "threshold": threshold,
                "num_dags": num_dags,
            }
        )
    return edges, {"threshold": threshold, "num_dags": num_dags}


def write_consensus_csv(
    group_dir: Path,
    output_path: Path,
    node_map: dict[int, str],
    threshold_mode: str | int = "auto",
) -> dict[str, int]:
    edges, metadata = build_consensus_edges(
        group_dir=group_dir,
        node_map=node_map,
        threshold_mode=threshold_mode,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["from", "to", "count", "threshold", "num_dags"],
        )
        writer.writeheader()
        writer.writerows(edges)
    metadata["num_edges"] = len(edges)
    metadata["num_nodes"] = len({row["from"] for row in edges}.union(row["to"] for row in edges))
    return metadata
