import os
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import networkx as nx
import pandas as pd


def map_node_id_to_gene(dags, id2gene):
    new_dags = []
    for dag_id, dag in dags:
        new_dag = nx.Graph()
        try:
            node_mapping = {node: id2gene[node] for node in dag.nodes()}
        except KeyError as e:
            print(f"Warning: id {e} not found in id2gene table")
            continue
        for old_node, new_node in node_mapping.items():
            new_dag.add_node(new_node, **dag.nodes[old_node])
        for u, v in dag.edges():
            new_dag.add_edge(node_mapping[u], node_mapping[v], **dag.edges[u, v])
        new_dags.append((dag_id, new_dag))
    return new_dags


def map_node_id_to_gene_directed(dags, id2gene):
    new_dags = []
    for dag_id, dag in dags:
        new_dag = nx.DiGraph()
        try:
            node_mapping = {node: id2gene[node] for node in dag.nodes()}
        except KeyError as e:
            print(f"Warning: id {e} not found in id2gene table")
            continue
        for old_node, new_node in node_mapping.items():
            new_dag.add_node(new_node, **dag.nodes[old_node])
        for u, v in dag.edges():
            new_dag.add_edge(node_mapping[u], node_mapping[v], **dag.edges[u, v])
        new_dags.append((dag_id, new_dag))
    return new_dags


def get_global_graph(dags, all_nodes=None):
    global_graph = nx.DiGraph()
    if all_nodes is not None:
        global_graph.add_nodes_from(all_nodes)
    for _, dag in dags:
        global_graph.add_nodes_from(dag.nodes())
        global_graph.add_edges_from(dag.edges())
    return global_graph


def draw_global_network(global_graph, save_path="results", title="Global Network of DAGs"):
    os.makedirs(save_path, exist_ok=True)
    pos = nx.circular_layout(global_graph)
    plt.figure(figsize=(12, 10))
    nx.draw_networkx_nodes(
        global_graph,
        pos,
        node_color="#A0CBE2",
        node_size=700,
        edgecolors="k",
        linewidths=1.5,
        alpha=0.9,
    )
    nx.draw_networkx_edges(
        global_graph,
        pos,
        width=1.5,
        alpha=0.6,
        arrowsize=22,
        edge_color="#696969",
        connectionstyle="arc3,rad=0.05",
    )
    nx.draw_networkx_labels(
        global_graph,
        pos,
        font_size=12,
        font_color="black",
        font_family="sans-serif",
        font_weight="bold",
    )
    plt.title(f"{title}", fontsize=22, fontweight="bold", y=1.02)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(f"{save_path}/{title}.png", dpi=300)
    plt.show()


def filter_graph_nodes(graph, node_list):
    if not node_list:
        return graph.__class__()
    existing_nodes = [node for node in node_list if node in graph]
    if not existing_nodes:
        return graph.__class__()
    return graph.subgraph(existing_nodes).copy()


def build_biomarker_dag(
    global_graph,
    highlighted_treatments,
    outcome_node="DISEASE",
    candidate_treatments=None,
    outcome_parent_nodes=None,
    confounder_method="classic",
    connect_treatments_to_outcome=True,
):
    candidate_treatments = list(
        candidate_treatments if candidate_treatments is not None else highlighted_treatments
    )
    outcome_parent_nodes = list(
        outcome_parent_nodes if outcome_parent_nodes is not None else highlighted_treatments
    )
    graph_with_outcome = add_sink_node_to_graph(global_graph.copy(), sink_node_name=outcome_node)
    filter_nodes = set(highlighted_treatments)
    confounders_by_treatment = {}

    for treatment in candidate_treatments:
        if treatment not in graph_with_outcome:
            confounders_by_treatment[treatment] = []
            continue
        confounders = find_confounders(
            graph_with_outcome,
            treatment,
            outcome_node,
            method=confounder_method,
        )
        confounders_by_treatment[treatment] = confounders
        filter_nodes.update(confounders)

    filtered_graph = filter_graph_nodes(global_graph, sorted(filter_nodes))
    filtered_graph.add_node(outcome_node)
    if connect_treatments_to_outcome:
        for treatment in outcome_parent_nodes:
            if treatment in filtered_graph:
                filtered_graph.add_edge(treatment, outcome_node)
    return filtered_graph, confounders_by_treatment


def _resolve_highlighted_dag_path(save_path):
    output_path = Path(save_path)
    if output_path.suffix.lower() != ".png":
        output_path = output_path.with_suffix(".png")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    return output_path


def _resolve_legend_path(output_path):
    return output_path.with_name(f"{output_path.stem}_legend{output_path.suffix}")


def _resolve_concentric_layout(global_graph, treatment_nodes, outcome_nodes):
    import math

    outcome_set = set(outcome_nodes)
    treatment_set = set(treatment_nodes) - outcome_set
    other_nodes = [
        node for node in global_graph.nodes() if node not in outcome_set and node not in treatment_set
    ]

    pos = {}
    for node in outcome_set:
        pos[node] = (0.0, 0.0)

    if treatment_set:
        treatment_nodes_sorted = sorted(treatment_set)
        radius = 1.4
        for idx, node in enumerate(treatment_nodes_sorted):
            angle = 2.0 * math.pi * idx / len(treatment_nodes_sorted)
            pos[node] = (radius * math.cos(angle), radius * math.sin(angle))

    if other_nodes:
        other_nodes_sorted = sorted(other_nodes)
        radius = 3.0
        for idx, node in enumerate(other_nodes_sorted):
            angle = 2.0 * math.pi * idx / len(other_nodes_sorted)
            pos[node] = (radius * math.cos(angle), radius * math.sin(angle))

    return pos


def _resolve_layout(
    global_graph,
    layout="spring",
    layout_seed=42,
    treatment_nodes=None,
    outcome_nodes=None,
):
    if layout == "spring":
        return nx.spring_layout(global_graph, seed=layout_seed)
    if layout == "kamada_kawai":
        return nx.kamada_kawai_layout(global_graph)
    if layout == "circular":
        return nx.circular_layout(global_graph)
    if layout == "concentric":
        return _resolve_concentric_layout(
            global_graph,
            treatment_nodes=treatment_nodes or [],
            outcome_nodes=outcome_nodes or [],
        )
    if layout == "dag":
        layout_graph = global_graph.copy()
        if nx.is_directed_acyclic_graph(global_graph):
            for layer, nodes in enumerate(nx.topological_generations(global_graph)):
                for node in nodes:
                    layout_graph.nodes[node]["layer"] = layer
            return nx.multipartite_layout(layout_graph, subset_key="layer", align="horizontal")
        return nx.kamada_kawai_layout(global_graph)
    raise ValueError(
        f"Unsupported layout '{layout}'. Use one of: spring, kamada_kawai, circular, concentric, dag."
    )


def _resolve_node_size_map(
    global_graph,
    size_by=None,
    min_node_size=700,
    max_node_size=1400,
):
    if size_by is None:
        return {}

    if size_by == "degree":
        metric = dict(global_graph.degree())
    elif size_by == "in_degree":
        metric = dict(global_graph.in_degree())
    elif size_by == "out_degree":
        metric = dict(global_graph.out_degree())
    else:
        raise ValueError(
            f"Unsupported size_by '{size_by}'. Use one of: degree, in_degree, out_degree."
        )

    values = list(metric.values())
    if not values:
        return {}
    min_value = min(values)
    max_value = max(values)
    if min_value == max_value:
        midpoint = int((min_node_size + max_node_size) / 2)
        return {node: midpoint for node in global_graph.nodes()}

    scale = max_node_size - min_node_size
    return {
        node: int(min_node_size + ((value - min_value) / (max_value - min_value)) * scale)
        for node, value in metric.items()
    }


def _resolve_node_metric(global_graph, size_by=None):
    if size_by is None:
        return {}
    if size_by == "degree":
        return dict(global_graph.degree())
    if size_by == "in_degree":
        return dict(global_graph.in_degree())
    if size_by == "out_degree":
        return dict(global_graph.out_degree())
    raise ValueError(
        f"Unsupported size_by '{size_by}'. Use one of: degree, in_degree, out_degree."
    )


def _legend_marker_size(node_size):
    return max(8, (float(node_size) ** 0.5) / 2.0)


def draw_highlighted_dag_legend(
    save_path,
    size_by=None,
    min_node_size=700,
    max_node_size=1400,
    metric_values=None,
):
    output_path = _resolve_highlighted_dag_path(save_path)
    legend_path = _resolve_legend_path(output_path)

    handles = [
        Patch(facecolor="lightcoral", edgecolor="firebrick", label="Biomarker"),
        Patch(facecolor="skyblue", edgecolor="navy", label="Outcome"),
        Patch(facecolor="lightgray", edgecolor="darkgray", label="Other"),
    ]

    if size_by is not None and metric_values:
        unique_values = sorted(set(metric_values.values()))
        low_value = unique_values[0]
        high_value = unique_values[-1]
        mid_value = unique_values[len(unique_values) // 2]
        scale_values = [low_value, mid_value, high_value]
        min_metric = low_value
        max_metric = high_value

        if min_metric == max_metric:
            scale_sizes = [int((min_node_size + max_node_size) / 2)] * 3
        else:
            scale_sizes = [
                int(
                    min_node_size
                    + ((value - min_metric) / (max_metric - min_metric))
                    * (max_node_size - min_node_size)
                )
                for value in scale_values
            ]

        for value, size in zip(scale_values, scale_sizes):
            handles.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    label=f"{size_by} = {value}",
                    markerfacecolor="gray",
                    markeredgecolor="gray",
                    markersize=_legend_marker_size(size),
                )
            )

    fig, ax = plt.subplots(figsize=(4.5, 3.5))
    ax.axis("off")
    ax.legend(handles=handles, loc="center", frameon=False, title="Legend")
    fig.tight_layout()
    fig.savefig(legend_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return legend_path


def draw_global_network_highlighted(
    global_graph,
    treatment_nodes,
    outcome_nodes,
    save_path="results",
    title="Biomarker DAG",
    layout="spring",
    layout_seed=42,
    label_scope="all",
    size_by=None,
    min_node_size=700,
    max_node_size=1400,
):
    output_path = _resolve_highlighted_dag_path(save_path)
    pos = _resolve_layout(
        global_graph,
        layout=layout,
        layout_seed=layout_seed,
        treatment_nodes=treatment_nodes,
        outcome_nodes=outcome_nodes,
    )
    plt.figure(figsize=(14, 10))

    node_colors = []
    node_sizes = []
    edge_colors = []
    linewidths = []
    treatment_set = set(treatment_nodes)
    outcome_set = set(outcome_nodes)
    metric_values = _resolve_node_metric(global_graph, size_by=size_by)
    size_map = _resolve_node_size_map(
        global_graph,
        size_by=size_by,
        min_node_size=min_node_size,
        max_node_size=max_node_size,
    )

    for node in global_graph.nodes():
        base_size = size_map.get(node)
        if node in outcome_set:
            node_colors.append("skyblue")
            node_sizes.append(base_size if base_size is not None else 1100)
            edge_colors.append("navy")
            linewidths.append(2.2)
        elif node in treatment_set:
            node_colors.append("lightcoral")
            node_sizes.append(base_size if base_size is not None else 950)
            edge_colors.append("firebrick")
            linewidths.append(2.0)
        else:
            node_colors.append("lightgray")
            node_sizes.append(base_size if base_size is not None else 700)
            edge_colors.append("darkgray")
            linewidths.append(1.0)

    nx.draw_networkx_nodes(
        global_graph,
        pos,
        node_color=node_colors,
        node_size=node_sizes,
        edgecolors=edge_colors,
        linewidths=linewidths,
    )
    nx.draw_networkx_edges(
        global_graph,
        pos,
        width=1.8,
        alpha=0.7,
        arrows=True,
        arrowsize=30,
        arrowstyle="-|>",
        edge_color="dimgray",
        connectionstyle="arc3,rad=0.04",
    )
    label_nodes = {}
    if label_scope == "all":
        label_nodes = {node: node for node in global_graph.nodes()}
    elif label_scope == "biomarkers":
        visible = set(treatment_nodes) | set(outcome_nodes)
        label_nodes = {node: node for node in global_graph.nodes() if node in visible}
    elif label_scope == "none":
        label_nodes = {}
    else:
        raise ValueError(
            f"Unsupported label_scope '{label_scope}'. Use one of: all, biomarkers, none."
        )

    if label_nodes:
        nx.draw_networkx_labels(
            global_graph,
            pos,
            labels=label_nodes,
            font_size=10,
            font_color="dimgray",
            font_family="sans-serif",
        )

    plt.title(title, fontsize=20, fontweight="bold", color="dimgray")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    legend_path = draw_highlighted_dag_legend(
        save_path=str(output_path),
        size_by=size_by,
        min_node_size=min_node_size,
        max_node_size=max_node_size,
        metric_values=metric_values,
    )
    return output_path, legend_path


def add_sink_node_to_graph(graph, sink_node_name="SINK"):
    leaf_nodes = [node for node in graph.nodes()]
    graph.add_node(sink_node_name)
    for leaf_node in leaf_nodes:
        graph.add_edge(leaf_node, sink_node_name)

    print(
        f"Added sink node '{sink_node_name}' and connected {len(leaf_nodes)} leaf nodes to it."
    )
    return graph


def find_confounders(dag, treatment_node, outcome_node, method="classic"):
    if treatment_node not in dag.nodes() or outcome_node not in dag.nodes():
        raise ValueError(
            f"Treatment node '{treatment_node}' or outcome node '{outcome_node}' not in DAG"
        )

    if method == "classic":
        return _find_classic_confounders(dag, treatment_node, outcome_node)
    if method == "backdoor":
        return _find_backdoor_confounders(dag, treatment_node, outcome_node)
    if method == "parents":
        return _find_parent_confounders(dag, treatment_node, outcome_node)
    if method == "all":
        return {
            "classic": _find_classic_confounders(dag, treatment_node, outcome_node),
            "backdoor": _find_backdoor_confounders(dag, treatment_node, outcome_node),
            "parents": _find_parent_confounders(dag, treatment_node, outcome_node),
        }
    raise ValueError("Method must be 'classic', 'backdoor', 'parents', or 'all'")


def _find_classic_confounders(dag, treatment_node, outcome_node):
    treatment_parents = set(dag.predecessors(treatment_node))
    outcome_parents = set(dag.predecessors(outcome_node))
    outcome_parents.discard(treatment_node)
    return list(treatment_parents.intersection(outcome_parents))


def _find_parent_confounders(dag, treatment_node, outcome_node):
    treatment_parents = set(dag.predecessors(treatment_node))
    outcome_parents = set(dag.predecessors(outcome_node))
    return list(treatment_parents & outcome_parents)


def _find_backdoor_confounders(dag, treatment_node, outcome_node):
    dag_copy = dag.copy()
    if dag_copy.has_edge(treatment_node, outcome_node):
        dag_copy.remove_edge(treatment_node, outcome_node)
    undirected_dag = dag_copy.to_undirected()
    backdoor_nodes = set()

    try:
        if nx.has_path(undirected_dag, treatment_node, outcome_node):
            all_paths = list(
                nx.all_simple_paths(undirected_dag, treatment_node, outcome_node)
            )
            for path in all_paths:
                if len(path) > 2:
                    first_edge_backwards = dag.has_edge(path[1], path[0])
                    if first_edge_backwards:
                        backdoor_nodes.update(path[1:-1])
    except Exception:
        pass
    return list(backdoor_nodes)


def compose_group_graphs(group_graphs):
    global_graph = nx.DiGraph()
    for graph in group_graphs.values():
        global_graph = nx.compose(global_graph, graph)
    return global_graph


def identify_biomarkers_from_group_graphs(
    group_graphs,
    expression_df,
    gene_names,
    outcome="DISEASE",
    sink_node_name="DISEASE",
    confounder_method="classic",
    include_n_confounders=False,
    sort_by_abs_ace=False,
    skip_missing_genes=False,
    run_causal_analysis_fn=None,
):
    if run_causal_analysis_fn is None:
        from .causal import run_causal_analysis as run_causal_analysis_fn

    global_graph = compose_group_graphs(group_graphs)
    global_graph_with_outcome = add_sink_node_to_graph(
        global_graph, sink_node_name=sink_node_name
    )

    biomarkers = []
    ace_values = []
    confounder_counts = []
    for treatment in gene_names:
        if skip_missing_genes and treatment not in global_graph_with_outcome:
            continue
        confounders = find_confounders(
            global_graph_with_outcome,
            treatment,
            outcome,
            method=confounder_method,
        )
        causal_results = run_causal_analysis_fn(
            dag=global_graph_with_outcome,
            data=expression_df,
            treatment=treatment,
            outcome=outcome,
            confounders=confounders,
        )
        if not causal_results.get("success", False):
            continue
        results = causal_results.get("results")
        if results is None:
            continue
        adjustment_formula = results["adjustment_formula"]
        if adjustment_formula != 0:
            biomarkers.append(treatment)
            ace_values.append(adjustment_formula)
            if include_n_confounders:
                confounder_counts.append(len(confounders))

    biomarker_data = {"gene": biomarkers, "ACE": ace_values}
    if include_n_confounders:
        biomarker_data["n_confounders"] = confounder_counts
    biomarkers_df = pd.DataFrame(biomarker_data)
    if sort_by_abs_ace and not biomarkers_df.empty:
        biomarkers_df = biomarkers_df.sort_values(
            by="ACE", key=lambda series: series.abs(), ascending=False
        ).reset_index(drop=True)
    return biomarkers_df
