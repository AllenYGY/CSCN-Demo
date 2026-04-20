import os

import matplotlib.pyplot as plt
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


def get_global_graph(dags):
    global_graph = nx.DiGraph()
    for _, dag in dags:
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
