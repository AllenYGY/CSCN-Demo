from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


DEFAULT_DATASET_NAME = "BreastTumer"
DEFAULT_GROUPS = ("epithelial_cells_normal", "epithelial_cells_cancer")
DEFAULT_GENE_LIST = "output_deseq/deseq2_cancer_vs_epi_top150_genes.csv"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Rebuild a biomarker DAG from saved CSCN objects and draw a highlighted "
            "DAG PNG. Defaults match evaluate_ATE_BreastTumer.ipynb."
        )
    )
    parser.add_argument(
        "--dataset-name",
        default=DEFAULT_DATASET_NAME,
        help="Dataset name used for logging and some default path resolution.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=REPO_ROOT / "data" / DEFAULT_DATASET_NAME,
        help="Dataset directory containing Biomarkers.csv and saved CSCN objects.",
    )
    parser.add_argument(
        "--groups",
        nargs="+",
        default=list(DEFAULT_GROUPS),
        help="Group suffixes used in <cscn-prefix>_<group>_cscn object names.",
    )
    parser.add_argument(
        "--cscn-prefix",
        default=None,
        help=(
            "Prefix used in saved CSCN object names. Final object path is "
            "<data-dir>/<cscn-prefix>_<group>_cscn. Defaults to --dataset-name."
        ),
    )
    parser.add_argument(
        "--gene-list-path",
        type=Path,
        default=None,
        help="Optional path to the top-gene CSV used for id->gene mapping.",
    )
    parser.add_argument(
        "--biomarkers-path",
        type=Path,
        default=None,
        help="Optional path to biomarker CSV. Defaults to an inferred Biomarkers*.csv.",
    )
    parser.add_argument(
        "--output-path",
        type=Path,
        default=None,
        help="Optional output PNG path. Defaults to <data-dir>/<cscn-prefix>_Global_Biomarker_DAG.png.",
    )
    parser.add_argument(
        "--outcome-node",
        default="DISEASE",
        help="Outcome node name to append to the graph. Default: DISEASE",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional plot title. Defaults to '<dataset-name> Biomarker DAG'.",
    )
    return parser.parse_args()


def log(dataset_name: str, message: str) -> None:
    print(f"[{dataset_name}] {message}", flush=True)


def resolve_cscn_prefix(dataset_name: str, cscn_prefix: str | None) -> str:
    return cscn_prefix or dataset_name


def resolve_default_gene_list_path(
    data_dir: Path,
    dataset_name: str,
    cscn_prefix: str,
) -> Path:
    candidates = []
    if dataset_name == DEFAULT_DATASET_NAME:
        candidates.append(data_dir / DEFAULT_GENE_LIST)
    if cscn_prefix.startswith(f"{dataset_name}_"):
        run_slug = cscn_prefix[len(dataset_name) + 1 :]
        candidates.append(data_dir / "output_deseq" / f"deseq2_{run_slug}_top150_genes.csv")
        candidates.append(data_dir / "output_deseq" / f"{cscn_prefix}_top150_genes_used.csv")
    candidates.append(data_dir / "output_deseq" / f"{cscn_prefix}_top150_genes_used.csv")

    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()

    matched = sorted((data_dir / "output_deseq").glob("*top150_genes*.csv"))
    if len(matched) == 1:
        return matched[0].resolve()

    raise FileNotFoundError(
        "Could not infer gene list CSV. Pass --gene-list-path explicitly."
    )


def resolve_default_biomarkers_path(data_dir: Path, cscn_prefix: str) -> Path:
    candidates = [
        data_dir / "Biomarkers.csv",
        data_dir / f"Biomarkers_{cscn_prefix}.csv",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()

    matched = sorted(data_dir.glob("Biomarkers*.csv"))
    if len(matched) == 1:
        return matched[0].resolve()

    raise FileNotFoundError(
        "Could not infer biomarker CSV. Pass --biomarkers-path explicitly."
    )


def load_gene_names(gene_list_path: Path) -> list[str]:
    import pandas as pd

    gene_info_df = pd.read_csv(gene_list_path)
    if "gene" in gene_info_df.columns:
        genes = gene_info_df["gene"]
    elif "Unnamed: 0" in gene_info_df.columns:
        genes = gene_info_df["Unnamed: 0"]
    else:
        raise ValueError(
            f"Missing gene column in {gene_list_path}. Expected `gene` or `Unnamed: 0`."
        )
    return genes.dropna().astype(str).tolist()


def load_biomarkers(biomarkers_path: Path) -> list[str]:
    import pandas as pd

    biomarkers_df = pd.read_csv(biomarkers_path)
    if "gene" not in biomarkers_df.columns:
        raise ValueError(f"Missing `gene` column in {biomarkers_path}.")
    return biomarkers_df["gene"].dropna().astype(str).tolist()


def load_group_graph(data_dir: Path, cscn_prefix: str, group: str, gene_names: list[str]):
    from biomarker.cscn import CSCN
    from biomarker.graph_utils import get_global_graph, map_node_id_to_gene_directed

    cscn_path = data_dir / f"{cscn_prefix}_{group}_cscn"
    if not cscn_path.is_file():
        raise FileNotFoundError(
            f"Missing CSCN object for group {group}: {cscn_path}"
        )
    cscn = CSCN.load_from_file(cscn_path)
    dags = cscn.load_all_dags()
    id2gene = {idx: gene for idx, gene in enumerate(gene_names)}
    directed_dags = map_node_id_to_gene_directed(dags, id2gene)
    return get_global_graph(directed_dags)


def main() -> None:
    args = parse_args()

    from biomarker.graph_utils import (
        build_biomarker_dag,
        draw_global_network_highlighted,
        get_global_graph,
    )

    data_dir = args.data_dir.resolve()
    dataset_name = args.dataset_name
    cscn_prefix = resolve_cscn_prefix(dataset_name, args.cscn_prefix)
    gene_list_path = (
        args.gene_list_path.resolve()
        if args.gene_list_path is not None
        else resolve_default_gene_list_path(data_dir, dataset_name, cscn_prefix)
    )
    biomarkers_path = (
        args.biomarkers_path.resolve()
        if args.biomarkers_path is not None
        else resolve_default_biomarkers_path(data_dir, cscn_prefix)
    )
    output_path = (
        args.output_path.resolve()
        if args.output_path is not None
        else (data_dir / f"{cscn_prefix}_Global_Biomarker_DAG.png").resolve()
    )
    title = args.title or f"{dataset_name} Biomarker DAG"

    log(dataset_name, f"dataset dir: {data_dir}")
    log(dataset_name, f"cscn prefix: {cscn_prefix}")
    log(dataset_name, f"groups: {args.groups}")
    log(dataset_name, f"gene list path: {gene_list_path}")
    log(dataset_name, f"biomarkers path: {biomarkers_path}")
    log(dataset_name, f"output path: {output_path}")

    if not gene_list_path.is_file():
        raise FileNotFoundError(f"Missing gene list file: {gene_list_path}")
    if not biomarkers_path.is_file():
        raise FileNotFoundError(f"Missing biomarkers file: {biomarkers_path}")

    gene_names = load_gene_names(gene_list_path)
    biomarkers = load_biomarkers(biomarkers_path)
    log(dataset_name, f"loaded {len(gene_names)} genes")
    log(dataset_name, f"loaded {len(biomarkers)} biomarkers")

    group_graphs = []
    for group in args.groups:
        graph = load_group_graph(data_dir, cscn_prefix, group, gene_names)
        log(
            dataset_name,
            f"group {group}: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges",
        )
        group_graphs.append(graph)

    global_graph = get_global_graph(list(enumerate(group_graphs)))
    log(
        dataset_name,
        f"global graph: {global_graph.number_of_nodes()} nodes, {global_graph.number_of_edges()} edges",
    )

    biomarker_dag, confounders_by_treatment = build_biomarker_dag(
        global_graph=global_graph,
        highlighted_treatments=biomarkers,
        outcome_node=args.outcome_node,
        candidate_treatments=gene_names,
    )
    confounder_count = len({node for nodes in confounders_by_treatment.values() for node in nodes})
    log(
        dataset_name,
        f"biomarker DAG: {biomarker_dag.number_of_nodes()} nodes, "
        f"{biomarker_dag.number_of_edges()} edges, {confounder_count} unique confounders"
    )

    saved_path = draw_global_network_highlighted(
        biomarker_dag,
        treatment_nodes=biomarkers,
        outcome_nodes=[args.outcome_node],
        save_path=str(output_path),
        title=title,
    )
    log(dataset_name, f"saved DAG plot to {saved_path}")


if __name__ == "__main__":
    main()
