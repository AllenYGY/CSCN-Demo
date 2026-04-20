from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


DATA_SET = "GSE138852"
GROUP_TO_LABEL = {"ct": 0, "AD": 1}
EXCLUDED_CELL_TYPES = {"doublet", "unID"}
DEFAULT_SAMPLE_SIZE = 3000
DEFAULT_RANDOM_SEED = 42
DEFAULT_MAX_WORKERS = min(8, os.cpu_count() or 1)
DEFAULT_USE_BITMAP = True


def log(message):
    print(f"[{DATA_SET}] {message}")


def log_stage(title):
    print()
    print(f"=== {title} ===")


def validate_required_file(path: Path, description: str):
    if not path.exists():
        raise FileNotFoundError(f"Missing {description}: {path}")
    log(f"{description}: {path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run CSCN biomarker analysis for GSE138852 using sampled AD/ct nuclei."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=REPO_ROOT / "data" / DATA_SET,
        help="Dataset directory that contains GSE138852_counts.csv.gz and output_deseq/.",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=DEFAULT_SAMPLE_SIZE,
        help="Number of nuclei to sample per group.",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=DEFAULT_RANDOM_SEED,
        help="Random seed for reproducible sampling.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=DEFAULT_MAX_WORKERS,
        help="Worker count passed to CSCN.run_pc_concurrently().",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Only prepare sampled .npy matrices and stop before CSCN.",
    )
    parser.add_argument(
        "--gene-limit",
        type=int,
        default=None,
        help="Only use the first N genes from the DESeq2 top-gene list.",
    )
    return parser.parse_args()


def run_group_cscn(data_dir: Path, group: str, matrix, max_workers: int):
    from biomarker.cscn import CSCN

    dag_dir = data_dir / "DAG" / group
    dag_dir.mkdir(parents=True, exist_ok=True)
    cscn_path = data_dir / f"{DATA_SET}_{group}_cscn"

    log_stage(f"CSCN {group}")
    log(f"group: {group}")
    log(f"matrix shape: {matrix.shape}")
    log(f"DAG output dir: {dag_dir}")
    log(f"CSCN object path: {cscn_path}")
    log(f"max_workers: {max_workers}")
    log(f"use bitmap: {DEFAULT_USE_BITMAP}")

    cscn = CSCN(
        output_dir=str(dag_dir),
        sigmoid_score=0.1,
        significance_level=0.01,
        max_cond_vars=20,
        use_bitmap=DEFAULT_USE_BITMAP,
        debug=False,
        show_progress=False,
        progress_interval=100,
    )
    log(f"running run_core() for {group}")
    cscn.run_core(matrix, usingNMF=False)
    log(
        f"running run_pc_concurrently() for {group}; "
        f"DAG files are written incrementally as each task completes"
    )
    cscn.run_pc_concurrently(
        max_workers=max_workers,
        progress_interval=100,
        progress_label=group,
    )
    CSCN.save_to_file(cscn, cscn_path)
    dag_count = len(list(dag_dir.glob("result_*.pkl")))
    log(f"saved CSCN object: {cscn_path}")
    log(f"DAG files written for {group}: {dag_count}")


def main():
    args = parse_args()

    from biomarker.causal import run_causal_analysis
    from biomarker.cscn import CSCN
    from biomarker.datasets import (
        build_expression_df,
        load_gene_names,
        load_gse138852_eligible_cells,
        load_saved_group_graphs,
        normalize_log1p,
        read_expression_for_sampled_cells,
        sample_cells_by_group,
        save_prepared_inputs,
    )
    from biomarker.graph_utils import (
        add_sink_node_to_graph,
        find_confounders,
        get_global_graph,
        identify_biomarkers_from_group_graphs,
        map_node_id_to_gene,
    )
    data_dir = args.data_dir.resolve()
    output_dir = data_dir / "output_deseq"

    counts_path = data_dir / "GSE138852_counts.csv.gz"
    covariates_path = data_dir / "GSE138852_covariates.csv.gz"
    gene_list_path = output_dir / "deseq2_ad_vs_ct_top150_genes.csv"
    biomarker_path = data_dir / "Biomarkers.csv"

    log_stage("Configuration")
    log(f"dataset dir: {data_dir}")
    log(f"output dir: {output_dir}")
    log(f"sample size: {args.sample_size}")
    log(f"random seed: {args.random_seed}")
    log(f"max workers: {args.max_workers}")
    log(f"prepare only: {args.prepare_only}")
    log(f"gene limit: {args.gene_limit}")
    validate_required_file(counts_path, "counts matrix")
    validate_required_file(covariates_path, "covariates table")
    validate_required_file(gene_list_path, "top-gene list")

    log_stage("Load Inputs")
    top_genes = load_gene_names(gene_list_path)
    if args.gene_limit is not None:
        if args.gene_limit <= 0:
            raise ValueError("--gene-limit must be a positive integer")
        requested_gene_count = min(args.gene_limit, len(top_genes))
        top_genes = top_genes[:requested_gene_count]
        log(f"using the first {requested_gene_count} genes for this run")
    log(f"top genes requested: {len(top_genes)}")
    log(f"first 5 genes: {top_genes[:5]}")

    log("loading eligible cells from covariates")
    eligible_cells = load_gse138852_eligible_cells(
        covariates_path=covariates_path,
        group_to_label=GROUP_TO_LABEL,
        excluded_cell_types=EXCLUDED_CELL_TYPES,
    )
    for group, cell_ids in eligible_cells.items():
        log(f"eligible cells for {group}: {len(cell_ids)}")

    log_stage("Sample Cells")
    sampled_cells = sample_cells_by_group(
        eligible_cells,
        sample_size=args.sample_size,
        random_seed=args.random_seed,
    )
    for group, cell_ids in sampled_cells.items():
        log(f"sampled cells for {group}: {len(cell_ids)}")
        log(f"first 3 sampled {group} cells: {cell_ids[:3]}")

    log_stage("Extract Expression")
    log("reading selected genes from counts matrix; this can take a while")
    raw_matrices, used_genes = read_expression_for_sampled_cells(
        counts_path=counts_path,
        sampled_cells=sampled_cells,
        top_genes=top_genes,
    )
    for group, matrix in raw_matrices.items():
        log(
            f"raw {group} matrix shape: {matrix.shape}, "
            f"min={matrix.min():.3f}, max={matrix.max():.3f}"
        )

    log_stage("Normalize And Save")
    matrices = {group: normalize_log1p(matrix) for group, matrix in raw_matrices.items()}
    save_prepared_inputs(
        output_dir=output_dir,
        dataset_name=DATA_SET,
        sampled_cells=sampled_cells,
        matrices=matrices,
        used_genes=used_genes,
        used_genes_filename=f"{DATA_SET}_top{len(used_genes)}_genes_used.csv",
    )

    for group, matrix in matrices.items():
        log(
            f"normalized {group} matrix shape: {matrix.shape}, "
            f"dtype={matrix.dtype}, min={matrix.min():.3f}, max={matrix.max():.3f}"
        )
        log(f"saved matrix: {output_dir / f'{DATA_SET}_{group}.npy'}")
        log(f"saved sampled cells: {output_dir / f'{DATA_SET}_{group}_sampled_cells.csv'}")
    log(f"genes used after extraction: {len(used_genes)}")
    log(f"saved used-gene list: {output_dir / f'{DATA_SET}_top{len(used_genes)}_genes_used.csv'}")

    if args.prepare_only:
        log("prepare-only mode enabled; stopping before CSCN")
        return

    log_stage("Run CSCN")
    for group in GROUP_TO_LABEL:
        run_group_cscn(data_dir, group, matrices[group], max_workers=args.max_workers)

    log_stage("Load Graphs")
    expression_df = build_expression_df(matrices, used_genes, GROUP_TO_LABEL)
    log(f"combined expression dataframe shape: {expression_df.shape}")
    group_graphs = load_saved_group_graphs(
        data_dir=data_dir,
        dataset_name=DATA_SET,
        groups=GROUP_TO_LABEL.keys(),
        gene_names=used_genes,
        cscn_cls=CSCN,
        map_node_id_to_gene_fn=map_node_id_to_gene,
        get_global_graph_fn=get_global_graph,
    )
    for group, graph in group_graphs.items():
        log(
            f"group graph {group}: nodes={graph.number_of_nodes()}, "
            f"edges={graph.number_of_edges()}"
        )

    log_stage("Identify Biomarkers")
    biomarkers_df = identify_biomarkers_from_group_graphs(
        group_graphs=group_graphs,
        expression_df=expression_df,
        gene_names=used_genes,
        outcome="DISEASE",
        sink_node_name="DISEASE",
        confounder_method="classic",
        include_n_confounders=True,
        sort_by_abs_ace=True,
        skip_missing_genes=True,
        run_causal_analysis_fn=run_causal_analysis,
    )
    biomarkers_df.to_csv(biomarker_path, index=False)
    log(f"saved biomarkers: {biomarker_path}")
    log(f"biomarker count: {len(biomarkers_df)}")
    if not biomarkers_df.empty:
        log("top 5 biomarkers:")
        print(biomarkers_df.head(5).to_string(index=False))


if __name__ == "__main__":
    main()
