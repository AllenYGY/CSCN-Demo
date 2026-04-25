from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


DATA_SET = "GSE159115"
RUN_SLUG = "ccrcc_tumor_vs_ptb_ptc_normal"
RUN_NAME = f"{DATA_SET}_{RUN_SLUG}"
GROUP_TO_LABEL = {"normal": 0, "tumor": 1}
DEFAULT_SAMPLE_SIZE = 100
DEFAULT_RANDOM_SEED = 42
DEFAULT_MAX_WORKERS = min(8, os.cpu_count() or 1)
DEFAULT_USE_BITMAP = True


def log(message):
    print(f"[{DATA_SET}] {message}", flush=True)


def log_stage(title):
    print(flush=True)
    print(f"=== {title} ===", flush=True)


def resolve_default_data_dir():
    candidates = [
        REPO_ROOT / "data" / DATA_SET,
        REPO_ROOT / "data" / DATA_SET.lower(),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def validate_required_file(path: Path, description: str):
    if not path.exists():
        raise FileNotFoundError(f"Missing {description}: {path}")
    log(f"{description}: {path}")


def default_gene_list_path(output_dir: Path):
    return output_dir / "deseq2_ccrcc_tumor_vs_ptb_ptc_normal_top150_genes.csv"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run CSCN biomarker analysis for GSE159115 using paired ccRCC tumor "
            "cells versus normal PT-B/PT-C cells."
        )
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=resolve_default_data_dir(),
        help="Dataset directory containing GSE159115 GEO files.",
    )
    parser.add_argument(
        "--gene-list-path",
        type=Path,
        default=None,
        help="Optional explicit path to the top-gene CSV used to seed CSCN.",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=DEFAULT_SAMPLE_SIZE,
        help="Number of cells to sample per group before CSCN.",
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
        help="Only prepare sampled matrices and stop before CSCN.",
    )
    parser.add_argument(
        "--run-name",
        default=RUN_NAME,
        help="Run prefix used for saved matrices, DAGs, and biomarkers.",
    )
    return parser.parse_args()


def run_group_cscn(data_dir: Path, run_name: str, group: str, matrix, max_workers: int):
    from biomarker.cscn import CSCN

    dag_dir = data_dir / "DAG" / run_name / group
    dag_dir.mkdir(parents=True, exist_ok=True)
    cscn_path = data_dir / f"{run_name}_{group}_cscn"

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


def prepare_expression_inputs(data_dir: Path, gene_list_path: Path, sample_size: int, random_seed: int, run_name: str):
    helper_path = REPO_ROOT / "scripts" / "biomarker" / "prepare_GSE159115_expression.R"
    cmd = [
        "Rscript",
        str(helper_path),
        "--data-dir",
        str(data_dir),
        "--gene-list-path",
        str(gene_list_path),
        "--sample-size",
        str(sample_size),
        "--random-seed",
        str(random_seed),
        "--run-name",
        run_name,
    ]
    log("preparing sampled expression matrices via R helper")
    subprocess.run(cmd, check=True)


def load_group_matrix(path: Path):
    df = pd.read_csv(path)
    if "cell_id" not in df.columns:
        raise ValueError(f"Missing cell_id column in {path}")
    gene_columns = [column for column in df.columns if column != "cell_id"]
    matrix = df[gene_columns].to_numpy(dtype=float)
    cell_ids = df["cell_id"].astype(str).tolist()
    return cell_ids, gene_columns, matrix


def main():
    args = parse_args()

    from biomarker.causal import run_causal_analysis
    from biomarker.datasets import build_expression_df, load_saved_group_graphs, save_prepared_inputs
    from biomarker.graph_utils import (
        get_global_graph,
        identify_biomarkers_from_group_graphs,
        map_node_id_to_gene,
    )
    from biomarker.cscn import CSCN

    data_dir = args.data_dir.resolve()
    output_dir = data_dir / "output_deseq"
    gene_list_path = (
        args.gene_list_path.resolve()
        if args.gene_list_path is not None
        else default_gene_list_path(output_dir)
    )
    run_name = args.run_name
    biomarker_path = data_dir / f"Biomarkers_{run_name}.csv"

    log_stage("Configuration")
    log(f"dataset dir: {data_dir}")
    log(f"output dir: {output_dir}")
    log(f"run name: {run_name}")
    log(f"sample size: {args.sample_size}")
    log(f"random seed: {args.random_seed}")
    log(f"max workers: {args.max_workers}")
    log(f"prepare only: {args.prepare_only}")
    validate_required_file(gene_list_path, "top-gene list")

    prepare_expression_inputs(
        data_dir=data_dir,
        gene_list_path=gene_list_path,
        sample_size=args.sample_size,
        random_seed=args.random_seed,
        run_name=run_name,
    )

    normal_expr_path = output_dir / f"{run_name}_normal_expression.csv.gz"
    tumor_expr_path = output_dir / f"{run_name}_tumor_expression.csv.gz"
    used_gene_paths = sorted(output_dir.glob(f"{run_name}_top*_genes_used.csv"))
    if len(used_gene_paths) != 1:
        raise FileNotFoundError(
            f"Expected exactly one used-gene file for {run_name}, found {len(used_gene_paths)}"
        )
    used_genes_path = used_gene_paths[0]

    validate_required_file(normal_expr_path, "prepared normal expression matrix")
    validate_required_file(tumor_expr_path, "prepared tumor expression matrix")
    validate_required_file(used_genes_path, "used-gene list")

    log_stage("Load Prepared Matrices")
    normal_cells, normal_genes, normal_matrix = load_group_matrix(normal_expr_path)
    tumor_cells, tumor_genes, tumor_matrix = load_group_matrix(tumor_expr_path)
    if normal_genes != tumor_genes:
        raise ValueError("Normal and tumor prepared matrices use different gene columns")
    used_genes = normal_genes
    matrices = {"normal": normal_matrix, "tumor": tumor_matrix}
    sampled_cells = {"normal": normal_cells, "tumor": tumor_cells}

    for group, matrix in matrices.items():
        log(
            f"prepared {group} matrix shape: {matrix.shape}, "
            f"min={matrix.min():.3f}, max={matrix.max():.3f}"
        )

    log_stage("Save Prepared Inputs")
    save_prepared_inputs(
        output_dir=output_dir,
        dataset_name=run_name,
        sampled_cells=sampled_cells,
        matrices=matrices,
        used_genes=used_genes,
        used_genes_filename=used_genes_path.name,
    )

    if args.prepare_only:
        log("prepare-only mode enabled; stopping before CSCN")
        return

    log_stage("Run CSCN")
    for group in GROUP_TO_LABEL:
        run_group_cscn(
            data_dir=data_dir,
            run_name=run_name,
            group=group,
            matrix=matrices[group],
            max_workers=args.max_workers,
        )

    log_stage("Load Graphs")
    expression_df = build_expression_df(matrices, used_genes, GROUP_TO_LABEL)
    log(f"combined expression dataframe shape: {expression_df.shape}")
    group_graphs = load_saved_group_graphs(
        data_dir=data_dir,
        dataset_name=run_name,
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
        print(biomarkers_df.head(5).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
