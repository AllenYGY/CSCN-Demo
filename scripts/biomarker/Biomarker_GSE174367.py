from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
os.environ.setdefault("MPLCONFIGDIR", "/tmp/cscn-matplotlib-cache")

from biomarker.causal import run_causal_analysis
from biomarker.cscn import CSCN
from biomarker.datasets import (
    build_expression_df,
    load_gene_names,
    load_saved_group_graphs,
    normalize_log1p,
    sample_cells_by_group,
    save_prepared_inputs,
)
from biomarker.graph_utils import (
    get_global_graph,
    identify_biomarkers_from_group_graphs,
    map_node_id_to_gene,
)
from gse174367_utils import (
    DATA_SET,
    DEFAULT_CASE_GROUP,
    DEFAULT_CONTROL_GROUP,
    build_run_slug,
    extract_sampled_expression_from_h5,
    filter_metadata_for_cell_type,
    load_gse174367_metadata,
    resolve_default_sample_size,
    summarize_group_cells,
)


DEFAULT_RANDOM_SEED = 42
DEFAULT_MAX_WORKERS = min(8, os.cpu_count() or 1)
DEFAULT_USE_BITMAP = True


def log(message: str) -> None:
    print(f"[{DATA_SET}] {message}", flush=True)


def log_stage(title: str) -> None:
    print(flush=True)
    print(f"=== {title} ===", flush=True)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run CSCN biomarker analysis for GSE174367 within one cell type, "
            "for example ODC AD vs Control."
        )
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=REPO_ROOT / "data" / DATA_SET,
        help="Dataset directory containing the GSE174367 snRNA-seq files.",
    )
    parser.add_argument(
        "--cell-type",
        default="ODC",
        help="Cell.Type value to analyze. Default: ODC",
    )
    parser.add_argument(
        "--case-group",
        default=DEFAULT_CASE_GROUP,
        help="Case diagnosis label. Default: AD",
    )
    parser.add_argument(
        "--control-group",
        default=DEFAULT_CONTROL_GROUP,
        help="Control diagnosis label. Default: Control",
    )
    parser.add_argument(
        "--gene-list-path",
        type=Path,
        default=None,
        help="Optional explicit path to the DESeq2 top-gene CSV used to seed CSCN.",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=None,
        help=(
            "Number of cells to sample per group. Default: auto, "
            "min(3000, smaller group size)."
        ),
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
    parser.add_argument(
        "--min-cells-per-sample",
        type=int,
        default=20,
        help="Minimum cells required per donor before the biomarker sampling step.",
    )
    return parser.parse_args()


def run_group_cscn(data_dir: Path, run_name: str, group: str, matrix, max_workers: int) -> None:
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
    log(f"running run_pc_concurrently() for {group}")
    cscn.run_pc_concurrently(
        max_workers=max_workers,
        progress_interval=100,
        progress_label=group,
    )
    CSCN.save_to_file(cscn, cscn_path)
    dag_count = len(list(dag_dir.glob("result_*.pkl")))
    log(f"saved CSCN object: {cscn_path}")
    log(f"DAG files written for {group}: {dag_count}")


def main() -> None:
    args = parse_args()
    if args.min_cells_per_sample <= 0:
        raise ValueError("--min-cells-per-sample must be a positive integer")

    run_slug = build_run_slug(
        cell_type=args.cell_type,
        case_group=args.case_group,
        control_group=args.control_group,
    )
    run_name = f"{DATA_SET}_{run_slug}"
    group_to_label = {
        args.control_group: 0,
        args.case_group: 1,
    }

    data_dir = args.data_dir.resolve()
    output_dir = data_dir / "output_deseq"
    h5_path = data_dir / "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5"
    metadata_path = data_dir / "GSE174367_snRNA-seq_cell_meta.csv.gz"
    gene_list_path = (
        args.gene_list_path.resolve()
        if args.gene_list_path is not None
        else output_dir / f"deseq2_{run_slug}_top150_genes.csv"
    )
    biomarker_path = data_dir / f"Biomarkers_{run_name}.csv"

    log_stage("Configuration")
    log(f"dataset dir: {data_dir}")
    log(f"cell type: {args.cell_type}")
    log(f"run name: {run_name}")
    log(f"sample size: {args.sample_size or 'auto'}")
    log(f"random seed: {args.random_seed}")
    log(f"max workers: {args.max_workers}")
    log(f"prepare only: {args.prepare_only}")
    log(f"gene limit: {args.gene_limit}")
    log(f"gene list path: {gene_list_path}")

    log_stage("Load Inputs")
    metadata_df = load_gse174367_metadata(metadata_path)
    selected_df = filter_metadata_for_cell_type(
        metadata_df,
        cell_type=args.cell_type,
        case_group=args.case_group,
        control_group=args.control_group,
        min_cells_per_sample=args.min_cells_per_sample,
    )
    cells_by_group = summarize_group_cells(
        selected_df,
        case_group=args.case_group,
        control_group=args.control_group,
    )
    for group, cells in cells_by_group.items():
        log(f"eligible cells for {group}: {len(cells)}")

    top_genes = load_gene_names(gene_list_path)
    if args.gene_limit is not None:
        if args.gene_limit <= 0:
            raise ValueError("--gene-limit must be a positive integer")
        requested_gene_count = min(args.gene_limit, len(top_genes))
        top_genes = top_genes[:requested_gene_count]
        log(f"using the first {requested_gene_count} genes for this run")
    log(f"top genes requested: {len(top_genes)}")
    log(f"first 5 genes: {top_genes[:5]}")

    log_stage("Sample Cells")
    sample_size = args.sample_size
    if sample_size is None:
        sample_size = resolve_default_sample_size(cells_by_group, cap=3000)
        log(f"auto sample size selected: {sample_size}")
    sampled_cells = sample_cells_by_group(
        cells_by_group,
        sample_size=sample_size,
        random_seed=args.random_seed,
    )
    for group, cell_ids in sampled_cells.items():
        log(f"sampled cells for {group}: {len(cell_ids)}")
        log(f"first 3 sampled {group} cells: {cell_ids[:3]}")

    log_stage("Extract Expression")
    raw_matrices, used_genes = extract_sampled_expression_from_h5(
        h5_path=h5_path,
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
        dataset_name=run_name,
        sampled_cells=sampled_cells,
        matrices=matrices,
        used_genes=used_genes,
        used_genes_filename=f"{run_name}_top{len(used_genes)}_genes_used.csv",
    )
    for group, matrix in matrices.items():
        log(
            f"normalized {group} matrix shape: {matrix.shape}, "
            f"dtype={matrix.dtype}, min={matrix.min():.3f}, max={matrix.max():.3f}"
        )
        log(f"saved matrix: {output_dir / f'{run_name}_{group}.npy'}")
        log(f"saved sampled cells: {output_dir / f'{run_name}_{group}_sampled_cells.csv'}")
    log(f"genes used after extraction: {len(used_genes)}")
    log(f"saved used-gene list: {output_dir / f'{run_name}_top{len(used_genes)}_genes_used.csv'}")

    if args.prepare_only:
        log("prepare-only mode enabled; stopping before CSCN")
        return

    log_stage("Run CSCN")
    for group in (args.control_group, args.case_group):
        run_group_cscn(
            data_dir=data_dir,
            run_name=run_name,
            group=group,
            matrix=matrices[group],
            max_workers=args.max_workers,
        )

    log_stage("Load Graphs")
    expression_df = build_expression_df(matrices, used_genes, group_to_label)
    log(f"combined expression dataframe shape: {expression_df.shape}")
    group_graphs = load_saved_group_graphs(
        data_dir=data_dir,
        dataset_name=run_name,
        groups=(args.control_group, args.case_group),
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
