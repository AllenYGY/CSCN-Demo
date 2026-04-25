from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

DATA_SET = "GSE132465"
RUN_NAME = f"{DATA_SET}_tumor_epithelial_vs_normal_epithelial"
GROUP_TO_LABEL = {"normal": 0, "tumor": 1}
GROUP_FILTERS = {
    "normal": {
        "class_labels": {"Normal"},
        "cell_types": {"Epithelial cells"},
    },
    "tumor": {
        "class_labels": {"Tumor"},
        "cell_types": {"Epithelial cells"},
    },
}
DEFAULT_SAMPLE_SIZE = 1000
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
        REPO_ROOT / "data" / "GSE132465",
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
    return output_dir / "deseq2_tumor_epithelial_vs_normal_epithelial_top150_genes.csv"


def parse_args():
    parser = argparse.ArgumentParser(description=(
        "Run CSCN biomarker analysis for GSE132465 using tumor epithelial cells "
        "versus normal epithelial cells."))
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=resolve_default_data_dir(),
        help=
        ("Dataset directory that contains GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz "
         "and GSE132465_GEO_processed_CRC_10X_natural_log_TPM_matrix.txt.gz."),
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
        help="Number of cells to sample per group.",
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
        help="Only use the first N genes from the top-gene list.",
    )
    parser.add_argument(
        "--run-name",
        default=RUN_NAME,
        help="Run prefix used for saved matrices, DAGs, and biomarkers.",
    )
    return parser.parse_args()


def run_group_cscn(data_dir: Path, run_name: str, group: str, matrix,
                   max_workers: int):
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
    log(f"running run_pc_concurrently() for {group}; "
        f"DAG files are written incrementally as each task completes")
    cscn.run_pc_concurrently(
        max_workers=max_workers,
        progress_interval=100,
        progress_label=group,
    )
    CSCN.save_to_file(cscn, cscn_path)
    dag_count = len(list(dag_dir.glob("result_*.pkl")))
    log(f"saved CSCN object: {cscn_path}")
    log(f"DAG files written for {group}: {dag_count}")


def log_group_summary(summary_by_group):
    for group, summary in summary_by_group.items():
        log(f"eligible cells for {group}: {summary['total_cells']}")
        sample_counts = summary["sample_counts"]
        if sample_counts:
            sorted_samples = sorted(sample_counts.items(),
                                    key=lambda item: (-item[1], item[0]))
            preview = ", ".join(f"{name}={count}"
                                for name, count in sorted_samples[:10])
            log(f"{group} sample counts: {preview}")
        subtype_counts = summary["cell_subtype_counts"]
        if subtype_counts:
            sorted_subtypes = sorted(
                subtype_counts.items(),
                key=lambda item: (-item[1], item[0]),
            )
            preview = ", ".join(f"{name}={count}"
                                for name, count in sorted_subtypes[:10])
            log(f"{group} subtype counts: {preview}")


def main():
    args = parse_args()

    from biomarker.datasets import (build_expression_df, load_gene_names,
                                    load_gse132465_grouped_cells,
                                    read_expression_for_sampled_cells,
                                    sample_cells_by_group,
                                    save_prepared_inputs)

    data_dir = args.data_dir.resolve()
    output_dir = data_dir / "output_deseq"
    annotation_path = data_dir / "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
    expression_matrix_path = data_dir / "GSE132465_GEO_processed_CRC_10X_natural_log_TPM_matrix.txt.gz"
    gene_list_path = (args.gene_list_path.resolve() if args.gene_list_path
                      is not None else default_gene_list_path(output_dir))
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
    log(f"gene limit: {args.gene_limit}")
    validate_required_file(annotation_path, "cell annotation")
    validate_required_file(expression_matrix_path, "natural log TPM matrix")
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

    log("loading grouped cells from cell annotation")
    eligible_cells, summary_by_group = load_gse132465_grouped_cells(
        annotation_path=annotation_path,
        group_filters=GROUP_FILTERS,
    )
    log_group_summary(summary_by_group)

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
    log("reading selected genes from natural log TPM matrix; this can take a while"
        )
    matrices, used_genes = read_expression_for_sampled_cells(
        counts_path=expression_matrix_path,
        sampled_cells=sampled_cells,
        top_genes=top_genes,
        delimiter="\t",
    )
    for group, matrix in matrices.items():
        log(f"natural log TPM {group} matrix shape: {matrix.shape}, "
            f"min={matrix.min():.3f}, max={matrix.max():.3f}")

    log_stage("Save Prepared Inputs")
    save_prepared_inputs(
        output_dir=output_dir,
        dataset_name=run_name,
        sampled_cells=sampled_cells,
        matrices=matrices,
        used_genes=used_genes,
        used_genes_filename=f"{run_name}_top{len(used_genes)}_genes_used.csv",
    )

    for group, matrix in matrices.items():
        log(f"saved {group} matrix shape: {matrix.shape}, "
            f"dtype={matrix.dtype}, min={matrix.min():.3f}, max={matrix.max():.3f}"
            )
        log(f"saved matrix: {output_dir / f'{run_name}_{group}.npy'}")
        log(f"saved sampled cells: {output_dir / f'{run_name}_{group}_sampled_cells.csv'}"
            )
    log(f"genes used after extraction: {len(used_genes)}")
    log(f"saved used-gene list: {output_dir / f'{run_name}_top{len(used_genes)}_genes_used.csv'}"
        )

    if args.prepare_only:
        log("prepare-only mode enabled; stopping before CSCN")
        return

    from biomarker.causal import run_causal_analysis
    from biomarker.cscn import CSCN
    from biomarker.datasets import load_saved_group_graphs
    from biomarker.graph_utils import (get_global_graph,
                                       identify_biomarkers_from_group_graphs,
                                       map_node_id_to_gene)

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
        log(f"group graph {group}: nodes={graph.number_of_nodes()}, "
            f"edges={graph.number_of_edges()}")

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
