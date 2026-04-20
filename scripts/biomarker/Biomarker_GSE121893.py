from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


BASE_DATA_SET = "GSE121893"
DEFAULT_CASE_GROUP = "dHF"
DEFAULT_CONTROL_GROUP = "N"
DEFAULT_TARGET_SAMPLE_SIZE = 500
DEFAULT_RANDOM_SEED = 42
DEFAULT_MAX_WORKERS = min(8, os.cpu_count() or 1)
DEFAULT_SIGNIFICANCE_LEVEL = 0.05
DEFAULT_USE_BITMAP = True
EXCLUDED_CELL_TYPES = set()


def log(message):
    print(f"[{BASE_DATA_SET}] {message}")


def log_stage(title):
    print()
    print(f"=== {title} ===")


def validate_required_file(path: Path, description: str):
    if not path.exists():
        raise FileNotFoundError(f"Missing {description}: {path}")
    log(f"{description}: {path}")


def slugify(value: str):
    return "".join(char.lower() if char.isalnum() else "_" for char in value)


def default_run_name(case_group: str, control_group: str, region: str, cell_compartment: str):
    return (
        f"{BASE_DATA_SET}_{slugify(case_group)}_vs_{slugify(control_group)}_"
        f"{region.lower()}_{cell_compartment.lower()}"
    )


def default_gene_list_path(output_dir: Path, case_group: str, control_group: str):
    return output_dir / f"deseq2_{slugify(case_group)}_vs_{slugify(control_group)}_top150_genes.csv"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run CSCN biomarker analysis for GSE121893 using sampled cells from a "
            "selected disease comparison."
        )
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=REPO_ROOT / "data" / BASE_DATA_SET,
        help="Dataset directory that contains the raw GEO files plus GSE121893_covariates.csv.gz.",
    )
    parser.add_argument(
        "--case-group",
        default=DEFAULT_CASE_GROUP,
        help="Disease group to treat as the case cohort. Typical values: dHF or cHF.",
    )
    parser.add_argument(
        "--control-group",
        default=DEFAULT_CONTROL_GROUP,
        help="Control group to compare against. Typical values: N.",
    )
    parser.add_argument(
        "--gene-list-path",
        type=Path,
        default=None,
        help="Optional explicit path to the top-gene CSV used to seed CSCN.",
    )
    parser.add_argument(
        "--region",
        choices=("all", "LA", "LV"),
        default="all",
        help="Restrict eligible cells to a specific chamber region.",
    )
    parser.add_argument(
        "--cell-compartment",
        choices=("all", "CM", "NCM"),
        default="all",
        help="Restrict eligible cells to cardiomyocytes or non-cardiomyocytes.",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=None,
        help=(
            "Number of cells to sample per group. If omitted, uses the smaller of "
            f"{DEFAULT_TARGET_SAMPLE_SIZE} and the smallest eligible group."
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
        "--significance-level",
        type=float,
        default=DEFAULT_SIGNIFICANCE_LEVEL,
        help="Conditional independence significance level passed to CSCN.",
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
        default=None,
        help="Optional override for the run prefix used for saved matrices, DAGs, and biomarkers.",
    )
    return parser.parse_args()


def resolve_sample_size(requested_size, eligible_cells):
    available_counts = {group: len(cell_ids) for group, cell_ids in eligible_cells.items()}
    min_available = min(available_counts.values())
    if min_available <= 0:
        raise ValueError(f"At least one group has no eligible cells: {available_counts}")

    if requested_size is None:
        return min(DEFAULT_TARGET_SAMPLE_SIZE, min_available)

    if requested_size <= 0:
        raise ValueError("--sample-size must be a positive integer")

    if requested_size > min_available:
        raise ValueError(
            f"--sample-size={requested_size} exceeds the smallest eligible group size "
            f"({min_available}); available counts: {available_counts}"
        )

    return requested_size


def run_group_cscn(
    data_dir: Path,
    run_name: str,
    group: str,
    matrix,
    max_workers: int,
    significance_level: float,
):
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
    log(f"significance level: {significance_level}")
    log(f"use bitmap: {DEFAULT_USE_BITMAP}")

    cscn = CSCN(
        output_dir=str(dag_dir),
        sigmoid_score=0.1,
        significance_level=significance_level,
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
    if args.case_group == args.control_group:
        raise ValueError("--case-group and --control-group must be different")

    from biomarker.causal import run_causal_analysis
    from biomarker.cscn import CSCN
    from biomarker.datasets import (
        build_expression_df,
        load_gene_names,
        load_gse121893_eligible_cells,
        load_saved_group_graphs,
        normalize_log1p,
        read_expression_for_sampled_cells,
        sample_cells_by_group,
        save_prepared_inputs,
    )
    from biomarker.graph_utils import (
        get_global_graph,
        identify_biomarkers_from_group_graphs,
        map_node_id_to_gene,
    )

    data_dir = args.data_dir.resolve()
    output_dir = data_dir / "output_deseq"
    counts_path = data_dir / "GSE121893_human_heart_sc_umi.csv.gz"
    covariates_path = data_dir / "GSE121893_covariates.csv.gz"
    gene_list_path = (
        args.gene_list_path.resolve()
        if args.gene_list_path is not None
        else default_gene_list_path(output_dir, args.case_group, args.control_group)
    )
    run_name = args.run_name or default_run_name(
        case_group=args.case_group,
        control_group=args.control_group,
        region=args.region,
        cell_compartment=args.cell_compartment,
    )
    biomarker_path = data_dir / f"Biomarkers_{run_name}.csv"

    group_to_label = {
        args.control_group: 0,
        args.case_group: 1,
    }

    log_stage("Configuration")
    log(f"dataset dir: {data_dir}")
    log(f"output dir: {output_dir}")
    log(f"run name: {run_name}")
    log(f"case group: {args.case_group}")
    log(f"control group: {args.control_group}")
    log(f"region filter: {args.region}")
    log(f"cell compartment filter: {args.cell_compartment}")
    log(f"random seed: {args.random_seed}")
    log(f"max workers: {args.max_workers}")
    log(f"significance level: {args.significance_level}")
    log(f"prepare only: {args.prepare_only}")
    log(f"gene limit: {args.gene_limit}")
    validate_required_file(counts_path, "counts matrix")
    validate_required_file(covariates_path, "prepared covariates table")
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

    allowed_regions = None if args.region == "all" else [args.region]
    log("loading eligible cells from prepared covariates")
    eligible_cells = load_gse121893_eligible_cells(
        covariates_path=covariates_path,
        group_to_label=group_to_label,
        allowed_regions=allowed_regions,
        cell_compartment=args.cell_compartment,
        excluded_cell_types=EXCLUDED_CELL_TYPES,
    )
    for group, cell_ids in eligible_cells.items():
        log(f"eligible cells for {group}: {len(cell_ids)}")

    sample_size = resolve_sample_size(args.sample_size, eligible_cells)
    log(f"sample size per group: {sample_size}")

    log_stage("Sample Cells")
    sampled_cells = sample_cells_by_group(
        eligible_cells,
        sample_size=sample_size,
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
    for group in group_to_label:
        run_group_cscn(
            data_dir,
            run_name,
            group,
            matrices[group],
            max_workers=args.max_workers,
            significance_level=args.significance_level,
        )

    log_stage("Load Graphs")
    expression_df = build_expression_df(matrices, used_genes, group_to_label)
    log(f"combined expression dataframe shape: {expression_df.shape}")
    group_graphs = load_saved_group_graphs(
        data_dir=data_dir,
        dataset_name=run_name,
        groups=group_to_label.keys(),
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
