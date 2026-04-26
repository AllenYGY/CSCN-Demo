from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


DATA_SET = "SCP259"
RUN_SLUG = "inflamed_vs_healthy_crypt_prolif_epi"
RUN_NAME = f"{DATA_SET}_{RUN_SLUG}"
GROUP_TO_LABEL = {"healthy": 0, "inflamed": 1}
DEFAULT_SAMPLE_SIZE = 4000
DEFAULT_RANDOM_SEED = 42
DEFAULT_MAX_WORKERS = min(8, os.cpu_count() or 1)
DEFAULT_USE_BITMAP = True
CRYPT_PROLIF_CLUSTERS = {
    "Stem",
    "Cycling TA",
    "TA 1",
    "TA 2",
    "Enterocyte Progenitors",
    "Secretory TA",
}


def log(message):
    print(f"[{DATA_SET}] {message}", flush=True)


def log_stage(title):
    print(flush=True)
    print(f"=== {title} ===", flush=True)


def resolve_default_data_dir():
    candidates = [
        REPO_ROOT / "data" / "scp259",
        REPO_ROOT / "data" / "SCP259",
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
    return output_dir / "deseq2_inflamed_vs_healthy_crypt_prolif_epi_top150_genes.csv"


def read_scp259_metadata_rows(metadata_path: Path):
    import csv
    from collections import Counter

    grouped = {"healthy": [], "inflamed": []}
    summary = {
        "healthy": {"total_cells": 0, "subject_counts": Counter(), "sample_counts": Counter(), "cluster_counts": Counter()},
        "inflamed": {"total_cells": 0, "subject_counts": Counter(), "sample_counts": Counter(), "cluster_counts": Counter()},
    }

    with open(metadata_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"NAME", "Cluster", "Subject", "Health", "Location", "Sample"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Missing required metadata columns in {metadata_path}: {sorted(missing)}")

        for row in reader:
            if row["NAME"] == "TYPE":
                continue
            if row["Location"] != "Epi":
                continue
            if row["Cluster"] not in CRYPT_PROLIF_CLUSTERS:
                continue
            health = row["Health"]
            if health == "Healthy":
                group = "healthy"
            elif health == "Inflamed":
                group = "inflamed"
            else:
                continue

            cell_id = row["NAME"]
            grouped[group].append(cell_id)
            info = summary[group]
            info["total_cells"] += 1
            info["subject_counts"][row["Subject"]] += 1
            info["sample_counts"][row["Sample"]] += 1
            info["cluster_counts"][row["Cluster"]] += 1

    return grouped, summary


def log_group_summary(summary_by_group):
    for group, summary in summary_by_group.items():
        log(f"eligible cells for {group}: {summary['total_cells']}")
        sample_counts = summary["sample_counts"]
        if sample_counts:
            preview = ", ".join(
                f"{name}={count}"
                for name, count in sorted(sample_counts.items(), key=lambda item: (-item[1], item[0]))[:10]
            )
            log(f"{group} sample counts: {preview}")
        cluster_counts = summary["cluster_counts"]
        if cluster_counts:
            preview = ", ".join(
                f"{name}={count}"
                for name, count in sorted(cluster_counts.items(), key=lambda item: (-item[1], item[0]))
            )
            log(f"{group} cluster counts: {preview}")


def sample_cells_by_group(cells_by_group, sample_size, random_seed):
    rng = np.random.default_rng(random_seed)
    sampled_cells = {}
    for group, cell_ids in cells_by_group.items():
        if len(cell_ids) < sample_size:
            raise ValueError(
                f"Group {group} has only {len(cell_ids)} eligible cells, fewer than sample size {sample_size}."
            )
        sampled = rng.choice(np.array(cell_ids), size=sample_size, replace=False)
        sampled_cells[group] = sampled.tolist()
    return sampled_cells


def load_scp259_barcodes(barcodes_path: Path):
    barcodes = []
    with open(barcodes_path) as handle:
        for line in handle:
            value = line.strip()
            if not value or value == "TYPE":
                continue
            barcodes.append(value)
    return barcodes


def load_scp259_genes(genes_path: Path):
    genes = []
    with open(genes_path) as handle:
        for line in handle:
            value = line.strip()
            if value:
                genes.append(value)
    return genes


def extract_scp259_sampled_expression(matrix_path: Path, barcodes_path: Path, genes_path: Path, sampled_cells, top_genes):
    barcodes = load_scp259_barcodes(barcodes_path)
    cell_to_col = {cell_id: idx + 1 for idx, cell_id in enumerate(barcodes)}

    group_col_lookup = {}
    for group, cell_ids in sampled_cells.items():
        missing = [cell_id for cell_id in cell_ids if cell_id not in cell_to_col]
        if missing:
            raise ValueError(
                f"Missing sampled cells from barcodes file for group {group}: {missing[:5]}"
            )
        group_col_lookup[group] = [cell_to_col[cell_id] for cell_id in cell_ids]

    genes = load_scp259_genes(genes_path)
    gene_to_row = {}
    for idx, gene in enumerate(genes, start=1):
        if gene not in gene_to_row:
            gene_to_row[gene] = idx

    used_genes = [gene for gene in top_genes if gene in gene_to_row]
    if not used_genes:
        raise ValueError("None of the requested top genes were found in the SCP259 Epi matrix.")

    row_target_map = {gene_to_row[gene]: pos for pos, gene in enumerate(used_genes)}

    selected_col_map = {}
    matrices = {}
    totals = {}
    for group, column_ids in group_col_lookup.items():
        matrices[group] = np.zeros((len(column_ids), len(used_genes)), dtype=np.float64)
        totals[group] = np.zeros(len(column_ids), dtype=np.float64)
        for row_idx, col_idx in enumerate(column_ids):
            selected_col_map[col_idx] = (group, row_idx)

    log("streaming SCP259 Epi matrix to extract sampled cells and top genes; this can take a while")
    with open(matrix_path) as handle:
        first = handle.readline().strip()
        if not first.startswith("%%MatrixMarket"):
            raise ValueError(f"Unexpected Matrix Market header in {matrix_path}: {first}")

        dims_line = None
        for line in handle:
            if not line.startswith("%"):
                dims_line = line.strip()
                break
        if dims_line is None:
            raise ValueError(f"Missing dimension line in {matrix_path}")
        n_rows, n_cols, n_entries = map(int, dims_line.split())
        log(f"SCP259 Epi matrix dims: genes={n_rows}, cells={n_cols}, nnz={n_entries}")

        for line_number, line in enumerate(handle, start=1):
            parts = line.split()
            if len(parts) != 3:
                continue
            row_idx = int(parts[0])
            col_idx = int(parts[1])
            value = float(parts[2])

            target = selected_col_map.get(col_idx)
            if target is None:
                continue
            group, cell_pos = target
            totals[group][cell_pos] += value

            gene_pos = row_target_map.get(row_idx)
            if gene_pos is not None:
                matrices[group][cell_pos, gene_pos] += value

            if line_number % 20000000 == 0:
                log(f"processed {line_number} matrix entries")

    for group in matrices:
        group_totals = totals[group]
        group_totals[group_totals == 0] = 1.0
        matrices[group] = np.log1p((matrices[group] / group_totals[:, None]) * 1e6)

    return matrices, used_genes


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run CSCN biomarker analysis for SCP259 using crypt/proliferative "
            "epithelial cells: Healthy vs Inflamed."
        )
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=resolve_default_data_dir(),
        help="Dataset directory containing SCP259 files.",
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


def main():
    args = parse_args()

    from biomarker.causal import run_causal_analysis
    from biomarker.datasets import (
        build_expression_df,
        load_gene_names,
        load_saved_group_graphs,
        save_prepared_inputs,
    )
    from biomarker.graph_utils import (
        get_global_graph,
        identify_biomarkers_from_group_graphs,
        map_node_id_to_gene,
    )
    from biomarker.cscn import CSCN

    data_dir = args.data_dir.resolve()
    output_dir = data_dir / "output_deseq"
    metadata_path = data_dir / "metadata" / "all.meta2.txt"
    expression_dir = data_dir / "expression" / "5cdc540d328cee7a2efc2348"
    matrix_path = expression_dir / "gene_sorted-Epi.matrix.mtx"
    barcodes_path = expression_dir / "Epi.barcodes2.tsv"
    genes_path = expression_dir / "Epi.genes.tsv"
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
    validate_required_file(metadata_path, "metadata table")
    validate_required_file(matrix_path, "Epi matrix")
    validate_required_file(barcodes_path, "Epi barcodes")
    validate_required_file(genes_path, "Epi genes")
    validate_required_file(gene_list_path, "top-gene list")

    log_stage("Load Inputs")
    top_genes = load_gene_names(gene_list_path)
    grouped_cells, summary_by_group = read_scp259_metadata_rows(metadata_path)
    log_group_summary(summary_by_group)

    log_stage("Sample Cells")
    sampled_cells = sample_cells_by_group(
        grouped_cells,
        sample_size=args.sample_size,
        random_seed=args.random_seed,
    )
    for group, cell_ids in sampled_cells.items():
        log(f"sampled cells for {group}: {len(cell_ids)}")
        log(f"first 3 sampled {group} cells: {cell_ids[:3]}")

    log_stage("Extract Expression")
    matrices, used_genes = extract_scp259_sampled_expression(
        matrix_path=matrix_path,
        barcodes_path=barcodes_path,
        genes_path=genes_path,
        sampled_cells=sampled_cells,
        top_genes=top_genes,
    )

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
        used_genes_filename=f"{run_name}_top{len(used_genes)}_genes_used.csv",
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
