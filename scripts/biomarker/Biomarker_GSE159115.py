from __future__ import annotations

import argparse
import gzip
import os
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
import h5py

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
TUMOR_ANNOS = {"Tumor"}
NORMAL_ANNOS = {"PT-B", "PT-C"}


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


def read_annotation_table(annotation_path: Path, disease: str, allowed_annos):
    import csv
    from collections import Counter

    rows = []
    summary = {
        "total_cells": 0,
        "patient_counts": Counter(),
        "sample_counts": Counter(),
        "anno_counts": Counter(),
    }

    with gzip.open(annotation_path, "rt", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"cell", "sample", "anno", "patient", "doublet"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Missing required columns in {annotation_path}: {sorted(missing)}")

        for row in reader:
            if row["anno"] not in allowed_annos:
                continue
            if row["doublet"] != "FALSE":
                continue
            item = {
                "cell": row["cell"],
                "sample": row["sample"],
                "patient": row["patient"],
                "anno": row["anno"],
                "disease": disease,
                "group": disease,
            }
            rows.append(item)
            summary["total_cells"] += 1
            summary["patient_counts"][row["patient"]] += 1
            summary["sample_counts"][row["sample"]] += 1
            summary["anno_counts"][row["anno"]] += 1

    return rows, summary


def select_paired_grouped_cells(ccrcc_annotation_path: Path, normal_annotation_path: Path):
    tumor_rows, tumor_summary = read_annotation_table(ccrcc_annotation_path, "tumor", TUMOR_ANNOS)
    normal_rows, normal_summary = read_annotation_table(normal_annotation_path, "normal", NORMAL_ANNOS)

    tumor_patients = {row["patient"] for row in tumor_rows}
    normal_patients = {row["patient"] for row in normal_rows}
    paired_patients = sorted(tumor_patients & normal_patients)

    tumor_rows = [row for row in tumor_rows if row["patient"] in paired_patients]
    normal_rows = [row for row in normal_rows if row["patient"] in paired_patients]

    grouped_cells = {
        "tumor": [row["cell"] for row in tumor_rows],
        "normal": [row["cell"] for row in normal_rows],
    }
    from collections import Counter

    filtered_tumor_samples = Counter(row["sample"] for row in tumor_rows)
    filtered_tumor_annos = Counter(row["anno"] for row in tumor_rows)
    filtered_tumor_patients = Counter(row["patient"] for row in tumor_rows)
    filtered_normal_samples = Counter(row["sample"] for row in normal_rows)
    filtered_normal_annos = Counter(row["anno"] for row in normal_rows)
    filtered_normal_patients = Counter(row["patient"] for row in normal_rows)
    summary_by_group = {
        "tumor": {
            "total_cells": len(tumor_rows),
            "patient_counts": filtered_tumor_patients,
            "sample_counts": filtered_tumor_samples,
            "anno_counts": filtered_tumor_annos,
            "paired_patients": paired_patients,
        },
        "normal": {
            "total_cells": len(normal_rows),
            "patient_counts": filtered_normal_patients,
            "sample_counts": filtered_normal_samples,
            "anno_counts": filtered_normal_annos,
            "paired_patients": paired_patients,
        },
    }
    return grouped_cells, summary_by_group


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
        anno_counts = summary["anno_counts"]
        if anno_counts:
            preview = ", ".join(
                f"{name}={count}"
                for name, count in sorted(anno_counts.items(), key=lambda item: (-item[1], item[0]))
            )
            log(f"{group} anno counts: {preview}")
        paired = summary.get("paired_patients")
        if paired:
            log(f"{group} paired patients: {', '.join(paired)}")


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


def extract_h5_cache(raw_tar_path: Path):
    out_dir = raw_tar_path.parent / ".gse159115_h5_cache_py"
    out_dir.mkdir(parents=True, exist_ok=True)
    h5_files = sorted(out_dir.glob("*.h5"))
    if not h5_files:
        with tarfile.open(raw_tar_path) as tf:
            tf.extractall(out_dir)
        h5_files = sorted(out_dir.glob("*.h5"))
    if not h5_files:
        raise FileNotFoundError(f"No .h5 files found after extracting {raw_tar_path}")
    sample_to_h5 = {}
    for path in h5_files:
        sample = path.name.split("_filtered_gene_bc_matrices_h5.h5")[0].split("_", 1)[1]
        sample_to_h5[sample] = path
    return sample_to_h5


def split_sample_barcode(cell_id: str):
    parts = cell_id.split("_", 2)
    if len(parts) != 3:
        raise ValueError(f"Unexpected GSE159115 cell id format: {cell_id}")
    return f"{parts[0]}_{parts[1]}", parts[2]


def read_10x_h5_old(h5_path: Path):
    with h5py.File(h5_path, "r") as h5:
        root_name = next(iter(h5.keys()))
        grp = h5[root_name]
        gene_names = [x.decode() if isinstance(x, bytes) else str(x) for x in grp["gene_names"][()]]
        barcodes = [x.decode() if isinstance(x, bytes) else str(x) for x in grp["barcodes"][()]]
        indices = grp["indices"][()]
        indptr = grp["indptr"][()]
        data = grp["data"][()]
        shape = tuple(int(x) for x in grp["shape"][()])
    mat = csc_matrix((data, indices, indptr), shape=shape)
    return gene_names, barcodes, mat


def extract_sampled_expression_from_h5(raw_tar_path: Path, sampled_cells, top_genes):
    sample_to_h5 = extract_h5_cache(raw_tar_path)
    top_genes = list(top_genes)
    matrices = {group: None for group in sampled_cells}
    used_genes = None

    for group, cell_ids in sampled_cells.items():
        by_sample = {}
        for cell_id in cell_ids:
            sample_id, barcode = split_sample_barcode(cell_id)
            by_sample.setdefault(sample_id, []).append((cell_id, barcode))

        group_blocks = []
        for sample_id, pairs in by_sample.items():
            h5_path = sample_to_h5.get(sample_id)
            if h5_path is None:
                raise FileNotFoundError(f"No H5 file found for sample {sample_id}")
            gene_names, barcodes, mat = read_10x_h5_old(h5_path)
            barcode_to_col = {barcode: idx for idx, barcode in enumerate(barcodes)}

            selected_cols = []
            ordered_cells = []
            for cell_id, barcode in pairs:
                col_idx = barcode_to_col.get(barcode)
                if col_idx is None:
                    raise ValueError(f"Missing cell barcode {barcode} in sample {sample_id}")
                selected_cols.append(col_idx)
                ordered_cells.append(cell_id)

            selected_cols = np.array(selected_cols, dtype=int)
            totals = np.asarray(mat[:, selected_cols].sum(axis=0)).ravel().astype(np.float64)
            totals[totals == 0] = 1.0

            gene_to_rows = {}
            for idx, gene in enumerate(gene_names):
                if gene in top_genes:
                    gene_to_rows.setdefault(gene, []).append(idx)

            if used_genes is None:
                used_genes = [gene for gene in top_genes if gene in gene_to_rows]
                if not used_genes:
                    raise ValueError("None of the requested top genes were found in the GSE159115 h5 files.")

            block = np.zeros((len(selected_cols), len(used_genes)), dtype=np.float64)
            for gene_pos, gene in enumerate(used_genes):
                rows = gene_to_rows.get(gene)
                if not rows:
                    continue
                vec = np.asarray(mat[rows, :][:, selected_cols].sum(axis=0)).ravel()
                block[:, gene_pos] = vec

            block = np.log1p((block / totals[:, None]) * 1e6)
            group_blocks.append((ordered_cells, block))

        ordered = []
        arrays = []
        for ordered_cells, block in group_blocks:
            ordered.extend(ordered_cells)
            arrays.append(block)
        matrices[group] = (ordered, np.vstack(arrays))

    sampled_cell_order = {group: cells for group, (cells, _) in matrices.items()}
    matrix_values = {group: arr for group, (_, arr) in matrices.items()}
    return sampled_cell_order, matrix_values, used_genes


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
    raw_tar_path = data_dir / "GSE159115_RAW.tar"
    ccrcc_annotation_path = data_dir / "GSE159115_ccRCC_anno.csv.gz"
    normal_annotation_path = data_dir / "GSE159115_normal_anno.csv.gz"
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
    validate_required_file(raw_tar_path, "RAW tar")
    validate_required_file(ccrcc_annotation_path, "ccRCC annotation")
    validate_required_file(normal_annotation_path, "normal annotation")
    validate_required_file(gene_list_path, "top-gene list")

    log_stage("Load Inputs")
    top_genes = load_gene_names(gene_list_path)
    grouped_cells, summary_by_group = select_paired_grouped_cells(
        ccrcc_annotation_path=ccrcc_annotation_path,
        normal_annotation_path=normal_annotation_path,
    )
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
    sampled_cells, matrices, used_genes = extract_sampled_expression_from_h5(
        raw_tar_path=raw_tar_path,
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
