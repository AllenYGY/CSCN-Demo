import csv
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def unique_preserve_order(values):
    seen = set()
    ordered = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def load_gene_list_df(gene_list_path, gene_column="gene", fallback_columns=("Unnamed: 0",)):
    gene_info_df = pd.read_csv(gene_list_path)
    if gene_column in gene_info_df.columns:
        return gene_info_df[[gene_column]].rename(columns={gene_column: "gene"})
    for fallback in fallback_columns:
        if fallback in gene_info_df.columns:
            return gene_info_df[[fallback]].rename(columns={fallback: "gene"})
    raise ValueError(f"Missing gene column in {gene_list_path}")


def load_gene_names(
    gene_list_path,
    gene_column="gene",
    fallback_columns=("Unnamed: 0",),
    deduplicate=True,
):
    gene_df = load_gene_list_df(
        gene_list_path,
        gene_column=gene_column,
        fallback_columns=fallback_columns,
    )
    genes = gene_df["gene"].dropna().astype(str).tolist()
    if deduplicate:
        genes = unique_preserve_order(genes)
    if not genes:
        raise ValueError("Top-gene list is empty.")
    return genes


def build_expression_df(matrices, gene_names, group_to_label):
    frames = []
    for group, label in group_to_label.items():
        frame = pd.DataFrame(matrices[group], columns=gene_names)
        frame["DISEASE"] = label
        frames.append(frame)
    return pd.concat(frames, axis=0, ignore_index=True)


def load_group_npy_expression_inputs(gene_list_path, group_npy_paths, group_to_label):
    gene_list_df = load_gene_list_df(gene_list_path)
    gene_names = gene_list_df["gene"].dropna().astype(str).tolist()
    matrices = {group: np.load(path) for group, path in group_npy_paths.items()}
    expression_df = build_expression_df(matrices, gene_names, group_to_label)
    id2gene = {idx: gene for idx, gene in enumerate(gene_names)}
    return gene_list_df, gene_names, matrices, expression_df, id2gene


def load_saved_group_graphs(
    data_dir,
    dataset_name,
    groups,
    gene_names,
    cscn_cls,
    map_node_id_to_gene_fn,
    get_global_graph_fn,
):
    id2gene = {idx: gene for idx, gene in enumerate(gene_names)}
    group_graphs = {}
    for group in groups:
        cscn = cscn_cls.load_from_file(data_dir / f"{dataset_name}_{group}_cscn")
        dags = cscn.load_all_dags()
        dags = map_node_id_to_gene_fn(dags, id2gene)
        group_graphs[group] = get_global_graph_fn(dags)
    return group_graphs


def load_gse138852_eligible_cells(covariates_path, group_to_label, excluded_cell_types):
    cells_by_group = {group: [] for group in group_to_label}
    with gzip.open(covariates_path, "rt", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            cell_id = row[""]
            condition = row["oupSample.batchCond"]
            cell_type = row["oupSample.cellType"]
            if condition not in group_to_label:
                continue
            if cell_type in excluded_cell_types:
                continue
            cells_by_group[condition].append(cell_id)
    return cells_by_group


def load_gse131907_grouped_cells(annotation_path, group_filters):
    annotation_df = pd.read_csv(annotation_path, sep="\t")
    required_columns = {"Index", "Sample_Origin", "Cell_type", "Cell_subtype"}
    missing_columns = required_columns - set(annotation_df.columns)
    if missing_columns:
        raise ValueError(
            f"Missing required annotation columns in {annotation_path}: {sorted(missing_columns)}"
        )

    cells_by_group = {}
    summary_by_group = {}
    for group, filters in group_filters.items():
        mask = pd.Series(True, index=annotation_df.index)

        sample_origins = filters.get("sample_origins")
        if sample_origins:
            mask &= annotation_df["Sample_Origin"].isin(sample_origins)

        cell_types = filters.get("cell_types")
        if cell_types:
            mask &= annotation_df["Cell_type"].isin(cell_types)

        cell_subtypes = filters.get("cell_subtypes")
        if cell_subtypes:
            mask &= annotation_df["Cell_subtype"].isin(cell_subtypes)

        selected_df = annotation_df.loc[mask].copy()
        cells_by_group[group] = selected_df["Index"].astype(str).tolist()
        summary_by_group[group] = {
            "total_cells": int(selected_df.shape[0]),
            "sample_origin_counts": selected_df["Sample_Origin"].value_counts().to_dict(),
            "cell_type_counts": selected_df["Cell_type"].value_counts().to_dict(),
            "cell_subtype_counts": selected_df["Cell_subtype"].value_counts().to_dict(),
        }

    return cells_by_group, summary_by_group


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


def load_gse121893_eligible_cells(
    covariates_path,
    group_to_label,
    allowed_regions=None,
    cell_compartment="all",
    excluded_cell_types=None,
):
    normalized_regions = None
    if allowed_regions:
        normalized_regions = {str(region).upper() for region in allowed_regions}
        if "ALL" in normalized_regions:
            normalized_regions = None

    normalized_compartment = None
    if cell_compartment and str(cell_compartment).lower() != "all":
        normalized_compartment = str(cell_compartment).upper()

    excluded_cell_types = set(excluded_cell_types or ())
    cells_by_group = {group: [] for group in group_to_label}

    with gzip.open(covariates_path, "rt", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            disease = row["disease"]
            if disease not in group_to_label:
                continue

            region = row["region"].upper()
            if normalized_regions and region not in normalized_regions:
                continue

            compartment = row["group"].upper()
            if normalized_compartment and compartment != normalized_compartment:
                continue

            cell_type = row["cell_type"]
            if cell_type in excluded_cell_types:
                continue

            cells_by_group[disease].append(row["cell_id"])

    return cells_by_group


def read_expression_for_sampled_cells(
    counts_path,
    sampled_cells,
    top_genes,
    delimiter=",",
):
    max_csv_field_size = sys.maxsize
    while True:
        try:
            csv.field_size_limit(max_csv_field_size)
            break
        except OverflowError:
            max_csv_field_size //= 10

    group_positions = {}
    gene_to_group_values = {group: {} for group in sampled_cells}
    gene_set = set(top_genes)

    with gzip.open(counts_path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader)
        cell_to_position = {cell_id: idx for idx, cell_id in enumerate(header[1:])}

        for group, cell_ids in sampled_cells.items():
            missing = [cell_id for cell_id in cell_ids if cell_id not in cell_to_position]
            if missing:
                raise ValueError(
                    f"Missing sampled cells from counts matrix for group {group}: {missing[:5]}"
                )
            group_positions[group] = [cell_to_position[cell_id] for cell_id in cell_ids]

        for row in reader:
            gene = row[0]
            if gene not in gene_set:
                continue
            for group, positions in group_positions.items():
                values = np.fromiter(
                    (float(row[position + 1]) for position in positions),
                    dtype=np.float64,
                    count=len(positions),
                )
                gene_to_group_values[group][gene] = values

    used_genes = [
        gene
        for gene in top_genes
        if all(gene in gene_to_group_values[group] for group in sampled_cells)
    ]
    if not used_genes:
        raise ValueError("None of the requested genes were found in the counts matrix.")

    matrices = {}
    for group in sampled_cells:
        matrices[group] = np.vstack(
            [gene_to_group_values[group][gene] for gene in used_genes]
        ).T
    return matrices, used_genes


def normalize_log1p(matrix, target_sum=1e6):
    totals = matrix.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    normalized = (matrix / totals) * target_sum
    return np.log1p(normalized)


def save_prepared_inputs(
    output_dir,
    dataset_name,
    sampled_cells,
    matrices,
    used_genes,
    used_genes_filename=None,
):
    output_dir.mkdir(parents=True, exist_ok=True)

    gene_df = pd.DataFrame({"gene": used_genes})
    if used_genes_filename is None:
        used_genes_filename = f"{dataset_name}_top150_genes_used.csv"
    gene_df.to_csv(output_dir / used_genes_filename, index=False)

    for group, matrix in matrices.items():
        np.save(output_dir / f"{dataset_name}_{group}.npy", matrix.astype(np.float32))
        pd.DataFrame({"cell_id": sampled_cells[group]}).to_csv(
            output_dir / f"{dataset_name}_{group}_sampled_cells.csv",
            index=False,
        )
