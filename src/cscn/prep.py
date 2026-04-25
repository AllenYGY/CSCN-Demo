from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from biomarker.datasets import load_gene_names

from .config import CSCNConfig, ConfigError
from .io import LoadedDataset
from .spatial_blocks import CellBlockAssignment, SpatialBlock, build_cell_block_assignments, generate_adaptive_blocks


@dataclass(frozen=True)
class PreparedAdaptiveBlockData:
    blocks: list[SpatialBlock]
    assignments: dict[str, CellBlockAssignment]
    block_gene_names: list[str]
    matrix_sum: np.ndarray
    matrix_mean: np.ndarray


@dataclass(frozen=True)
class PreparedGroupData:
    group_key: str
    cell_ids: list[str]
    matrix: np.ndarray
    raw_matrix: np.ndarray
    spatial_coords: np.ndarray | None = None
    adaptive_blocks: PreparedAdaptiveBlockData | None = None


@dataclass(frozen=True)
class PreparedRunData:
    gene_names: list[str]
    groups: dict[str, PreparedGroupData]
    sampled_metadata: pd.DataFrame


def normalize_counts(matrix: np.ndarray, target_sum: float = 1e6) -> np.ndarray:
    totals = matrix.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    return (matrix / totals) * target_sum


def preprocess_matrix(
    matrix: np.ndarray,
    *,
    normalize: bool,
    log1p: bool,
) -> np.ndarray:
    processed = np.asarray(matrix, dtype=np.float64, copy=True)
    if normalize:
        processed = normalize_counts(processed)
    if log1p:
        processed = np.log1p(processed)
    return processed


def _resolve_gene_names(dataset: LoadedDataset, config: CSCNConfig) -> list[str]:
    selection = config.preprocess.gene_selection
    if selection.gene_list_path is not None:
        requested = load_gene_names(selection.gene_list_path)
        gene_lookup = {gene: gene for gene in dataset.gene_names}
        selected = [gene_lookup[gene] for gene in requested if gene in gene_lookup]
        if not selected:
            raise ConfigError("No requested genes were found in the loaded dataset.")
        return selected

    variances = dataset.expression.var(axis=0)
    top_n = min(selection.top_n, len(variances))
    if top_n <= 0:
        raise ConfigError("No genes are available after loading the dataset.")
    ranked = variances.sort_values(ascending=False)
    return ranked.head(top_n).index.astype(str).tolist()


def _resolve_block_gene_names(
    frame: pd.DataFrame,
    gene_names: list[str],
    config: CSCNConfig,
) -> list[str]:
    top_n = min(config.run.spatial.block_gene_top_n, len(gene_names))
    ranked = frame.var(axis=0).sort_values(ascending=False)
    return ranked.head(top_n).index.astype(str).tolist()


def _resolve_group_labels(dataset: LoadedDataset, config: CSCNConfig) -> dict[str, list[str]]:
    group_key = config.input.obs_group_key
    if not group_key:
        return {"all": dataset.cell_ids}
    if group_key not in dataset.metadata.columns:
        raise ConfigError(f"Metadata is missing the configured group column: {group_key}")

    metadata = dataset.metadata.copy()
    metadata[group_key] = metadata[group_key].astype(str)
    allowed_groups = set(config.run.groups or [])
    groups: dict[str, list[str]] = {}
    for label, frame in metadata.groupby(group_key, sort=False):
        if allowed_groups and label not in allowed_groups:
            continue
        groups[str(label)] = frame.index.astype(str).tolist()

    if allowed_groups:
        missing = allowed_groups.difference(groups)
        if missing:
            raise ConfigError(f"Configured groups were not found in metadata: {sorted(missing)}")
    if not groups:
        raise ConfigError("No groups are available after applying the configured filters.")
    return groups


def _sample_cells_per_group(
    group_cells: dict[str, list[str]],
    sample_per_group: int | None,
    random_seed: int,
) -> dict[str, list[str]]:
    if sample_per_group is None:
        return group_cells

    rng = np.random.default_rng(random_seed)
    sampled: dict[str, list[str]] = {}
    for group_key, cell_ids in group_cells.items():
        if len(cell_ids) < sample_per_group:
            raise ConfigError(
                f"Group {group_key} has only {len(cell_ids)} cells, fewer than "
                f"sample_per_group={sample_per_group}."
            )
        picks = rng.choice(np.asarray(cell_ids), size=sample_per_group, replace=False)
        sampled[group_key] = [str(item) for item in picks.tolist()]
    return sampled


def _build_adaptive_blocks(
    *,
    cell_ids: list[str],
    raw_frame: pd.DataFrame,
    processed_frame: pd.DataFrame,
    spatial_coords: np.ndarray,
    gene_names: list[str],
    config: CSCNConfig,
) -> PreparedAdaptiveBlockData:
    density = config.run.spatial.density_clustering or {}
    blocks = generate_adaptive_blocks(
        cell_ids=cell_ids,
        coords=np.asarray(spatial_coords[:, :2], dtype=np.float64),
        eps=float(density.get("eps", 1.5)),
        min_samples=int(density.get("min_samples", 2)),
        min_cells_per_block=config.run.spatial.min_cells_per_block,
        overlap_min=config.run.spatial.block_overlap_min,
    )
    assignments = build_cell_block_assignments(
        cell_ids=cell_ids,
        coords=np.asarray(spatial_coords[:, :2], dtype=np.float64),
        blocks=blocks,
        halo_neighbor_blocks=config.run.spatial.halo_neighbor_blocks,
    )
    block_gene_names = _resolve_block_gene_names(processed_frame, gene_names, config)
    block_sum_rows: list[np.ndarray] = []
    block_mean_rows: list[np.ndarray] = []
    for block in blocks:
        member_ids = list(block.cell_ids)
        raw_block = raw_frame.loc[member_ids, block_gene_names].to_numpy(dtype=np.float64, copy=True)
        processed_block = processed_frame.loc[member_ids, block_gene_names].to_numpy(
            dtype=np.float64,
            copy=True,
        )
        sum_row = preprocess_matrix(
            raw_block.sum(axis=0, keepdims=True),
            normalize=config.preprocess.normalize,
            log1p=config.preprocess.log1p,
        )[0]
        mean_row = processed_block.mean(axis=0)
        block_sum_rows.append(sum_row.astype(np.float32))
        block_mean_rows.append(np.asarray(mean_row, dtype=np.float32))
    if not block_sum_rows:
        empty = np.zeros((0, len(block_gene_names)), dtype=np.float32)
        return PreparedAdaptiveBlockData(
            blocks=blocks,
            assignments=assignments,
            block_gene_names=block_gene_names,
            matrix_sum=empty,
            matrix_mean=empty.copy(),
        )
    return PreparedAdaptiveBlockData(
        blocks=blocks,
        assignments=assignments,
        block_gene_names=block_gene_names,
        matrix_sum=np.vstack(block_sum_rows).astype(np.float32),
        matrix_mean=np.vstack(block_mean_rows).astype(np.float32),
    )


def prepare_run_inputs(dataset: LoadedDataset, config: CSCNConfig) -> PreparedRunData:
    gene_names = _resolve_gene_names(dataset, config)
    grouped_cells = _resolve_group_labels(dataset, config)
    sampled_cells = _sample_cells_per_group(
        group_cells=grouped_cells,
        sample_per_group=config.preprocess.sample_per_group,
        random_seed=config.preprocess.random_seed,
    )

    groups: dict[str, PreparedGroupData] = {}
    sampled_metadata_frames: list[pd.DataFrame] = []
    spatial_columns = [
        key
        for key in (
            config.input.spatial_x_key,
            config.input.spatial_y_key,
            config.input.spatial_z_key,
        )
        if key
    ]
    if spatial_columns:
        missing_columns = [
            column for column in spatial_columns if column not in dataset.metadata.columns
        ]
        if missing_columns:
            raise ConfigError(
                f"Metadata is missing the configured spatial columns: {missing_columns}"
            )
    for group_key, cell_ids in sampled_cells.items():
        raw_frame = dataset.expression.loc[cell_ids, gene_names].copy()
        raw_matrix = raw_frame.to_numpy(dtype=np.float64, copy=True)
        matrix = preprocess_matrix(
            raw_matrix,
            normalize=config.preprocess.normalize,
            log1p=config.preprocess.log1p,
        )
        processed_frame = pd.DataFrame(matrix, index=cell_ids, columns=gene_names)
        spatial_coords = None
        adaptive_blocks = None
        if config.input.spatial_x_key and config.input.spatial_y_key:
            coords_frame = dataset.metadata.loc[cell_ids, spatial_columns].copy()
            coords_frame = coords_frame.apply(pd.to_numeric, errors="coerce")
            if coords_frame.isnull().any().any():
                raise ConfigError(
                    f"Spatial coordinates contain missing or non-numeric values for group "
                    f"{group_key}."
                )
            spatial_coords = coords_frame.to_numpy(dtype=np.float32, copy=True)
            if config.run.spatial.enabled and config.run.spatial.strategy == "adaptive_block_prior":
                adaptive_blocks = _build_adaptive_blocks(
                    cell_ids=[str(cell_id) for cell_id in cell_ids],
                    raw_frame=raw_frame,
                    processed_frame=processed_frame,
                    spatial_coords=spatial_coords,
                    gene_names=gene_names,
                    config=config,
                )
        groups[group_key] = PreparedGroupData(
            group_key=group_key,
            cell_ids=[str(cell_id) for cell_id in cell_ids],
            matrix=matrix.astype(np.float32),
            raw_matrix=raw_matrix.astype(np.float32),
            spatial_coords=spatial_coords,
            adaptive_blocks=adaptive_blocks,
        )

        group_metadata = dataset.metadata.loc[cell_ids].copy()
        group_metadata.insert(0, "cscn_group", group_key)
        sampled_metadata_frames.append(group_metadata)

    sampled_metadata = pd.concat(sampled_metadata_frames, axis=0)
    if "cell_id" in sampled_metadata.columns:
        sampled_metadata = sampled_metadata.drop(columns=["cell_id"])
    sampled_metadata.insert(0, "cell_id", sampled_metadata.index.astype(str))
    return PreparedRunData(
        gene_names=gene_names,
        groups=groups,
        sampled_metadata=sampled_metadata.reset_index(drop=True),
    )
