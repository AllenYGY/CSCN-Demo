from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from biomarker.datasets import load_gene_names

from .config import CSCNConfig, ConfigError
from .io import LoadedDataset


@dataclass(frozen=True)
class PreparedGroupData:
    group_key: str
    cell_ids: list[str]
    matrix: np.ndarray


@dataclass(frozen=True)
class PreparedRunData:
    gene_names: list[str]
    groups: dict[str, PreparedGroupData]
    sampled_metadata: pd.DataFrame


def _normalize_counts(matrix: np.ndarray, target_sum: float = 1e6) -> np.ndarray:
    totals = matrix.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    return (matrix / totals) * target_sum


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
                f"Group {group_key} has only {len(cell_ids)} cells, fewer than sample_per_group={sample_per_group}."
            )
        picks = rng.choice(np.asarray(cell_ids), size=sample_per_group, replace=False)
        sampled[group_key] = [str(item) for item in picks.tolist()]
    return sampled


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
    for group_key, cell_ids in sampled_cells.items():
        frame = dataset.expression.loc[cell_ids, gene_names]
        matrix = frame.to_numpy(dtype=np.float64, copy=True)
        if config.preprocess.normalize:
            matrix = _normalize_counts(matrix)
        if config.preprocess.log1p:
            matrix = np.log1p(matrix)
        groups[group_key] = PreparedGroupData(
            group_key=group_key,
            cell_ids=[str(cell_id) for cell_id in cell_ids],
            matrix=matrix.astype(np.float32),
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
