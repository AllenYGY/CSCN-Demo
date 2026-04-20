from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .config import CSCNConfig, ConfigError


@dataclass(frozen=True)
class LoadedDataset:
    expression: pd.DataFrame
    metadata: pd.DataFrame
    cell_ids: list[str]
    gene_names: list[str]


def _infer_delimiter(path: Path, explicit: str | None) -> str | None:
    if explicit:
        return explicit
    if path.suffix.lower() in {".tsv", ".txt"} or path.name.endswith(".tsv.gz"):
        return "\t"
    return ","


def _read_table(path: Path, delimiter: str | None) -> pd.DataFrame:
    if not path.is_file():
        raise ConfigError(f"Missing table file: {path}")
    sep = _infer_delimiter(path, delimiter)
    return pd.read_csv(path, sep=sep)


def _deduplicate_names(names: list[str]) -> list[str]:
    counts: dict[str, int] = {}
    deduped: list[str] = []
    for name in names:
        base = str(name)
        seen = counts.get(base, 0)
        counts[base] = seen + 1
        if seen == 0:
            deduped.append(base)
        else:
            deduped.append(f"{base}__{seen}")
    return deduped


def _coerce_numeric_frame(frame: pd.DataFrame) -> pd.DataFrame:
    numeric = frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return numeric.astype(np.float64)


def _align_metadata(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    metadata_cell_id_column: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if (
        metadata_cell_id_column == ""
        and metadata_cell_id_column not in metadata.columns
    ):
        unnamed_columns = [
            column
            for column in metadata.columns
            if str(column).startswith("Unnamed:")
        ]
        if unnamed_columns:
            metadata_cell_id_column = str(unnamed_columns[0])
    if metadata_cell_id_column not in metadata.columns:
        raise ConfigError(
            f"Metadata is missing the required cell id column: {metadata_cell_id_column}"
        )
    metadata = metadata.copy()
    metadata[metadata_cell_id_column] = metadata[metadata_cell_id_column].astype(str)
    metadata = metadata.drop_duplicates(subset=[metadata_cell_id_column]).set_index(
        metadata_cell_id_column
    )
    expression = expression.copy()
    expression.index = expression.index.astype(str)
    shared = [cell_id for cell_id in expression.index if cell_id in metadata.index]
    if not shared:
        raise ConfigError("Expression matrix and metadata do not share any cell ids.")
    return expression.loc[shared], metadata.loc[shared]


def _load_h5ad(config: CSCNConfig) -> LoadedDataset:
    try:
        import anndata as ad
    except ModuleNotFoundError as exc:
        raise ConfigError(
            "h5ad input requires `anndata`. Install the project dependencies first."
        ) from exc

    path = config.input.path
    assert path is not None
    adata = ad.read_h5ad(path)
    matrix = adata.layers[config.input.layer] if config.input.layer else adata.X
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix, dtype=np.float64)

    if config.input.gene_key and config.input.gene_key in adata.var.columns:
        gene_names = adata.var[config.input.gene_key].astype(str).tolist()
    else:
        gene_names = adata.var_names.astype(str).tolist()
    gene_names = _deduplicate_names(gene_names)

    if config.input.obs_cell_id_key and config.input.obs_cell_id_key in adata.obs.columns:
        cell_ids = adata.obs[config.input.obs_cell_id_key].astype(str).tolist()
    else:
        cell_ids = adata.obs_names.astype(str).tolist()

    expression = pd.DataFrame(matrix, index=cell_ids, columns=gene_names)
    metadata = adata.obs.copy()
    metadata.index = cell_ids
    return LoadedDataset(
        expression=expression,
        metadata=metadata,
        cell_ids=cell_ids,
        gene_names=gene_names,
    )


def _load_tables(config: CSCNConfig) -> LoadedDataset:
    expr_path = config.input.expr_path
    metadata_path = config.input.metadata_path
    assert expr_path is not None
    assert metadata_path is not None

    expr = _read_table(expr_path, config.input.expr_delimiter)
    metadata = _read_table(metadata_path, config.input.metadata_delimiter)

    if config.input.expr_orientation == "genes_by_cells":
        gene_index_column = config.input.gene_key or expr.columns[0]
        if gene_index_column not in expr.columns:
            raise ConfigError(
                f"Expression matrix is missing the gene column: {gene_index_column}"
            )
        expr = expr.set_index(gene_index_column)
        expr.columns = expr.columns.astype(str)
        expr = _coerce_numeric_frame(expr).T
        expr.index = expr.index.astype(str)
        expr.columns = _deduplicate_names(expr.columns.astype(str).tolist())
    else:
        if config.input.expr_cell_id_column and config.input.expr_cell_id_column in expr.columns:
            expr = expr.set_index(config.input.expr_cell_id_column)
        else:
            first_column = expr.columns[0]
            metadata_ids = set(metadata[config.input.metadata_cell_id_column].astype(str))
            overlap_ratio = expr[first_column].astype(str).isin(metadata_ids).mean()
            if overlap_ratio > 0.5:
                expr = expr.set_index(first_column)
        expr = _coerce_numeric_frame(expr)
        expr.index = expr.index.astype(str)
        expr.columns = _deduplicate_names(expr.columns.astype(str).tolist())

    expr, metadata = _align_metadata(
        expression=expr,
        metadata=metadata,
        metadata_cell_id_column=config.input.metadata_cell_id_column,
    )
    return LoadedDataset(
        expression=expr,
        metadata=metadata,
        cell_ids=expr.index.astype(str).tolist(),
        gene_names=expr.columns.astype(str).tolist(),
    )


def load_dataset(config: CSCNConfig) -> LoadedDataset:
    if config.input.format == "h5ad":
        return _load_h5ad(config)
    if config.input.format == "tables":
        return _load_tables(config)
    raise ConfigError(f"Unsupported input format: {config.input.format}")
