from __future__ import annotations

from dataclasses import asdict, dataclass, field
from numbers import Integral
from pathlib import Path
from typing import Any

import yaml


class ConfigError(ValueError):
    """Raised when a CSCN config file is invalid."""


def _resolve_path(base_dir: Path, value: str | None) -> Path | None:
    if value is None:
        return None
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return path


def _as_list(value: Any) -> list[str] | None:
    if value is None:
        return None
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ConfigError("Expected a list of strings.")
    return value


def _as_string_list(
    mapping: dict[str, Any],
    key: str,
    *,
    default: list[str] | None = None,
) -> list[str]:
    value = mapping.get(key, default)
    if value is None:
        return []
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ConfigError(f"`{key}` must be a list of strings.")
    return [str(item).strip() for item in value]


def _get_string(
    mapping: dict[str, Any],
    key: str,
    *,
    default: str | None = None,
    preserve_blank: bool = False,
) -> str | None:
    if key not in mapping:
        return default
    value = mapping.get(key)
    if value is None:
        return default
    text = str(value)
    if text == "" and not preserve_blank:
        return default
    return text


def _infer_run_name(config_path: Path, raw: dict[str, Any]) -> str:
    run_name = raw.get("run_name")
    if run_name is None:
        return config_path.stem
    if not isinstance(run_name, str) or not run_name.strip():
        raise ConfigError("`run_name` must be a non-empty string.")
    return run_name.strip()


def _get_optional_float(mapping: dict[str, Any], key: str) -> float | None:
    if key not in mapping or mapping.get(key) in (None, ""):
        return None
    return float(mapping[key])


@dataclass(frozen=True)
class InputConfig:
    format: str
    path: Path | None = None
    expr_path: Path | None = None
    metadata_path: Path | None = None
    spatial_x_key: str | None = None
    spatial_y_key: str | None = None
    spatial_z_key: str | None = None
    layer: str | None = None
    gene_key: str | None = None
    obs_group_key: str | None = None
    expr_orientation: str = "cells_by_genes"
    expr_cell_id_column: str | None = None
    metadata_cell_id_column: str = "cell_id"
    expr_delimiter: str | None = None
    metadata_delimiter: str | None = None
    obs_cell_id_key: str | None = None


@dataclass(frozen=True)
class GeneSelectionConfig:
    strategy: str = "variance"
    top_n: int = 150
    gene_list_path: Path | None = None


@dataclass(frozen=True)
class PreprocessConfig:
    normalize: bool = True
    log1p: bool = True
    sample_per_group: int | None = None
    random_seed: int = 42
    gene_selection: GeneSelectionConfig = field(default_factory=GeneSelectionConfig)


@dataclass(frozen=True)
class RunConfig:
    output_dir: Path
    groups: list[str] | None = None
    max_workers: int | None = None
    sigmoid_score: float = 0.1
    significance_level: float = 0.01
    max_cond_vars: int = 20
    use_bitmap: bool = True
    using_nmf: bool = False
    show_progress: bool = False
    progress_interval: int = 100
    spatial: "SpatialConfig" = field(default_factory=lambda: SpatialConfig())


@dataclass(frozen=True)
class SpatialConfig:
    enabled: bool = False
    strategy: str = "weighted_counts"
    mode: str = "knn"
    k: int = 8
    radius: float | None = None
    kernel: str = "gaussian"
    bandwidth: float | None = None
    lambda_expr: float = 0.2
    min_effective_neighbors: int = 15
    block_generator: str = "density_cluster_bbox"
    assignment: str = "nearest_block"
    local_subset: str = "block_plus_halo"
    prior_mode: str = "skeleton_only"
    aggregation_modes: list[str] = field(
        default_factory=lambda: ["sum_then_normalize", "mean"]
    )
    default_aggregation: str = "sum_then_normalize"
    block_gene_top_n: int = 50
    min_cells_per_block: int = 5
    halo_neighbor_blocks: int = 1
    block_overlap_min: float = 0.25
    prior_consensus_threshold: str | int = "auto"
    fallback_to_local_knn_subset: bool = True
    density_clustering: dict[str, Any] = field(
        default_factory=lambda: {"eps": 1.5, "min_samples": 2}
    )


@dataclass(frozen=True)
class AggregateConfig:
    consensus: bool = True
    consensus_threshold_mode: str | int = "auto"


@dataclass(frozen=True)
class BiomarkerConfig:
    enabled: bool = False
    case_group: str | None = None
    control_group: str | None = None
    outcome_column: str = "DISEASE"
    confounder_method: str = "classic"


@dataclass(frozen=True)
class ViewerConfig:
    enabled: bool = False
    host: str = "127.0.0.1"
    port: int = 8000
    reload: bool = False


@dataclass(frozen=True)
class CSCNConfig:
    config_path: Path
    run_name: str
    input: InputConfig
    preprocess: PreprocessConfig
    run: RunConfig
    aggregate: AggregateConfig
    biomarker: BiomarkerConfig
    viewer: ViewerConfig


def load_config(config_path: str | Path) -> CSCNConfig:
    path = Path(config_path).expanduser().resolve()
    if not path.is_file():
        raise ConfigError(f"Config file does not exist: {path}")

    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        raise ConfigError("The config file must contain a top-level mapping.")

    base_dir = path.parent
    run_name = _infer_run_name(path, raw)

    input_raw = raw.get("input") or {}
    if not isinstance(input_raw, dict):
        raise ConfigError("`input` must be a mapping.")
    input_format = str(input_raw.get("format") or "").strip().lower()
    if input_format not in {"h5ad", "tables"}:
        raise ConfigError("`input.format` must be `h5ad` or `tables`.")
    input_config = InputConfig(
        format=input_format,
        path=_resolve_path(base_dir, input_raw.get("path")),
        expr_path=_resolve_path(base_dir, input_raw.get("expr_path")),
        metadata_path=_resolve_path(base_dir, input_raw.get("metadata_path")),
        spatial_x_key=_get_string(input_raw, "spatial_x_key"),
        spatial_y_key=_get_string(input_raw, "spatial_y_key"),
        spatial_z_key=_get_string(input_raw, "spatial_z_key"),
        layer=input_raw.get("layer"),
        gene_key=input_raw.get("gene_key"),
        obs_group_key=input_raw.get("obs_group_key"),
        expr_orientation=str(input_raw.get("expr_orientation") or "cells_by_genes"),
        expr_cell_id_column=_get_string(input_raw, "expr_cell_id_column"),
        metadata_cell_id_column=_get_string(
            input_raw,
            "metadata_cell_id_column",
            default="cell_id",
            preserve_blank=True,
        )
        or "",
        expr_delimiter=input_raw.get("expr_delimiter"),
        metadata_delimiter=input_raw.get("metadata_delimiter"),
        obs_cell_id_key=input_raw.get("obs_cell_id_key"),
    )
    if input_config.format == "h5ad" and input_config.path is None:
        raise ConfigError("`input.path` is required for `h5ad` input.")
    if input_config.format == "tables":
        if input_config.expr_path is None or input_config.metadata_path is None:
            raise ConfigError(
                "`input.expr_path` and `input.metadata_path` are required for `tables` input."
            )
        if input_config.expr_orientation not in {"cells_by_genes", "genes_by_cells"}:
            raise ConfigError(
                "`input.expr_orientation` must be `cells_by_genes` or `genes_by_cells`."
            )
    if bool(input_config.spatial_x_key) != bool(input_config.spatial_y_key):
        raise ConfigError("`input.spatial_x_key` and `input.spatial_y_key` must be set together.")

    preprocess_raw = raw.get("preprocess") or {}
    if not isinstance(preprocess_raw, dict):
        raise ConfigError("`preprocess` must be a mapping.")
    gene_selection_raw = preprocess_raw.get("gene_selection") or {}
    if not isinstance(gene_selection_raw, dict):
        raise ConfigError("`preprocess.gene_selection` must be a mapping.")
    gene_selection = GeneSelectionConfig(
        strategy=str(gene_selection_raw.get("strategy") or "variance"),
        top_n=int(gene_selection_raw.get("top_n") or 150),
        gene_list_path=_resolve_path(base_dir, gene_selection_raw.get("gene_list_path")),
    )
    preprocess = PreprocessConfig(
        normalize=bool(preprocess_raw.get("normalize", True)),
        log1p=bool(preprocess_raw.get("log1p", True)),
        sample_per_group=(
            None
            if preprocess_raw.get("sample_per_group") in (None, "")
            else int(preprocess_raw["sample_per_group"])
        ),
        random_seed=int(preprocess_raw.get("random_seed") or 42),
        gene_selection=gene_selection,
    )
    if preprocess.sample_per_group is not None and preprocess.sample_per_group <= 0:
        raise ConfigError("`preprocess.sample_per_group` must be positive.")
    if preprocess.gene_selection.top_n <= 0:
        raise ConfigError("`preprocess.gene_selection.top_n` must be positive.")

    run_raw = raw.get("run") or {}
    if not isinstance(run_raw, dict):
        raise ConfigError("`run` must be a mapping.")
    spatial_raw = run_raw.get("spatial") or {}
    if not isinstance(spatial_raw, dict):
        raise ConfigError("`run.spatial` must be a mapping.")
    output_dir = _resolve_path(base_dir, run_raw.get("output_dir") or "runs")
    assert output_dir is not None
    run = RunConfig(
        output_dir=output_dir,
        groups=_as_list(run_raw.get("groups")),
        max_workers=(
            None if run_raw.get("max_workers") in (None, "") else int(run_raw["max_workers"])
        ),
        sigmoid_score=float(run_raw.get("sigmoid_score") or 0.1),
        significance_level=float(run_raw.get("significance_level") or 0.01),
        max_cond_vars=int(run_raw.get("max_cond_vars") or 20),
        use_bitmap=bool(run_raw.get("use_bitmap", True)),
        using_nmf=bool(run_raw.get("using_nmf", False)),
        show_progress=bool(run_raw.get("show_progress", False)),
        progress_interval=int(run_raw.get("progress_interval") or 100),
        spatial=SpatialConfig(
            enabled=bool(spatial_raw.get("enabled", False)),
            strategy=str(spatial_raw.get("strategy") or "weighted_counts").strip().lower(),
            mode=str(spatial_raw.get("mode") or "knn").strip().lower(),
            k=int(spatial_raw.get("k", 8)),
            radius=_get_optional_float(spatial_raw, "radius"),
            kernel=str(spatial_raw.get("kernel") or "gaussian").strip().lower(),
            bandwidth=_get_optional_float(spatial_raw, "bandwidth"),
            lambda_expr=float(spatial_raw.get("lambda_expr", 0.2)),
            min_effective_neighbors=int(spatial_raw.get("min_effective_neighbors", 15)),
            block_generator=str(spatial_raw.get("block_generator") or "density_cluster_bbox")
            .strip()
            .lower(),
            assignment=str(spatial_raw.get("assignment") or "nearest_block").strip().lower(),
            local_subset=str(spatial_raw.get("local_subset") or "block_plus_halo")
            .strip()
            .lower(),
            prior_mode=str(spatial_raw.get("prior_mode") or "skeleton_only").strip().lower(),
            aggregation_modes=_as_string_list(
                spatial_raw,
                "aggregation_modes",
                default=["sum_then_normalize", "mean"],
            ),
            default_aggregation=str(
                spatial_raw.get("default_aggregation") or "sum_then_normalize"
            )
            .strip()
            .lower(),
            block_gene_top_n=int(spatial_raw.get("block_gene_top_n") or 50),
            min_cells_per_block=int(spatial_raw.get("min_cells_per_block") or 5),
            halo_neighbor_blocks=int(spatial_raw.get("halo_neighbor_blocks") or 1),
            block_overlap_min=float(spatial_raw.get("block_overlap_min", 0.25)),
            prior_consensus_threshold=spatial_raw.get("prior_consensus_threshold", "auto"),
            fallback_to_local_knn_subset=bool(
                spatial_raw.get("fallback_to_local_knn_subset", True)
            ),
            density_clustering=dict(spatial_raw.get("density_clustering") or {}),
        ),
    )
    if run.spatial.strategy not in {
        "weighted_counts",
        "local_knn_subset",
        "adaptive_block_prior",
    }:
        raise ConfigError(
            "`run.spatial.strategy` must be `weighted_counts`, `local_knn_subset`, or "
            "`adaptive_block_prior`."
        )
    if run.spatial.mode not in {"knn", "radius"}:
        raise ConfigError("`run.spatial.mode` must be `knn` or `radius`.")
    if run.spatial.kernel not in {"binary", "gaussian"}:
        raise ConfigError("`run.spatial.kernel` must be `binary` or `gaussian`.")
    if run.spatial.k <= 0:
        raise ConfigError("`run.spatial.k` must be positive.")
    if not 0 <= run.spatial.lambda_expr <= 1:
        raise ConfigError("`run.spatial.lambda_expr` must be between 0 and 1.")
    if run.spatial.min_effective_neighbors <= 0:
        raise ConfigError("`run.spatial.min_effective_neighbors` must be positive.")
    allowed_aggregation_modes = {"sum_then_normalize", "mean"}
    if not run.spatial.aggregation_modes:
        raise ConfigError("`run.spatial.aggregation_modes` must include at least one mode.")
    if any(mode not in allowed_aggregation_modes for mode in run.spatial.aggregation_modes):
        raise ConfigError(
            "`run.spatial.aggregation_modes` must only contain `sum_then_normalize` or `mean`."
        )
    if run.spatial.default_aggregation not in allowed_aggregation_modes:
        raise ConfigError(
            "`run.spatial.default_aggregation` must be `sum_then_normalize` or `mean`."
        )
    if run.spatial.default_aggregation not in run.spatial.aggregation_modes:
        raise ConfigError(
            "`run.spatial.default_aggregation` must appear in `run.spatial.aggregation_modes`."
        )
    if run.spatial.block_gene_top_n <= 0:
        raise ConfigError("`run.spatial.block_gene_top_n` must be positive.")
    if run.spatial.min_cells_per_block <= 0:
        raise ConfigError("`run.spatial.min_cells_per_block` must be positive.")
    if run.spatial.halo_neighbor_blocks < 0:
        raise ConfigError("`run.spatial.halo_neighbor_blocks` must be non-negative.")
    if not 0 <= run.spatial.block_overlap_min < 1:
        raise ConfigError("`run.spatial.block_overlap_min` must be in [0, 1).")
    if not (
        isinstance(run.spatial.prior_consensus_threshold, str)
        or isinstance(run.spatial.prior_consensus_threshold, Integral)
    ):
        raise ConfigError(
            "`run.spatial.prior_consensus_threshold` must be `auto` or an integer."
        )
    if run.spatial.block_generator != "density_cluster_bbox":
        raise ConfigError(
            "`run.spatial.block_generator` currently only supports `density_cluster_bbox`."
        )
    if run.spatial.assignment != "nearest_block":
        raise ConfigError("`run.spatial.assignment` currently only supports `nearest_block`.")
    if run.spatial.local_subset != "block_plus_halo":
        raise ConfigError("`run.spatial.local_subset` currently only supports `block_plus_halo`.")
    if run.spatial.prior_mode != "skeleton_only":
        raise ConfigError("`run.spatial.prior_mode` currently only supports `skeleton_only`.")
    density_eps = float((run.spatial.density_clustering or {}).get("eps", 1.5))
    density_min_samples = int((run.spatial.density_clustering or {}).get("min_samples", 2))
    if density_eps <= 0:
        raise ConfigError("`run.spatial.density_clustering.eps` must be positive.")
    if density_min_samples <= 0:
        raise ConfigError("`run.spatial.density_clustering.min_samples` must be positive.")
    if run.spatial.enabled and run.spatial.mode == "radius" and (
        run.spatial.radius is None or run.spatial.radius <= 0
    ):
        raise ConfigError("`run.spatial.radius` must be positive when spatial radius mode is enabled.")
    if run.spatial.enabled and run.spatial.strategy == "local_knn_subset":
        if not input_config.spatial_x_key or not input_config.spatial_y_key:
            raise ConfigError(
                "`run.spatial.strategy=local_knn_subset` requires spatial coordinate keys."
            )
        if run.spatial.mode != "knn":
            raise ConfigError(
                "`run.spatial.strategy=local_knn_subset` currently requires `run.spatial.mode=knn`."
            )
    if run.spatial.enabled and run.spatial.strategy == "adaptive_block_prior":
        if not input_config.spatial_x_key or not input_config.spatial_y_key:
            raise ConfigError(
                "`run.spatial.strategy=adaptive_block_prior` requires spatial coordinate keys."
            )
        if input_config.spatial_z_key:
            raise ConfigError(
                "`run.spatial.strategy=adaptive_block_prior` currently supports only 2D coordinates."
            )
        if run.using_nmf:
            raise ConfigError(
                "`run.using_nmf=true` is not supported with `run.spatial.strategy=adaptive_block_prior`."
            )

    aggregate_raw = raw.get("aggregate") or {}
    if not isinstance(aggregate_raw, dict):
        raise ConfigError("`aggregate` must be a mapping.")
    threshold_mode = aggregate_raw.get("consensus_threshold_mode", "auto")
    if not (
        isinstance(threshold_mode, str) or isinstance(threshold_mode, Integral)
    ):
        raise ConfigError(
            "`aggregate.consensus_threshold_mode` must be `auto` or an integer."
        )
    aggregate = AggregateConfig(
        consensus=bool(aggregate_raw.get("consensus", True)),
        consensus_threshold_mode=threshold_mode,
    )

    biomarker_raw = raw.get("biomarker") or {}
    if not isinstance(biomarker_raw, dict):
        raise ConfigError("`biomarker` must be a mapping.")
    biomarker = BiomarkerConfig(
        enabled=bool(biomarker_raw.get("enabled", False)),
        case_group=biomarker_raw.get("case_group"),
        control_group=biomarker_raw.get("control_group"),
        outcome_column=str(biomarker_raw.get("outcome_column") or "DISEASE"),
        confounder_method=str(biomarker_raw.get("confounder_method") or "classic"),
    )

    viewer_raw = raw.get("viewer") or {}
    if not isinstance(viewer_raw, dict):
        raise ConfigError("`viewer` must be a mapping.")
    viewer = ViewerConfig(
        enabled=bool(viewer_raw.get("enabled", False)),
        host=str(viewer_raw.get("host") or "127.0.0.1"),
        port=int(viewer_raw.get("port") or 8000),
        reload=bool(viewer_raw.get("reload", False)),
    )

    return CSCNConfig(
        config_path=path,
        run_name=run_name,
        input=input_config,
        preprocess=preprocess,
        run=run,
        aggregate=aggregate,
        biomarker=biomarker,
        viewer=viewer,
    )


def _serialize_value(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, list):
        return [_serialize_value(item) for item in value]
    if isinstance(value, dict):
        return {key: _serialize_value(item) for key, item in value.items()}
    return value


def serialize_config(config: CSCNConfig) -> dict[str, Any]:
    raw = asdict(config)
    return _serialize_value(raw)
