from __future__ import annotations

import csv
import pickle
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from cscn.config import load_config
from cscn.config import ConfigError
from cscn.layout import RunLayout
from cscn.postprocess import run_ckm
from cscn.viewer import discover_run_roots, load_run_graph, scan_run_root
from cscn.workflow import aggregate_run, load_prepared_run, prepare_run, run_biomarker_workflow, run_cscn
from scripts.analysis.compare_seqfish_ckm_clustering import run_analysis
from tests.support.fake_graph import FakeGraph


def test_prepare_run_writes_standard_layout_from_table_inputs(tmp_path):
    config_path = build_table_config(tmp_path, sample_per_group=None, top_n=2)

    config = load_config(config_path)
    summary = prepare_run(config)
    layout = RunLayout.from_config(config)

    assert summary.run_dir == layout.run_dir
    assert layout.genes_path.is_file()
    assert layout.group_manifest_path.is_file()
    assert layout.cell_metadata_path.is_file()
    assert layout.matrix_path("A").is_file()
    assert layout.matrix_path("B").is_file()

    genes = pd.read_csv(layout.genes_path)["gene_name"].tolist()
    assert genes == ["G2", "G3"]

    groups = pd.read_csv(layout.group_manifest_path)
    assert groups["group_key"].tolist() == ["A", "B"]
    assert groups["n_cells"].tolist() == [2, 2]


def test_aggregate_run_writes_consensus_csv_and_viewer_can_read_run_layout(tmp_path):
    config_path = build_table_config(tmp_path, sample_per_group=None, top_n=2)
    config = load_config(config_path)
    prepare_run(config)
    layout = RunLayout.from_config(config)
    create_fake_dags(layout, "A", [FakeGraph(nodes=[0, 1], edges=[(0, 1)]), FakeGraph(nodes=[0, 1], edges=[(0, 1)])])
    create_fake_dags(layout, "B", [FakeGraph(nodes=[0, 1], edges=[(1, 0)]), FakeGraph(nodes=[0, 1], edges=[(1, 0)])])

    aggregate_run(config)

    consensus_path = layout.consensus_csv_path("A")
    assert consensus_path.is_file()
    rows = list(csv.DictReader(consensus_path.open("r", encoding="utf-8", newline="")))
    assert rows[0]["from"] == "G2"
    assert rows[0]["to"] == "G3"
    assert rows[0]["threshold"] == "2"

    roots = discover_run_roots(layout.base_output_dir, default_root=layout.run_dir)
    assert roots["defaultRootPath"] == str(layout.run_dir)

    scan_payload = scan_run_root(layout.run_dir)
    assert scan_payload["rootKind"] == "run"
    assert scan_payload["collectionKind"] == "single"
    assert scan_payload["groups"][0]["method"] == layout.run_name

    single_payload = load_run_graph(
        root_path=layout.run_dir,
        method=layout.run_name,
        group_key="A",
        mode="single",
        result_id=0,
    )
    assert single_payload["metadata"]["cellRunId"] == "c1"
    assert {node["label"] for node in single_payload["nodes"]} == {"G2", "G3"}

    consensus_payload = load_run_graph(
        root_path=layout.run_dir,
        method=layout.run_name,
        group_key="A",
        mode="consensus",
    )
    assert consensus_payload["metadata"]["threshold"] == 2
    assert consensus_payload["edges"][0]["count"] == 2


def test_biomarker_requires_case_and_control_groups(tmp_path):
    config_path = build_table_config(tmp_path, sample_per_group=None, top_n=2)
    config = load_config(config_path)
    prepare_run(config)

    try:
        run_biomarker_workflow(config)
    except ValueError as exc:
        assert "case_group" in str(exc)
    else:
        raise AssertionError("Expected biomarker workflow to reject missing case/control groups")


def test_prepare_run_accepts_blank_metadata_cell_id_column(tmp_path):
    expr_path = tmp_path / "expression.csv"
    metadata_path = tmp_path / "metadata.csv"
    expr_path.write_text(
        "\n".join(
            [
                "cell_id,G1,G2",
                "c1,1,10",
                "c2,2,20",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    metadata_path.write_text(
        "\n".join(
            [
                ",condition",
                "c1,A",
                "c2,A",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "\n".join(
            [
                "run_name: blank_metadata_id",
                "input:",
                "  format: tables",
                "  expr_path: expression.csv",
                "  metadata_path: metadata.csv",
                "  expr_orientation: cells_by_genes",
                "  expr_cell_id_column: cell_id",
                '  metadata_cell_id_column: ""',
                "  obs_group_key: condition",
                "preprocess:",
                "  gene_selection:",
                "    top_n: 2",
                "run:",
                "  output_dir: runs",
                "aggregate:",
                "  consensus: false",
                "biomarker:",
                "  enabled: false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    config = load_config(config_path)
    summary = prepare_run(config)

    assert config.input.metadata_cell_id_column == ""
    assert summary.groups == {"A": 2}


def test_spatial_config_parses_and_prepare_persists_coords(tmp_path):
    config_path = build_table_config(tmp_path, sample_per_group=None, top_n=2, include_spatial=True)

    config = load_config(config_path)
    assert config.input.spatial_x_key == "array_row"
    assert config.input.spatial_y_key == "array_col"
    assert config.run.spatial.enabled is True
    assert config.run.spatial.strategy == "weighted_counts"
    assert config.run.spatial.mode == "knn"
    assert config.run.spatial.k == 2
    assert config.run.spatial.lambda_expr == 0.0

    prepare_run(config)
    layout = RunLayout.from_config(config)
    assert layout.spatial_coords_path("A").is_file()

    coords = pd.read_csv(layout.spatial_coords_path("A"))
    assert coords.columns.tolist() == ["cell_id", "spatial_x", "spatial_y"]
    assert coords["cell_id"].tolist() == ["c1", "c2"]

    prepared = load_prepared_run(layout)
    assert prepared.groups["A"].spatial_coords is not None
    assert prepared.groups["A"].spatial_coords.shape == (2, 2)


def test_spatial_config_rejects_invalid_values(tmp_path):
    config_path = build_table_config(
        tmp_path,
        sample_per_group=None,
        top_n=2,
        include_spatial=True,
        spatial_overrides=[
            "    mode: invalid",
            "    k: 2",
            "    kernel: gaussian",
            "    lambda_expr: 0.2",
            "    min_effective_neighbors: 2",
        ],
    )
    try:
        load_config(config_path)
    except ConfigError as exc:
        assert "run.spatial.mode" in str(exc)
    else:
        raise AssertionError("Expected invalid spatial mode to be rejected")

    bad_lambda_path = build_table_config(
        tmp_path / "bad_lambda",
        sample_per_group=None,
        top_n=2,
        include_spatial=True,
        spatial_overrides=[
            "    mode: knn",
            "    k: 2",
            "    kernel: gaussian",
            "    lambda_expr: 1.5",
            "    min_effective_neighbors: 2",
        ],
    )
    try:
        load_config(bad_lambda_path)
    except ConfigError as exc:
        assert "lambda_expr" in str(exc)
    else:
        raise AssertionError("Expected invalid lambda_expr to be rejected")

    bad_strategy_path = build_table_config(
        tmp_path / "bad_strategy",
        sample_per_group=None,
        top_n=2,
        include_spatial=True,
        spatial_overrides=[
            "    strategy: unsupported",
            "    mode: knn",
            "    k: 2",
        ],
    )
    try:
        load_config(bad_strategy_path)
    except ConfigError as exc:
        assert "strategy" in str(exc)
    else:
        raise AssertionError("Expected invalid strategy to be rejected")


def test_local_knn_subset_requires_spatial_keys_and_knn_mode(tmp_path):
    no_coord_config = build_table_config(
        tmp_path / "local_no_coord",
        sample_per_group=None,
        top_n=2,
        include_spatial=False,
    )
    no_coord_path = Path(no_coord_config)
    no_coord_path.write_text(
        no_coord_path.read_text(encoding="utf-8").replace(
            "run:\n  output_dir: runs\n",
            "run:\n  output_dir: runs\n  spatial:\n    enabled: true\n    strategy: local_knn_subset\n    mode: knn\n    k: 2\n",
        ),
        encoding="utf-8",
    )
    try:
        load_config(no_coord_path)
    except ConfigError as exc:
        assert "requires spatial coordinate keys" in str(exc)
    else:
        raise AssertionError("Expected local_knn_subset to require spatial keys")

    radius_path = build_table_config(
        tmp_path / "local_radius",
        sample_per_group=None,
        top_n=2,
        include_spatial=True,
        spatial_overrides=[
            "    enabled: true",
            "    strategy: local_knn_subset",
            "    mode: radius",
            "    radius: 1.5",
            "    k: 2",
        ],
    )
    try:
        load_config(radius_path)
    except ConfigError as exc:
        assert "requires `run.spatial.mode=knn`" in str(exc)
    else:
        raise AssertionError("Expected local_knn_subset to reject radius mode")


def test_spatial_input_requires_paired_keys(tmp_path):
    expr_path = tmp_path / "expression.csv"
    metadata_path = tmp_path / "metadata.csv"
    expr_path.write_text(
        "\n".join(
            [
                "cell_id,G1,G2",
                "c1,1,10",
                "c2,2,20",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    metadata_path.write_text(
        "\n".join(
            [
                "cell_id,condition,array_row",
                "c1,A,0",
                "c2,A,1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "\n".join(
            [
                "run_name: missing_spatial_pair",
                "input:",
                "  format: tables",
                "  expr_path: expression.csv",
                "  metadata_path: metadata.csv",
                "  expr_orientation: cells_by_genes",
                "  expr_cell_id_column: cell_id",
                "  metadata_cell_id_column: cell_id",
                "  obs_group_key: condition",
                "  spatial_x_key: array_row",
                "preprocess:",
                "  gene_selection:",
                "    top_n: 2",
                "run:",
                "  output_dir: runs",
                "aggregate:",
                "  consensus: false",
                "biomarker:",
                "  enabled: false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    try:
        load_config(config_path)
    except ConfigError as exc:
        assert "spatial_x_key" in str(exc)
    else:
        raise AssertionError("Expected missing spatial key pair to be rejected")


def test_run_cscn_with_spatial_coords_produces_dags(tmp_path):
    pytest.importorskip("scipy")
    pytest.importorskip("sklearn")
    pytest.importorskip("pgmpy")
    config_path = build_table_config(tmp_path, sample_per_group=None, top_n=2, include_spatial=True)
    config = load_config(config_path)
    prepare_run(config)
    summary = run_cscn(config)
    layout = RunLayout.from_config(config)

    assert summary.groups == {"A": 2, "B": 2}
    assert len(list(layout.dag_group_dir("A").glob("result_*.pkl"))) == 2
    assert len(list(layout.dag_group_dir("B").glob("result_*.pkl"))) == 2


def test_spatial_weighted_counts_support_knn_and_radius():
    pytest.importorskip("scipy")
    pytest.importorskip("sklearn")
    from cscn.core import CSCN

    matrix = pd.DataFrame(
        [
            [0.0, 0.0],
            [1.0, 1.0],
            [2.0, 2.0],
        ]
    ).to_numpy()
    coords = pd.DataFrame(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [5.0, 0.0],
        ]
    ).to_numpy()

    knn_cscn = CSCN(
        sigmoid_score=1.0,
        spatial_enabled=True,
        spatial_mode="knn",
        spatial_k=2,
        spatial_kernel="gaussian",
        spatial_lambda_expr=0.0,
        spatial_min_effective_neighbors=1,
    )
    knn_cscn.run_core(matrix, spatial_coords=coords)
    knn_count = knn_cscn.get_weighted_conditional_counts({0}, key_cell_idx=0)
    assert abs(knn_count - (1.0 + 0.6065306597)) < 1e-6

    radius_cscn = CSCN(
        sigmoid_score=1.0,
        spatial_enabled=True,
        spatial_mode="radius",
        spatial_radius=1.5,
        spatial_kernel="binary",
        spatial_lambda_expr=0.0,
        spatial_min_effective_neighbors=1,
    )
    radius_cscn.run_core(matrix, spatial_coords=coords)
    radius_count = radius_cscn.get_weighted_conditional_counts({0}, key_cell_idx=0)
    assert radius_count == 2.0


def test_local_knn_subset_indices_and_local_df_are_correct():
    pytest.importorskip("scipy")
    pytest.importorskip("sklearn")
    from cscn.core import CSCN

    matrix = pd.DataFrame(
        [
            [10.0, 0.0],
            [20.0, 1.0],
            [30.0, 2.0],
            [40.0, 3.0],
        ]
    ).to_numpy()
    coords = pd.DataFrame(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [5.0, 0.0],
            [6.0, 0.0],
        ]
    ).to_numpy()

    cscn = CSCN(
        spatial_enabled=True,
        spatial_strategy="local_knn_subset",
        spatial_mode="knn",
        spatial_k=2,
    )
    cscn.run_core(matrix, spatial_coords=coords)
    subset_indices = cscn.get_local_subset_indices(0)
    assert subset_indices.tolist() == [0, 1]

    local_df, local_key_idx = cscn.build_local_df_for_key_cell(0)
    assert local_key_idx == 0
    assert local_df.shape == (2, 2)
    assert local_df.iloc[0, 0] == 10.0
    assert local_df.iloc[1, 0] == 20.0


def test_spatial_fallback_triggers_for_too_few_effective_neighbors():
    pytest.importorskip("scipy")
    pytest.importorskip("sklearn")
    from cscn.core import CSCN

    matrix = pd.DataFrame(
        [
            [0.0, 0.0],
            [1.0, 1.0],
            [2.0, 2.0],
        ]
    ).to_numpy()
    coords = pd.DataFrame(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [5.0, 0.0],
        ]
    ).to_numpy()

    cscn = CSCN(
        sigmoid_score=1.0,
        spatial_enabled=True,
        spatial_mode="knn",
        spatial_k=1,
        spatial_kernel="binary",
        spatial_lambda_expr=0.0,
        spatial_min_effective_neighbors=4,
    )
    cscn.run_core(matrix, spatial_coords=coords)
    use_spatial, _, _ = cscn._should_use_spatial_counts(0)
    assert use_spatial is False


def test_run_cscn_with_local_knn_subset_produces_dags(tmp_path):
    pytest.importorskip("scipy")
    pytest.importorskip("sklearn")
    pytest.importorskip("pgmpy")
    config_path = build_table_config(
        tmp_path,
        sample_per_group=None,
        top_n=2,
        include_spatial=True,
        spatial_overrides=[
            "    enabled: true",
            "    strategy: local_knn_subset",
            "    mode: knn",
            "    k: 2",
            "    kernel: gaussian",
            "    lambda_expr: 0.0",
            "    min_effective_neighbors: 1",
        ],
    )
    config = load_config(config_path)
    prepare_run(config)
    summary = run_cscn(config)
    layout = RunLayout.from_config(config)

    assert summary.groups == {"A": 2, "B": 2}
    assert len(list(layout.dag_group_dir("A").glob("result_*.pkl"))) == 2
    assert len(list(layout.dag_group_dir("B").glob("result_*.pkl"))) == 2


def test_run_ckm_writes_ckm_matrix(tmp_path):
    pytest.importorskip("scipy")
    pytest.importorskip("sklearn")
    pytest.importorskip("pgmpy")
    config_path = build_table_config(tmp_path, sample_per_group=None, top_n=2)
    config = load_config(config_path)
    prepare_run(config)
    run_cscn(config)
    layout = RunLayout.from_config(config)

    ckm = run_ckm(layout, "A")

    assert ckm.shape == (2, 2)
    assert (layout.ckm_dir / "A_ckm.npy").is_file()


def test_compute_ckm_auto_beta_transform_handles_centered_inputs(tmp_path):
    pytest.importorskip("scipy")
    pytest.importorskip("sklearn")
    pytest.importorskip("pgmpy")

    import networkx as nx
    import numpy as np

    from cscn.core import CSCN

    cscn = CSCN(output_dir=str(tmp_path))
    cscn.run_core(
        np.asarray(
            [
                [-1.2, 0.2],
                [0.3, 1.1],
            ],
            dtype=float,
        )
    )

    dags = [
        (0, nx.DiGraph([(0, 1)])),
        (1, nx.DiGraph([(1, 0)])),
    ]
    ckm = cscn.compute_ckm(dags=dags, beta_transform="auto")

    assert ckm.shape == (2, 2)
    assert np.isfinite(ckm).all()


def test_compare_seqfish_ckm_clustering_script_outputs_tables(tmp_path):
    pytest.importorskip("scipy")
    pytest.importorskip("sklearn")
    pytest.importorskip("pgmpy")

    weighted_config_path = build_table_config(
        tmp_path / "weighted",
        sample_per_group=None,
        top_n=2,
        include_spatial=True,
    )
    local_config_path = build_table_config(
        tmp_path / "local",
        sample_per_group=None,
        top_n=2,
        include_spatial=True,
        spatial_overrides=[
            "    enabled: true",
            "    strategy: local_knn_subset",
            "    mode: knn",
            "    k: 2",
            "    kernel: gaussian",
            "    lambda_expr: 0.0",
            "    min_effective_neighbors: 1",
        ],
    )

    weighted_config = load_config(weighted_config_path)
    local_config = load_config(local_config_path)
    prepare_run(weighted_config)
    prepare_run(local_config)
    run_cscn(weighted_config)
    run_cscn(local_config)

    weighted_layout = RunLayout.from_config(weighted_config)
    local_layout = RunLayout.from_config(local_config)
    output_dir = tmp_path / "analysis"
    run_analysis(
        expr_path=Path(weighted_config_path).parent / "expression.csv",
        metadata_path=Path(weighted_config_path).parent / "metadata.csv",
        weighted_run_dir=weighted_layout.run_dir,
        local_knn_run_dir=local_layout.run_dir,
        output_dir=output_dir,
        enable_umap=False,
    )

    metrics = pd.read_csv(output_dir / "clustering_metrics.csv")
    assignments = pd.read_csv(output_dir / "cell_assignments.csv")
    assert set(metrics["representation"]) == {"expr", "ckm_weighted", "ckm_local_knn"}
    assert assignments.columns.tolist() == [
        "cell_id",
        "cell_class_name",
        "expr_cluster",
        "ckm_weighted_cluster",
        "ckm_local_knn_cluster",
    ]
    assert (output_dir / "expr_features.csv").is_file()
    assert (output_dir / "ckm_weighted_features.csv").is_file()
    assert (output_dir / "ckm_local_knn_features.csv").is_file()


def build_table_config(
    tmp_path: Path,
    *,
    sample_per_group: int | None,
    top_n: int,
    include_spatial: bool = False,
    spatial_overrides: list[str] | None = None,
) -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    expr_path = tmp_path / "expression.csv"
    metadata_path = tmp_path / "metadata.csv"
    expr_path.write_text(
        "\n".join(
            [
                "cell_id,G1,G2,G3",
                "c1,1,10,1",
                "c2,2,20,1",
                "c3,3,30,10",
                "c4,4,40,10",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    metadata_path.write_text(
        "\n".join(
            (
                [
                    "cell_id,condition,cell_class_name,array_row,array_col",
                    "c1,A,A,0,0",
                    "c2,A,A,1,0",
                    "c3,B,B,0,1",
                    "c4,B,B,1,1",
                ]
                if include_spatial
                else [
                    "cell_id,condition,cell_class_name",
                    "c1,A,A",
                    "c2,A,A",
                    "c3,B,B",
                    "c4,B,B",
                ]
            )
        )
        + "\n",
        encoding="utf-8",
    )
    sample_line = "null" if sample_per_group is None else str(sample_per_group)
    config_lines = [
        "run_name: test_run",
        "input:",
        "  format: tables",
        "  expr_path: expression.csv",
        "  metadata_path: metadata.csv",
        "  expr_orientation: cells_by_genes",
        "  expr_cell_id_column: cell_id",
        "  metadata_cell_id_column: cell_id",
        "  obs_group_key: condition",
    ]
    if include_spatial:
        config_lines.extend(
            [
                "  spatial_x_key: array_row",
                "  spatial_y_key: array_col",
            ]
        )
    config_lines.extend(
        [
            "preprocess:",
            f"  sample_per_group: {sample_line}",
            "  random_seed: 7",
            "  gene_selection:",
            f"    top_n: {top_n}",
            "run:",
            "  output_dir: runs",
        ]
    )
    if include_spatial:
        config_lines.extend(
            [
                "  spatial:",
                *(
                    spatial_overrides
                    if spatial_overrides is not None
                    else [
                        "    enabled: true",
                        "    mode: knn",
                        "    k: 2",
                        "    kernel: gaussian",
                        "    lambda_expr: 0.0",
                        "    min_effective_neighbors: 1",
                    ]
                ),
            ]
        )
    config_lines.extend(
        [
            "aggregate:",
            "  consensus: true",
            "  consensus_threshold_mode: auto",
            "biomarker:",
            "  enabled: false",
        ]
    )
    config_path = tmp_path / "config.yaml"
    config_path.write_text("\n".join(config_lines) + "\n", encoding="utf-8")
    return config_path


def create_fake_dags(layout: RunLayout, group_key: str, graphs: list[FakeGraph]) -> None:
    dag_dir = layout.dag_group_dir(group_key)
    dag_dir.mkdir(parents=True, exist_ok=True)
    for index, graph in enumerate(graphs):
        with (dag_dir / f"result_{index}.pkl").open("wb") as handle:
            pickle.dump(graph, handle)
