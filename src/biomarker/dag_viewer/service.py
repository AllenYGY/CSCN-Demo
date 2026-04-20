from __future__ import annotations

import csv
import importlib.util
import math
import pickle
import re
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from numbers import Integral
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_DATA_ROOT = REPO_ROOT / "data"
KNOWN_METHODS = ("kTotal", "kWithin", "Module_Correlation")
RESULT_PATTERN = re.compile(r"^result_(\d+)\.pkl$")
DAY_GROUP_PATTERN = re.compile(r"^Day(?P<time>\d+)_(?P<cell_type>.+)$")


class DagViewerError(ValueError):
    """Raised when the DAG viewer receives an invalid path or graph request."""


@dataclass(frozen=True)
class ViewerRoot:
    root_path: Path
    root_kind: str
    dag_root: Path
    dataset_root: Path | None


def dependency_status() -> dict[str, Any]:
    names = ("fastapi", "uvicorn", "pgmpy", "networkx", "pandas")
    statuses = {name: importlib.util.find_spec(name) is not None for name in names}
    return {
        "dependencies": statuses,
        "supportedRootKinds": ["dataset", "dag"],
    }


def discover_root_options(base_dir: Path | None = None) -> dict[str, Any]:
    data_root = (base_dir or DEFAULT_DATA_ROOT).resolve()
    roots: list[dict[str, Any]] = []
    default_root_path: str | None = None

    if not data_root.is_dir():
        return {"roots": roots, "defaultRootPath": None, "dataRoot": str(data_root)}

    for dataset_dir in sorted(
        [path for path in data_root.iterdir() if path.is_dir() and not path.name.startswith(".")],
        key=lambda path: path.name,
    ):
        try:
            scan = scan_root(dataset_dir)
        except DagViewerError as exc:
            roots.append(
                {
                    "path": str(dataset_dir),
                    "label": dataset_dir.name,
                    "compatible": False,
                    "reason": str(exc),
                }
            )
            continue

        root_option = {
            "path": str(dataset_dir),
            "label": dataset_dir.name,
            "compatible": True,
            "rootKind": scan["rootKind"],
            "collectionKind": scan["collectionKind"],
            "collectionLabel": scan["collectionLabel"],
            "collectionCount": len(scan["collections"]),
        }
        roots.append(root_option)

        if default_root_path is None or dataset_dir.name == "E-GEOD-93593":
            default_root_path = str(dataset_dir)
            if dataset_dir.name == "E-GEOD-93593":
                # Prefer the main showcase dataset when present.
                continue

    return {
        "roots": roots,
        "defaultRootPath": default_root_path,
        "dataRoot": str(data_root),
    }


def resolve_viewer_root(root_path: str | Path) -> ViewerRoot:
    candidate = Path(root_path).expanduser().resolve()
    if not candidate.exists():
        raise DagViewerError(f"路径不存在: {candidate}")
    if not candidate.is_dir():
        raise DagViewerError(f"路径不是目录: {candidate}")

    dataset_dag_root = candidate / "DAG"
    if dataset_dag_root.is_dir():
        collections = _discover_collection_dirs(dataset_dag_root)
        if collections:
            return ViewerRoot(
                root_path=candidate,
                root_kind="dataset",
                dag_root=dataset_dag_root,
                dataset_root=candidate,
            )

    collections = _discover_collection_dirs(candidate)
    if collections:
        dataset_root = _discover_dataset_root(candidate)
        return ViewerRoot(
            root_path=candidate,
            root_kind="dag",
            dag_root=candidate,
            dataset_root=dataset_root,
        )

    raise DagViewerError(
        "目录中未发现兼容的 DAG 结构。需要传入包含 `DAG/` 的数据集根目录，"
        "或直接传入包含方法目录的 DAG 根目录。"
    )


def scan_root(root_path: str | Path) -> dict[str, Any]:
    viewer_root = resolve_viewer_root(root_path)
    collection_dirs = _discover_collection_dirs(viewer_root.dag_root)
    if not collection_dirs:
        raise DagViewerError(f"DAG 根目录为空: {viewer_root.dag_root}")

    collection_meta = _describe_collections(collection_dirs)
    groups: list[dict[str, Any]] = []
    for collection_dir in collection_dirs:
        collection_key = collection_dir.name
        node_map_available = bool(load_node_map(viewer_root.dataset_root, collection_key))
        for group_dir in sorted(
            [path for path in collection_dir.iterdir() if path.is_dir()],
            key=lambda path: _group_sort_key(path.name),
        ):
            result_ids = _discover_result_ids(group_dir)
            if not result_ids:
                continue

            time_value, cell_type = parse_group_name(group_dir.name)
            groups.append(
                {
                    "method": collection_key,
                    "groupKey": group_dir.name,
                    "time": time_value,
                    "cellType": cell_type,
                    "numDags": len(result_ids),
                    "hasNodeMap": node_map_available,
                    "hasConsensusCsv": _consensus_csv_path(
                        viewer_root.dataset_root, collection_key, group_dir.name
                    ).is_file(),
                    "resultIds": result_ids,
                }
            )

    if not groups:
        raise DagViewerError(f"未在 {viewer_root.dag_root} 中找到任何 `result_*.pkl` 文件。")

    return {
        "rootPath": str(viewer_root.root_path),
        "rootKind": viewer_root.root_kind,
        "rootLabel": (viewer_root.dataset_root or viewer_root.root_path).name,
        "methods": [path.name for path in collection_dirs],
        "collections": collection_meta["collections"],
        "collectionKind": collection_meta["kind"],
        "collectionLabel": collection_meta["label"],
        "showCollectionSelector": len(collection_dirs) > 1,
        "groups": groups,
    }


def load_graph(
    root_path: str | Path,
    method: str,
    group_key: str,
    mode: str,
    result_id: int | None = None,
) -> dict[str, Any]:
    viewer_root = resolve_viewer_root(root_path)
    method_dir = viewer_root.dag_root / method
    if not method_dir.is_dir():
        raise DagViewerError(f"方法目录不存在: {method_dir}")

    group_dir = method_dir / group_key
    if not group_dir.is_dir():
        raise DagViewerError(f"group 目录不存在: {group_dir}")

    node_map = load_node_map(viewer_root.dataset_root, method)
    if mode == "single":
        if result_id is None:
            raise DagViewerError("单细胞模式必须提供 resultId。")
        graph_path = group_dir / f"result_{result_id}.pkl"
        if not graph_path.is_file():
            raise DagViewerError(f"找不到 DAG 文件: {graph_path}")
        graph = _load_pickle_graph(graph_path)
        payload = _serialize_graph(
            graph=_normalize_graph(graph, node_map),
            mode=mode,
            method=method,
            group_key=group_key,
            result_id=result_id,
            cell_run_id=_resolve_cell_run_id(
                viewer_root.dataset_root,
                method,
                group_key,
                result_id,
            ),
        )
        return payload

    if mode != "consensus":
        raise DagViewerError("mode 只能是 `single` 或 `consensus`。")

    consensus_csv = _consensus_csv_path(viewer_root.dataset_root, method, group_key)
    if consensus_csv.is_file():
        consensus_graph, metadata = _load_consensus_from_csv(consensus_csv)
    else:
        consensus_graph, metadata = _build_consensus_from_pickles(group_dir, node_map)

    return _serialize_graph(
        graph=consensus_graph,
        mode=mode,
        method=method,
        group_key=group_key,
        threshold=metadata.get("threshold"),
        num_dags=metadata.get("num_dags"),
    )


def parse_group_name(group_key: str) -> tuple[str, str]:
    match = DAY_GROUP_PATTERN.match(group_key)
    if match:
        time_value = match.group("time")
        cell_type = match.group("cell_type").replace("_", " ")
        return time_value, cell_type
    if group_key.isdigit():
        return group_key, ""
    return "", group_key.replace("_", " ")


def safe_group_name(time_value: str, cell_type: str) -> str:
    return (
        f"Day{time_value}_{cell_type}".replace(" ", "_").replace("/", "_").replace("\\", "_")
    )


def load_node_map(dataset_root: Path | None, method: str) -> dict[int, str]:
    path = _node_map_path(dataset_root, method)
    if path.is_file():
        mapping: dict[int, str] = {}
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                try:
                    node_index = int((row.get("node_index") or "").strip())
                except ValueError:
                    continue
                mapping[node_index] = (row.get("gene_name") or "").strip()
        if mapping:
            return mapping

    return _load_node_map_from_used_genes(dataset_root, method)


def _discover_collection_dirs(dag_root: Path) -> list[Path]:
    if not dag_root.is_dir():
        return []

    candidates = [path for path in dag_root.iterdir() if path.is_dir()]
    collection_dirs = [path for path in candidates if _looks_like_collection_dir(path)]
    return sorted(collection_dirs, key=lambda path: (path.name not in KNOWN_METHODS, path.name))


def _looks_like_collection_dir(path: Path) -> bool:
    for group_dir in path.iterdir():
        if not group_dir.is_dir():
            continue
        if _discover_result_ids(group_dir):
            return True
    return False


def _describe_collections(collection_dirs: list[Path]) -> dict[str, Any]:
    keys = [path.name for path in collection_dirs]
    if keys and all(key in KNOWN_METHODS for key in keys):
        kind = "method"
        label = "方法"
    elif len(keys) == 1:
        kind = "single"
        label = "当前视图"
    elif any(key.startswith(("GSE", "E-GEOD")) for key in keys):
        kind = "run"
        label = "分析运行"
    else:
        kind = "collection"
        label = "DAG 视图"

    collections = []
    for key in keys:
        if kind == "single":
            display_label = "默认视图"
        elif kind == "method":
            display_label = key
        else:
            display_label = key.replace("_", " ")
        collections.append(
            {
                "key": key,
                "displayLabel": display_label,
                "rawLabel": key,
                "kind": kind,
            }
        )

    return {"kind": kind, "label": label, "collections": collections}


def _discover_result_ids(group_dir: Path) -> list[int]:
    result_ids: list[int] = []
    for path in group_dir.glob("result_*.pkl"):
        match = RESULT_PATTERN.match(path.name)
        if match:
            result_ids.append(int(match.group(1)))
    return sorted(result_ids)


def _discover_dataset_root(dag_root: Path) -> Path | None:
    candidate = dag_root.parent
    if (candidate / "gene").is_dir() or (candidate / "visualizations").is_dir():
        return candidate
    if (candidate / "expression_matrix.csv").is_file():
        return candidate
    if next(candidate.glob("*.sdrf.txt"), None):
        return candidate
    return None


def _node_map_path(dataset_root: Path | None, method: str) -> Path:
    if dataset_root is None:
        return Path("__missing__")
    return dataset_root / "gene" / f"{method}_node_map.csv"


def _consensus_csv_path(dataset_root: Path | None, method: str, group_key: str) -> Path:
    if dataset_root is None:
        return Path("__missing__")
    return dataset_root / "visualizations" / method / "consensus_edges" / f"{group_key}.csv"


def _group_sort_key(group_key: str) -> tuple[int, str, str]:
    time_value, cell_type = parse_group_name(group_key)
    time_sort = int(time_value) if time_value.isdigit() else 999999
    return time_sort, cell_type, group_key


def _load_pickle_graph(graph_path: Path) -> Any:
    try:
        with graph_path.open("rb") as handle:
            return pickle.load(handle)
    except ModuleNotFoundError as exc:
        raise DagViewerError(
            f"无法读取 {graph_path.name}：缺少依赖模块 `{exc.name}`。"
            " 对原始 CSCN pickle 通常需要先安装 `pgmpy`。"
        ) from exc
    except pickle.UnpicklingError as exc:
        raise DagViewerError(f"无法反序列化 {graph_path.name}：文件已损坏。") from exc
    except EOFError as exc:
        raise DagViewerError(f"无法读取 {graph_path.name}：pickle 文件不完整。") from exc
    except AttributeError as exc:
        raise DagViewerError(f"无法读取 {graph_path.name}：pickle 依赖的类定义不可用。") from exc


def _load_consensus_from_csv(csv_path: Path) -> tuple[dict[str, Any], dict[str, int | None]]:
    nodes: dict[str, dict[str, Any]] = {}
    edges: list[dict[str, Any]] = []
    threshold: int | None = None
    num_dags: int | None = None

    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            source = (row.get("from") or "").strip()
            target = (row.get("to") or "").strip()
            if not source or not target or source == target:
                continue
            count = _safe_int(row.get("count"), default=1)
            threshold = _safe_int(row.get("threshold"), default=threshold)
            num_dags = _safe_int(row.get("num_dags"), default=num_dags)

            for node_id in (source, target):
                nodes.setdefault(
                    node_id,
                    {
                        "id": node_id,
                        "label": node_id,
                        "rawId": node_id,
                        "mapped": True,
                    },
                )

            edges.append(
                {
                    "id": f"{source}->{target}",
                    "source": source,
                    "target": target,
                    "count": count,
                    "weight": count,
                }
            )

    graph = {"nodes": list(nodes.values()), "edges": edges}
    return graph, {"threshold": threshold, "num_dags": num_dags}


def _build_consensus_from_pickles(
    group_dir: Path,
    node_map: dict[int, str],
) -> tuple[dict[str, Any], dict[str, int]]:
    result_paths = sorted(group_dir.glob("result_*.pkl"))
    if not result_paths:
        raise DagViewerError(f"group 目录中没有找到 DAG 文件: {group_dir}")

    edge_counts: defaultdict[tuple[str, str], int] = defaultdict(int)
    node_info: dict[str, dict[str, Any]] = {}
    num_dags = 0

    for graph_path in result_paths:
        graph = _load_pickle_graph(graph_path)
        normalized = _normalize_graph(graph, node_map)
        num_dags += 1
        for node in normalized["nodes"]:
            node_info.setdefault(node["id"], node)
        for edge in normalized["edges"]:
            if edge["source"] == edge["target"]:
                continue
            edge_counts[(edge["source"], edge["target"])] += 1

    threshold = max(2, int(math.ceil(0.05 * num_dags)))
    edges: list[dict[str, Any]] = []
    for (source, target), count in sorted(edge_counts.items()):
        if count < threshold:
            continue
        edges.append(
            {
                "id": f"{source}->{target}",
                "source": source,
                "target": target,
                "count": count,
                "weight": count,
            }
        )

    referenced_nodes: dict[str, dict[str, Any]] = {}
    for edge in edges:
        referenced_nodes[edge["source"]] = node_info[edge["source"]]
        referenced_nodes[edge["target"]] = node_info[edge["target"]]

    graph = {"nodes": list(referenced_nodes.values()), "edges": edges}
    return graph, {"threshold": threshold, "num_dags": num_dags}


def _normalize_graph(graph: Any, node_map: dict[int, str]) -> dict[str, Any]:
    raw_nodes = [_serialize_node_ref(node) for node in _graph_nodes(graph)]
    node_by_raw = {node["rawRef"]: node for node in raw_nodes}
    edges: list[dict[str, Any]] = []

    for source, target in _graph_edges(graph):
        source_node = node_by_raw[_raw_node_key(source)]
        target_node = node_by_raw[_raw_node_key(target)]
        edges.append(
            {
                "id": f"{source_node['id']}->{target_node['id']}",
                "source": source_node["id"],
                "target": target_node["id"],
                "count": 1,
                "weight": 1,
            }
        )

    for node in raw_nodes:
        label, mapped = _display_label(node["rawValue"], node_map)
        node["label"] = label
        node["mapped"] = mapped

    return {"nodes": raw_nodes, "edges": edges}


def _serialize_graph(
    graph: dict[str, Any],
    mode: str,
    method: str,
    group_key: str,
    result_id: int | None = None,
    cell_run_id: str | None = None,
    threshold: int | None = None,
    num_dags: int | None = None,
) -> dict[str, Any]:
    node_lookup = {node["id"]: node for node in graph["nodes"]}
    in_degree: defaultdict[str, int] = defaultdict(int)
    out_degree: defaultdict[str, int] = defaultdict(int)
    for edge in graph["edges"]:
        out_degree[edge["source"]] += 1
        in_degree[edge["target"]] += 1

    nodes = []
    for node in sorted(graph["nodes"], key=lambda item: item["label"].lower()):
        node_id = node["id"]
        indeg = in_degree[node_id]
        outdeg = out_degree[node_id]
        nodes.append(
            {
                "id": node_id,
                "label": node["label"],
                "rawId": node["rawId"],
                "mapped": node["mapped"],
                "degree": indeg + outdeg,
                "inDegree": indeg,
                "outDegree": outdeg,
            }
        )

    edges = sorted(
        graph["edges"],
        key=lambda item: (item["source"], item["target"], item["weight"]),
    )

    metadata: dict[str, Any] = {
        "method": method,
        "groupKey": group_key,
        "numNodes": len(nodes),
        "numEdges": len(edges),
    }
    if result_id is not None:
        metadata["resultId"] = result_id
    if cell_run_id:
        metadata["cellRunId"] = cell_run_id
    if threshold is not None:
        metadata["threshold"] = threshold
    if num_dags is not None:
        metadata["numDags"] = num_dags

    return {
        "mode": mode,
        "nodes": nodes,
        "edges": edges,
        "metadata": metadata,
    }


def _serialize_node_ref(raw_value: Any) -> dict[str, Any]:
    raw_id = str(raw_value)
    return {
        "id": raw_id,
        "label": raw_id,
        "rawId": raw_id,
        "mapped": False,
        "rawRef": _raw_node_key(raw_value),
        "rawValue": raw_value,
    }


def _display_label(raw_value: Any, node_map: dict[int, str]) -> tuple[str, bool]:
    if isinstance(raw_value, Integral) and not isinstance(raw_value, bool):
        mapped = node_map.get(int(raw_value))
        if _is_valid_label(mapped):
            return str(mapped), True
    raw_text = str(raw_value)
    if _is_valid_label(raw_text):
        return raw_text, False
    return raw_text, False


def _is_valid_label(value: Any) -> bool:
    if value is None:
        return False
    text = str(value).strip()
    return text != "" and text.lower() != "nan"


def _raw_node_key(raw_value: Any) -> str:
    return f"{type(raw_value).__name__}:{raw_value!r}"


def _graph_nodes(graph: Any) -> list[Any]:
    nodes_attr = getattr(graph, "nodes", None)
    if callable(nodes_attr):
        return list(nodes_attr())
    if nodes_attr is not None:
        return list(nodes_attr)
    raise DagViewerError("图对象缺少 `nodes()` 接口，无法序列化。")


def _graph_edges(graph: Any) -> list[tuple[Any, Any]]:
    edges_attr = getattr(graph, "edges", None)
    if callable(edges_attr):
        return list(edges_attr())
    if edges_attr is not None:
        return list(edges_attr)
    raise DagViewerError("图对象缺少 `edges()` 接口，无法序列化。")


def _safe_int(raw_value: Any, default: int | None) -> int | None:
    if raw_value in (None, ""):
        return default
    try:
        return int(raw_value)
    except (TypeError, ValueError):
        return default


def _resolve_cell_run_id(
    dataset_root: Path | None,
    method: str,
    group_key: str,
    result_id: int,
) -> str | None:
    if dataset_root is None:
        return None

    group_runs = _load_group_run_ids(str(dataset_root))
    runs = group_runs.get(group_key)
    if runs is not None:
        if result_id < 0 or result_id >= len(runs):
            return None
        resolved = runs[result_id]
        if resolved:
            return resolved

    sampled_runs = _load_sampled_group_run_ids(str(dataset_root), method, group_key)
    if result_id < 0 or result_id >= len(sampled_runs):
        return None
    return sampled_runs[result_id]


@lru_cache(maxsize=8)
def _load_group_run_ids(dataset_root: str) -> dict[str, tuple[str, ...]]:
    dataset_dir = Path(dataset_root)
    expression_path = dataset_dir / "expression_matrix.csv"
    sdrf_path = _discover_sdrf_path(dataset_dir)
    if not expression_path.is_file() or sdrf_path is None:
        return {}

    ordered_runs = _read_expression_runs(expression_path)
    run_meta = _read_sdrf_metadata(sdrf_path)
    groups: defaultdict[str, list[str]] = defaultdict(list)

    for run_id in ordered_runs:
        time_and_cell = run_meta.get(run_id)
        if time_and_cell is None:
            continue
        time_value, cell_type = time_and_cell
        groups[time_value].append(run_id)
        groups[safe_group_name(time_value, cell_type)].append(run_id)

    return {key: tuple(value) for key, value in groups.items()}


@lru_cache(maxsize=64)
def _load_sampled_group_run_ids(dataset_root: str, method: str, group_key: str) -> tuple[str, ...]:
    dataset_dir = Path(dataset_root)
    sampled_path = dataset_dir / "output_deseq" / f"{method}_{group_key}_sampled_cells.csv"
    if not sampled_path.is_file():
        return ()

    values: list[str] = []
    with sampled_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return ()
        cell_column = _first_present(reader.fieldnames, ["cell_id", "cell", "run_id"])
        if cell_column is None:
            return ()
        for row in reader:
            value = (row.get(cell_column) or "").strip()
            if value:
                values.append(value)
    return tuple(values)


def _load_node_map_from_used_genes(dataset_root: Path | None, method: str) -> dict[int, str]:
    if dataset_root is None:
        return {}

    output_dir = dataset_root / "output_deseq"
    if not output_dir.is_dir():
        return {}

    candidates = sorted(output_dir.glob(f"{method}_top*_genes_used.csv"))
    if not candidates:
        return {}

    mapping: dict[int, str] = {}
    with candidates[0].open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return {}
        gene_column = _first_present(reader.fieldnames, ["gene", "Gene", "Unnamed: 0"])
        if gene_column is None:
            return {}
        for index, row in enumerate(reader):
            mapping[index] = (row.get(gene_column) or "").strip()
    return mapping


def _discover_sdrf_path(dataset_dir: Path) -> Path | None:
    matches = sorted(dataset_dir.glob("*.sdrf.txt"))
    return matches[0] if matches else None


def _read_expression_runs(expression_path: Path) -> tuple[str, ...]:
    runs: list[str] = []
    with expression_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        next(reader, None)
        for row in reader:
            if row and row[0].strip():
                runs.append(row[0].strip())
    return tuple(runs)


def _read_sdrf_metadata(sdrf_path: Path) -> dict[str, tuple[str, str]]:
    with sdrf_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return {}
        run_col = _first_present(reader.fieldnames, ["Comment [ENA_RUN]", "Comment[ENA_RUN]"])
        label_col = _first_present(
            reader.fieldnames,
            [
                "Factor Value[inferred cell type - ontology labels]",
                "Factor Value[inferred cell type - authors labels]",
                "Characteristics [cell type]",
            ],
        )
        if run_col is None or label_col is None:
            return {}

        mapping: dict[str, tuple[str, str]] = {}
        for row in reader:
            run_id = (row.get(run_col) or "").strip()
            time_value = (row.get("Factor Value[time]") or "").strip()
            cell_type = (row.get(label_col) or "").strip()
            if not run_id or not time_value or not cell_type:
                continue
            mapping[run_id] = (time_value, cell_type)
        return mapping


def _first_present(values: list[str], candidates: list[str]) -> str | None:
    for candidate in candidates:
        if candidate in values:
            return candidate
    return None
