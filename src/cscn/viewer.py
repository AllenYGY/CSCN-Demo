from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

try:
    from fastapi import FastAPI, HTTPException
    from fastapi.middleware.cors import CORSMiddleware
    from fastapi.responses import FileResponse, HTMLResponse
    from fastapi.staticfiles import StaticFiles
except ModuleNotFoundError as exc:  # pragma: no cover
    FastAPI = None
    HTTPException = None
    CORSMiddleware = None
    FileResponse = None
    HTMLResponse = None
    StaticFiles = None
    FASTAPI_IMPORT_ERROR = exc
else:
    FASTAPI_IMPORT_ERROR = None

from biomarker.dag_viewer.service import (
    DagViewerError,
    _load_consensus_from_csv,
    _load_pickle_graph,
    _normalize_graph,
    _serialize_graph,
    dependency_status,
)

from .aggregate import build_consensus_edges, load_node_map
from .layout import safe_component

REPO_ROOT = Path(__file__).resolve().parents[2]
FRONTEND_DIST = REPO_ROOT / "apps" / "dag-viewer" / "dist"


def _is_run_root(path: Path) -> bool:
    return path.is_dir() and (path / "dags").is_dir() and (path / "inputs" / "genes.csv").is_file()


def discover_run_roots(base_dir: Path, default_root: Path | None = None) -> dict[str, object]:
    base_dir = base_dir.resolve()
    roots = []
    candidates = []
    if _is_run_root(base_dir):
        candidates = [base_dir]
    elif base_dir.is_dir():
        candidates = [
            path
            for path in sorted(base_dir.iterdir(), key=lambda item: item.name)
            if _is_run_root(path)
        ]

    default_root_path = None
    for candidate in candidates:
        roots.append(
            {
                "path": str(candidate),
                "label": candidate.name,
                "compatible": True,
                "rootKind": "run",
                "collectionKind": "single",
                "collectionLabel": "当前视图",
                "collectionCount": 1,
            }
        )
    if default_root and _is_run_root(default_root):
        default_root_path = str(default_root.resolve())
    elif roots:
        default_root_path = str(candidates[0])
    return {
        "roots": roots,
        "defaultRootPath": default_root_path,
        "dataRoot": str(base_dir),
    }


def _load_group_labels(run_root: Path) -> dict[str, str]:
    groups_path = run_root / "inputs" / "groups.csv"
    if not groups_path.is_file():
        return {}
    labels = {}
    with groups_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            key = str(row.get("group_key") or "")
            label = str(row.get("group_label") or key)
            if key:
                labels[key] = label
    return labels


def _discover_result_ids(group_dir: Path) -> list[int]:
    result_ids = []
    for path in group_dir.glob("result_*.pkl"):
        try:
            result_ids.append(int(path.stem.split("_", 1)[1]))
        except (IndexError, ValueError):
            continue
    return sorted(result_ids)


def _block_debug_method(run_root: Path) -> str:
    return f"{safe_component(run_root.name)}__blocks"


def _load_block_node_map(run_root: Path, group_key: str) -> dict[int, str]:
    path = run_root / "blocks" / f"{safe_component(group_key)}_block_genes.csv"
    if not path.is_file():
        return {}
    frame = pd.read_csv(path)
    mapping: dict[int, str] = {}
    for row in frame.to_dict(orient="records"):
        try:
            node_index = int(row["node_index"])
        except (KeyError, TypeError, ValueError):
            continue
        mapping[node_index] = str(row.get("gene_name") or node_index)
    return mapping


def _resolve_block_run_id(run_root: Path, group_key: str, result_id: int | None) -> str | None:
    if result_id is None:
        return None
    manifest_path = run_root / "blocks" / f"{safe_component(group_key)}_block_manifest.csv"
    if not manifest_path.is_file():
        return None
    frame = pd.read_csv(manifest_path)
    match = frame.loc[frame["block_index"] == result_id]
    if match.empty:
        return None
    return str(match.iloc[0]["block_id"])


def scan_run_root(root_path: str | Path) -> dict[str, object]:
    run_root = Path(root_path).expanduser().resolve()
    if not _is_run_root(run_root):
        raise DagViewerError(f"Run 目录不兼容或不存在: {run_root}")

    labels = _load_group_labels(run_root)
    groups = []
    for group_dir in sorted((run_root / "dags").iterdir(), key=lambda item: item.name):
        if not group_dir.is_dir():
            continue
        result_ids = _discover_result_ids(group_dir)
        if not result_ids:
            continue
        group_key = group_dir.name
        groups.append(
            {
                "method": safe_component(run_root.name),
                "groupKey": group_key,
                "time": "",
                "cellType": labels.get(group_key, group_key),
                "numDags": len(result_ids),
                "hasNodeMap": (run_root / "inputs" / "genes.csv").is_file(),
                "hasConsensusCsv": (run_root / "consensus" / f"{group_key}.csv").is_file(),
                "resultIds": result_ids,
            }
        )

    block_debug_dir = run_root / "blocks" / "block_dags"
    if block_debug_dir.is_dir():
        debug_method = _block_debug_method(run_root)
        for group_dir in sorted(block_debug_dir.iterdir(), key=lambda item: item.name):
            if not group_dir.is_dir():
                continue
            result_ids = _discover_result_ids(group_dir)
            if not result_ids:
                continue
            group_key = group_dir.name
            groups.append(
                {
                    "method": debug_method,
                    "groupKey": group_key,
                    "time": "",
                    "cellType": f"{labels.get(group_key, group_key)} (block debug)",
                    "numDags": len(result_ids),
                    "hasNodeMap": (run_root / "blocks" / f"{group_key}_block_genes.csv").is_file(),
                    "hasConsensusCsv": False,
                    "resultIds": result_ids,
                }
            )

    if not groups:
        raise DagViewerError(f"Run 目录中未找到任何 `result_*.pkl` 文件: {run_root}")

    collection_key = safe_component(run_root.name)
    return {
        "rootPath": str(run_root),
        "rootKind": "run",
        "rootLabel": run_root.name,
        "methods": [collection_key, *([_block_debug_method(run_root)] if block_debug_dir.is_dir() else [])],
        "collections": [
            {
                "key": collection_key,
                "displayLabel": "默认视图",
                "rawLabel": run_root.name,
                "kind": "single",
            },
            *(
                [
                    {
                        "key": _block_debug_method(run_root),
                        "displayLabel": "Block Debug",
                        "rawLabel": f"{run_root.name}__blocks",
                        "kind": "single",
                    }
                ]
                if block_debug_dir.is_dir()
                else []
            ),
        ],
        "collectionKind": "single",
        "collectionLabel": "当前视图",
        "showCollectionSelector": block_debug_dir.is_dir(),
        "groups": groups,
    }


def _resolve_cell_run_id(run_root: Path, group_key: str, result_id: int | None) -> str | None:
    if result_id is None:
        return None
    cells_path = run_root / "matrices" / f"{group_key}_cells.csv"
    if not cells_path.is_file():
        return None
    frame = []
    with cells_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        frame = list(reader)
    if result_id < 0 or result_id >= len(frame):
        return None
    return str(frame[result_id].get("cell_id") or "")


def load_run_graph(
    root_path: str | Path,
    method: str,
    group_key: str,
    mode: str,
    result_id: int | None = None,
) -> dict[str, object]:
    run_root = Path(root_path).expanduser().resolve()
    if not _is_run_root(run_root):
        raise DagViewerError(f"Run 目录不兼容或不存在: {run_root}")

    expected_method = safe_component(run_root.name)
    debug_method = _block_debug_method(run_root)
    is_block_debug = method == debug_method
    if method and method not in {expected_method, debug_method}:
        raise DagViewerError(f"方法不存在: {method}")

    group_dir = (
        run_root / "blocks" / "block_dags" / safe_component(group_key)
        if is_block_debug
        else run_root / "dags" / safe_component(group_key)
    )
    if not group_dir.is_dir():
        raise DagViewerError(f"group 目录不存在: {group_dir}")

    node_map = _load_block_node_map(run_root, group_key) if is_block_debug else load_node_map(run_root)
    if mode == "single":
        if result_id is None:
            raise DagViewerError("单细胞模式必须提供 resultId。")
        graph_path = group_dir / f"result_{result_id}.pkl"
        if not graph_path.is_file():
            raise DagViewerError(f"找不到 DAG 文件: {graph_path}")
        graph = _load_pickle_graph(graph_path)
        return _serialize_graph(
            graph=_normalize_graph(graph, node_map),
            mode="single",
            method=debug_method if is_block_debug else expected_method,
            group_key=group_key,
            result_id=result_id,
            cell_run_id=(
                _resolve_block_run_id(run_root, group_key, result_id)
                if is_block_debug
                else _resolve_cell_run_id(run_root, group_key, result_id)
            ),
        )

    if mode != "consensus":
        raise DagViewerError("mode 只能是 `single` 或 `consensus`。")

    consensus_csv = run_root / "consensus" / f"{safe_component(group_key)}.csv"
    if consensus_csv.is_file():
        consensus_graph, metadata = _load_consensus_from_csv(consensus_csv)
    else:
        edges, metadata = build_consensus_edges(
            group_dir=group_dir,
            node_map=node_map,
            threshold_mode="auto",
        )
        nodes = {}
        serialized_edges = []
        for edge in edges:
            source = edge["from"]
            target = edge["to"]
            nodes.setdefault(
                source,
                {"id": source, "label": source, "rawId": source, "mapped": True},
            )
            nodes.setdefault(
                target,
                {"id": target, "label": target, "rawId": target, "mapped": True},
            )
            serialized_edges.append(
                {
                    "id": f"{source}->{target}",
                    "source": source,
                    "target": target,
                    "count": edge["count"],
                    "weight": edge["count"],
                }
            )
        consensus_graph = {"nodes": list(nodes.values()), "edges": serialized_edges}

    return _serialize_graph(
        graph=consensus_graph,
        mode="consensus",
        method=debug_method if is_block_debug else expected_method,
        group_key=group_key,
        threshold=metadata.get("threshold"),
        num_dags=metadata.get("num_dags"),
    )


def create_app(
    base_dir: Path,
    default_root: Path | None = None,
    frontend_dist: Path | None = None,
):
    if FastAPI is None:
        raise RuntimeError(
            "FastAPI 未安装，无法启动 CSCN Viewer。请先安装 `fastapi` 和 `uvicorn`。"
        ) from FASTAPI_IMPORT_ERROR

    dist_dir = frontend_dist or FRONTEND_DIST
    app = FastAPI(title="CSCN Viewer", version="0.1.0")
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    @app.get("/api/health")
    def health():
        status = dependency_status()
        status["frontendBuilt"] = dist_dir.is_dir() and (dist_dir / "index.html").is_file()
        status["frontendDist"] = str(dist_dir)
        return status

    @app.get("/api/roots")
    def roots():
        return discover_run_roots(base_dir=base_dir, default_root=default_root)

    @app.post("/api/scan")
    def scan(payload: dict[str, str]):
        root_path = (payload.get("rootPath") or "").strip()
        if not root_path:
            raise HTTPException(status_code=400, detail="缺少 rootPath。")
        try:
            return scan_run_root(root_path)
        except DagViewerError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc

    @app.post("/api/graph")
    def graph(payload: dict[str, object]):
        root_path = str(payload.get("rootPath") or "").strip()
        method = str(payload.get("method") or "").strip()
        group_key = str(payload.get("groupKey") or "").strip()
        mode = str(payload.get("mode") or "").strip()
        result_id = payload.get("resultId")
        if not root_path or not method or not group_key or not mode:
            raise HTTPException(status_code=400, detail="缺少必要参数。")
        parsed_result_id = None
        if result_id is not None:
            try:
                parsed_result_id = int(result_id)
            except (TypeError, ValueError) as exc:
                raise HTTPException(status_code=400, detail="resultId 必须是整数。") from exc
        try:
            return load_run_graph(
                root_path=root_path,
                method=method,
                group_key=group_key,
                mode=mode,
                result_id=parsed_result_id,
            )
        except DagViewerError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc

    if dist_dir.is_dir() and (dist_dir / "assets").is_dir():
        app.mount("/assets", StaticFiles(directory=dist_dir / "assets"), name="assets")

    @app.get("/", include_in_schema=False)
    def index():
        if dist_dir.is_dir() and (dist_dir / "index.html").is_file():
            return FileResponse(dist_dir / "index.html")
        return HTMLResponse(_missing_frontend_html(dist_dir))

    @app.get("/{full_path:path}", include_in_schema=False)
    def spa_fallback(full_path: str):
        if full_path.startswith("api/"):
            raise HTTPException(status_code=404, detail="未找到 API 路径。")
        candidate = dist_dir / full_path
        if candidate.is_file():
            return FileResponse(candidate)
        if dist_dir.is_dir() and (dist_dir / "index.html").is_file():
            return FileResponse(dist_dir / "index.html")
        return HTMLResponse(_missing_frontend_html(dist_dir))

    return app


def _missing_frontend_html(dist_dir: Path) -> str:
    return f"""
<!doctype html>
<html lang="zh-CN">
  <head>
    <meta charset="utf-8" />
    <title>CSCN Viewer</title>
    <style>
      body {{
        margin: 0;
        font-family: "Avenir Next", "Segoe UI", sans-serif;
        background: linear-gradient(180deg, #f7f4ea 0%, #efe6d2 100%);
        color: #18313a;
      }}
      main {{
        max-width: 840px;
        margin: 64px auto;
        padding: 32px;
        background: rgba(255, 255, 255, 0.78);
        border: 1px solid rgba(24, 49, 58, 0.12);
        border-radius: 24px;
        box-shadow: 0 18px 40px rgba(24, 49, 58, 0.08);
      }}
      code {{
        padding: 2px 6px;
        border-radius: 6px;
        background: #e8efe9;
      }}
    </style>
  </head>
  <body>
    <main>
      <h1>CSCN Viewer</h1>
      <p>后端 API 已就绪，但前端构建产物不存在。</p>
      <p>请在仓库根目录运行：</p>
      <p><code>cd apps/dag-viewer &amp;&amp; npm install &amp;&amp; npm run build</code></p>
      <p>然后重新刷新本页。当前预期的静态目录：</p>
      <p><code>{dist_dir}</code></p>
    </main>
  </body>
</html>
""".strip()
