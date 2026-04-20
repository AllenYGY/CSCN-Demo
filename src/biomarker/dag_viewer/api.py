from __future__ import annotations

from pathlib import Path

try:
    from fastapi import FastAPI, HTTPException
    from fastapi.middleware.cors import CORSMiddleware
    from fastapi.responses import FileResponse, HTMLResponse
    from fastapi.staticfiles import StaticFiles
except ModuleNotFoundError as exc:  # pragma: no cover - exercised only without FastAPI.
    FastAPI = None
    HTTPException = None
    CORSMiddleware = None
    FileResponse = None
    HTMLResponse = None
    StaticFiles = None
    FASTAPI_IMPORT_ERROR = exc
else:
    FASTAPI_IMPORT_ERROR = None

from .service import (
    DagViewerError,
    dependency_status,
    discover_root_options,
    load_graph,
    scan_root,
)

REPO_ROOT = Path(__file__).resolve().parents[3]
FRONTEND_DIST = REPO_ROOT / "apps" / "dag-viewer" / "dist"


def create_app(frontend_dist: Path | None = None):
    if FastAPI is None:
        raise RuntimeError(
            "FastAPI 未安装，无法启动 DAG Viewer API。"
            " 请先安装 `fastapi` 和 `uvicorn`。"
        ) from FASTAPI_IMPORT_ERROR

    dist_dir = frontend_dist or FRONTEND_DIST
    app = FastAPI(title="CSCN DAG Viewer", version="0.1.0")
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
        return discover_root_options()

    @app.post("/api/scan")
    def scan(payload: dict[str, str]):
        root_path = (payload.get("rootPath") or "").strip()
        if not root_path:
            raise HTTPException(status_code=400, detail="缺少 rootPath。")
        try:
            return scan_root(root_path)
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
            return load_graph(
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
    <title>CSCN DAG Viewer</title>
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
      <h1>CSCN DAG Viewer</h1>
      <p>后端 API 已就绪，但前端构建产物不存在。</p>
      <p>请在仓库根目录运行：</p>
      <p><code>cd apps/dag-viewer &amp;&amp; npm install &amp;&amp; npm run build</code></p>
      <p>然后重新刷新本页。当前预期的静态目录：</p>
      <p><code>{dist_dir}</code></p>
    </main>
  </body>
</html>
""".strip()
