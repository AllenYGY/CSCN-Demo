from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run the local CSCN DAG Viewer API and static frontend host."
    )
    parser.add_argument("--host", default="127.0.0.1", help="Bind address for the DAG viewer.")
    parser.add_argument("--port", type=int, default=8000, help="Bind port for the DAG viewer.")
    parser.add_argument(
        "--reload",
        action="store_true",
        help="Enable auto reload for local development.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        import uvicorn
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "缺少 `uvicorn`。请先安装 DAG Viewer 依赖，再运行 `scripts/biomarker/run_dag_viewer.py`。"
        ) from exc

    uvicorn.run(
        "biomarker.dag_viewer.api:create_app",
        factory=True,
        host=args.host,
        port=args.port,
        reload=args.reload,
    )


if __name__ == "__main__":
    main()
