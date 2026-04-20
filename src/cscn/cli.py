from __future__ import annotations

import argparse

from .config import load_config
from .layout import RunLayout
from .viewer import create_app
from .workflow import (
    aggregate_run,
    prepare_run,
    run_all,
    run_biomarker_workflow,
    run_cscn,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generic CSCN workflows for scRNA-seq.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    for command in ("prepare", "run", "aggregate", "biomarker", "run-all", "viewer"):
        sub = subparsers.add_parser(command)
        sub.add_argument("--config", required=True, help="Path to a CSCN YAML config.")

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    config = load_config(args.config)

    if args.command == "prepare":
        prepare_run(config)
        return 0
    if args.command == "run":
        run_cscn(config)
        return 0
    if args.command == "aggregate":
        aggregate_run(config)
        return 0
    if args.command == "biomarker":
        run_biomarker_workflow(config)
        return 0
    if args.command == "run-all":
        run_all(config)
        return 0
    if args.command == "viewer":
        import uvicorn

        layout = RunLayout.from_config(config)
        base_dir = config.run.output_dir
        app = create_app(base_dir=base_dir, default_root=layout.run_dir)
        uvicorn.run(
            app,
            host=config.viewer.host,
            port=config.viewer.port,
            reload=config.viewer.reload,
        )
        return 0

    parser.error(f"Unsupported command: {args.command}")
    return 2
