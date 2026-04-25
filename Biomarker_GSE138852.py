from __future__ import annotations

import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


def main() -> int:
    from cscn.cli import main as cscn_main

    print(
        "Biomarker_GSE138852.py is now a compatibility shim. "
        "Use `cscn run-all --config configs/GSE138852/config.yaml` for new runs.",
        file=sys.stderr,
    )
    return cscn_main(["run-all", "--config", "configs/GSE138852/config.yaml"])


if __name__ == "__main__":
    raise SystemExit(main())
