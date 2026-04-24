from __future__ import annotations

from pathlib import Path

import numpy as np

from ..layout import RunLayout
from ..workflow import load_prepared_run


def run_ckm(
    layout: RunLayout,
    group_key: str,
    *,
    alpha: float = 0.05,
    beta_transform: str = "log1p",
    strict: bool = True,
    save: bool = True,
) -> np.ndarray:
    from ..core import CSCN

    prepared = load_prepared_run(layout)
    if group_key not in prepared.groups:
        raise ValueError(
            f"CKM group must exist in the prepared run. Available groups: {sorted(prepared.groups)}"
        )

    cscn_path = layout.cscn_object_path(group_key)
    if not cscn_path.is_file():
        raise FileNotFoundError(f"Missing CSCN object for group {group_key}: {cscn_path}")

    cscn = CSCN.load_from_file(cscn_path)
    dags = [(idx, dag) for idx, dag in cscn.load_all_dags()]
    save_path: str | None = None
    if save:
        layout.ckm_dir.mkdir(parents=True, exist_ok=True)
        save_path = str(layout.ckm_dir / f"{group_key}_ckm.npy")

    return cscn.compute_ckm(
        dags=dags,
        alpha=alpha,
        beta_transform=beta_transform,
        save_path=save_path,
        strict=strict,
    )
