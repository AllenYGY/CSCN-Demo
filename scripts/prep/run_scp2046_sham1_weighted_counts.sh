#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export PYTHONPATH="${REPO_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/cscn_mpl}"
mkdir -p "${MPLCONFIGDIR}"

TMP_DIR="${TMPDIR:-/tmp}/cscn_scp2046_sham1_weighted_counts"
mkdir -p "${TMP_DIR}"
export CSCN_REPO_ROOT="${REPO_ROOT}"
export CSCN_TMP_DIR="${TMP_DIR}"
export CSCN_RUN_NAME="scp2046_sham1_weighted_counts"
export CSCN_SPATIAL_STRATEGY="weighted_counts"

python - <<'PY'
from pathlib import Path
import gzip
import os

import numpy as np
import pandas as pd
from scipy.io import mmread

base = Path(os.environ["CSCN_REPO_ROOT"])
data = base / "data" / "SCP2046"
tmp = Path(os.environ["CSCN_TMP_DIR"])
run_name = os.environ["CSCN_RUN_NAME"]
strategy = os.environ["CSCN_SPATIAL_STRATEGY"]
tmp.mkdir(parents=True, exist_ok=True)

spatial = pd.read_csv(data / "cluster" / "spatial_s1.csv", skiprows=[1])
spatial_names = spatial["NAME"].astype(str).tolist()
coord_map = {
    str(row["NAME"]): (float(row["X"]), float(row["Y"]))
    for _, row in spatial.iterrows()
}

def read_barcodes(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]

def add_suffix(barcodes: list[str], suffix: str) -> list[str]:
    if not suffix:
        return barcodes
    return [b if b.endswith(f"_{suffix}") else f"{b}_{suffix}" for b in barcodes]

def resolve_10x_triplet(root: Path):
    files = {path.name: path for path in root.iterdir() if path.is_file()}
    if {"matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"} <= set(files):
        return {
            "matrix": files["matrix.mtx.gz"],
            "features": files["features.tsv.gz"],
            "barcodes": files["barcodes.tsv.gz"],
        }
    for path in files.values():
        name = path.name
        if not name.endswith("_matrix.mtx.gz"):
            continue
        prefix = name[: -len("_matrix.mtx.gz")]
        features = files.get(f"{prefix}_features.tsv.gz")
        barcodes = files.get(f"{prefix}_barcodes.tsv.gz")
        if features is not None and barcodes is not None:
            return {"matrix": path, "features": features, "barcodes": barcodes}
    return None

candidates = []
for root, _, files in os.walk(data / "expression"):
    root = Path(root)
    triplet = resolve_10x_triplet(root)
    if triplet is not None:
        candidates.append({"root": root, "triplet": triplet})

if not candidates:
    raise SystemExit("No 10x expression folder found under data/SCP2046/expression")

best = None
for candidate in candidates:
    root = candidate["root"]
    triplet = candidate["triplet"]
    raw_barcodes = read_barcodes(triplet["barcodes"])
    for suffix in ("1", "2", "3", "4", ""):
        cell_ids = add_suffix(raw_barcodes, suffix)
        overlap = len(set(cell_ids) & set(spatial_names))
        if best is None or overlap > best["overlap"]:
            best = {
                "root": root,
                "triplet": triplet,
                "suffix": suffix,
                "overlap": overlap,
                "cell_ids": cell_ids,
            }

if best is None or best["overlap"] == 0:
    raise SystemExit("Could not match any expression folder to spatial_s1.csv")

features = pd.read_csv(best["triplet"]["features"], sep="\t", header=None)
if features.shape[1] >= 2:
    genes = features.iloc[:, 1].astype(str).fillna("")
    genes = genes.where(genes != "", features.iloc[:, 0].astype(str))
else:
    genes = features.iloc[:, 0].astype(str)
genes = genes.tolist()

with gzip.open(best["triplet"]["matrix"], "rb") as handle:
    mat = mmread(handle).tocsr()

id_to_col = {cid: i for i, cid in enumerate(best["cell_ids"])}
ordered_cell_ids = [cid for cid in spatial_names if cid in id_to_col]
ordered_cols = [id_to_col[cid] for cid in ordered_cell_ids]
if not ordered_cell_ids:
    raise SystemExit("No overlapping cell ids between expression and spatial_s1.csv")

k_neighbors = max(1, int(np.ceil(len(ordered_cell_ids) * 0.2)))
print(f"[SCP2046 weighted_counts] aligned cells={len(ordered_cell_ids)} k_neighbors={k_neighbors}")

mat = mat[:, ordered_cols]
mean = np.asarray(mat.mean(axis=1)).ravel()
mean_sq = np.asarray(mat.power(2).mean(axis=1)).ravel()
var = mean_sq - mean**2
top_idx = np.argsort(var)[::-1][:min(150, mat.shape[0])]
top_genes = [genes[i] for i in top_idx]

seen = {}
deduped = []
for gene in top_genes:
    n = seen.get(gene, 0)
    seen[gene] = n + 1
    deduped.append(gene if n == 0 else f"{gene}__{n}")

expr = pd.DataFrame(
    mat[top_idx, :].toarray(),
    index=deduped,
    columns=ordered_cell_ids,
)
expr.index.name = "GENE"
expr_path = tmp / "expression_s1_top150.csv.gz"
expr.to_csv(expr_path, compression="gzip")

meta = pd.DataFrame(
    {
        "cell_id": ordered_cell_ids,
        "group": ["sham1"] * len(ordered_cell_ids),
        "X": [coord_map[cid][0] for cid in ordered_cell_ids],
        "Y": [coord_map[cid][1] for cid in ordered_cell_ids],
    }
)
meta_path = tmp / "metadata_s1.csv"
meta.to_csv(meta_path, index=False)

config_path = tmp / "config.yaml"
config_path.write_text(
    f"""run_name: {run_name}

input:
  format: tables
  expr_path: {expr_path}
  metadata_path: {meta_path}
  expr_orientation: genes_by_cells
  gene_key: GENE
  metadata_cell_id_column: cell_id
  obs_group_key: group
  spatial_x_key: X
  spatial_y_key: Y

preprocess:
  normalize: false
  log1p: false
  sample_per_group: null
  random_seed: 42
  gene_selection:
    top_n: 150

run:
  output_dir: {base / "runs"}
  max_workers: 4
  sigmoid_score: 0.1
  significance_level: 0.05
  max_cond_vars: 10
  use_bitmap: true
  using_nmf: false
  show_progress: false
  progress_interval: 100
  spatial:
    enabled: true
    strategy: {strategy}
    mode: knn
    k: {k_neighbors}
    kernel: gaussian
    lambda_expr: 0.0
    min_effective_neighbors: 15

aggregate:
  consensus: true
  consensus_threshold_mode: auto

biomarker:
  enabled: false
""",
    encoding="utf-8",
)
print(config_path)
PY

CONFIG_PATH="${TMP_DIR}/config.yaml"
python -c "from cscn.cli import main; raise SystemExit(main())" prepare --config "${CONFIG_PATH}"
python -c "from cscn.cli import main; raise SystemExit(main())" run --config "${CONFIG_PATH}"
python -c "from cscn.cli import main; raise SystemExit(main())" aggregate --config "${CONFIG_PATH}"

echo "Done: ${REPO_ROOT}/runs/${CSCN_RUN_NAME}"
