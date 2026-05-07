#!/usr/bin/env bash
set -euo pipefail

cd /Users/allenygy/Research/CSCN

export PYTHONPATH="$PWD/src"
export MPLCONFIGDIR="$PWD/.mplcache"
mkdir -p "$MPLCONFIGDIR"

echo "[1/4] prepare GSE128639 MNC inputs (including shared3000 subset)"
./.venv/bin/python scripts/prep/prepare_GSE128639.py

echo "[2/4] run shared3000 CSCN cases in parallel"
./.venv/bin/python -c "from cscn.cli import main; raise SystemExit(main())" run-all --config configs/GSE128639/shared3000_rna_only.yaml \
  > gse128639_shared3000_rna_only.log 2>&1 &
PID_RNA=$!

./.venv/bin/python -c "from cscn.cli import main; raise SystemExit(main())" run-all --config configs/GSE128639/shared3000_adt_only.yaml \
  > gse128639_shared3000_adt_only.log 2>&1 &
PID_ADT=$!

./.venv/bin/python -c "from cscn.cli import main; raise SystemExit(main())" run-all --config configs/GSE128639/shared3000_rna_adt_joint.yaml \
  > gse128639_shared3000_rna_adt_joint.log 2>&1 &
PID_JOINT=$!

echo "  RNA-only PID:      $PID_RNA"
echo "  ADT-only PID:      $PID_ADT"
echo "  RNA+ADT joint PID: $PID_JOINT"

FAIL=0
for LABEL in RNA ADT JOINT; do
  PID_VAR="PID_${LABEL}"
  PID="${!PID_VAR}"
  if ! wait "$PID"; then
    echo "[ERROR] ${LABEL} case failed. Check its log file." >&2
    FAIL=1
  fi
done

if [[ "$FAIL" -ne 0 ]]; then
  echo "[ERROR] At least one CSCN case failed." >&2
  exit 1
fi

echo "[3/4] run modality clustering comparison"
./.venv/bin/python scripts/analysis/compare_gse128639_modalities.py --umap

echo "[4/4] finished all GSE128639 shared3000 cases"
echo
echo "Done. Results:"
echo "  CSCN runs:"
echo "    runs/gse128639_mnc_shared3000_rna_only"
echo "    runs/gse128639_mnc_shared3000_adt_only"
echo "    runs/gse128639_mnc_shared3000_rna_adt_joint"
echo "  Clustering comparison:"
echo "    results/gse128639_modality_clustering_shared3000"
echo
echo "Shared subset inputs:"
echo "  data/GSE128639/cscn_inputs/gse128639_mnc_shared3000_metadata.csv.gz"
echo "  data/GSE128639/cscn_inputs/gse128639_mnc_shared3000_cells.txt"
echo
echo "Key files:"
echo "  results/gse128639_modality_clustering_shared3000/clustering_metrics.csv"
echo "  results/gse128639_modality_clustering_shared3000/cell_assignments.csv"
echo "  results/gse128639_modality_clustering_shared3000/shared_cells.csv"
echo
echo "Per-case logs:"
echo "  gse128639_shared3000_rna_only.log"
echo "  gse128639_shared3000_adt_only.log"
echo "  gse128639_shared3000_rna_adt_joint.log"
