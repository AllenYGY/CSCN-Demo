#!/usr/bin/env bash
set -euo pipefail

cd /home/jovyan/work/CSCN

export PYTHONPATH="$PWD/src"
export MPLCONFIGDIR="$PWD/.mplcache"
mkdir -p "$MPLCONFIGDIR"

echo "[1/5] prepare GSE164378 shared2000 inputs"
python scripts/prep/prepare_GSE164378.py

echo "[2/5] run CSCN cases in parallel"
python -c "from cscn.cli import main; raise SystemExit(main())" run-all --config configs/GSE164378/rna_only.yaml \
  > gse164378_shared2000_rna_only.log 2>&1 &
PID_RNA=$!

python -c "from cscn.cli import main; raise SystemExit(main())" run-all --config configs/GSE164378/adt_only.yaml \
  > gse164378_shared2000_adt_only.log 2>&1 &
PID_ADT=$!

python -c "from cscn.cli import main; raise SystemExit(main())" run-all --config configs/GSE164378/rna_adt_joint.yaml \
  > gse164378_shared2000_rna_adt_joint.log 2>&1 &
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
  echo "[ERROR] At least one CSCN case failed. Skipping modality comparison." >&2
  exit 1
fi

echo "[5/5] run modality clustering comparison"
python scripts/analysis/compare_gse164378_modalities.py --umap

echo
echo "Done. Results:"
echo "  CSCN runs:"
echo "    runs/gse164378_3p_shared2000_rna_only"
echo "    runs/gse164378_3p_shared2000_adt_only"
echo "    runs/gse164378_3p_shared2000_rna_adt_joint"
echo "  Clustering comparison:"
echo "    results/gse164378_modality_clustering_shared2000"
echo
echo "Key files:"
echo "  results/gse164378_modality_clustering_shared2000/clustering_metrics.csv"
echo "  results/gse164378_modality_clustering_shared2000/cell_assignments.csv"
echo "  results/gse164378_modality_clustering_shared2000/shared_cells.csv"
echo
echo "Per-case logs:"
echo "  gse164378_shared2000_rna_only.log"
echo "  gse164378_shared2000_adt_only.log"
echo "  gse164378_shared2000_rna_adt_joint.log"
