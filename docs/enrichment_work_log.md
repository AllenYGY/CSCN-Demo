# Enrichment Workflow — Detailed Summary

This log summarizes everything done **from the enrichment stage onward**, including new outputs and script changes.

## 1) Prepare enrichment inputs (Python)
**Script:** `visualize_celltype_time.py`

### Additions

- **Filtered specific edges** (pseudogene/invalid nodes removed):
  - `data/E-GEOD-93593/visualizations/{method}/specific_edges_filtered/Day{time}_{cell_type}.csv`
  - `data/E-GEOD-93593/visualizations/{method}/specific_edges_filtered/specific_edges_summary.csv`
- **Sender/receiver gene lists** derived from filtered specific edges:
  - `data/E-GEOD-93593/visualizations/{method}/enrichment_inputs/Day{time}_{cell_type}__sender.csv`
  - `data/E-GEOD-93593/visualizations/{method}/enrichment_inputs/Day{time}_{cell_type}__receiver.csv`
  - `data/E-GEOD-93593/visualizations/{method}/enrichment_inputs/enrichment_inputs_summary.csv`

### Pseudogene filtering rules (current)
Filtered out genes matching patterns:

- `LOC*`, `LINC*`
- `RPLxxPxx`, `RPSxxPxx`, `MT*Pxx`
- `*P[0-9]+` (typical pseudogene suffix)
- `nan` / empty

### Rationale
These filters reduce noise in enrichment and improve biological interpretability.

---

## 2) GO / KEGG enrichment (R)
**Script:** `run_enrichment_sender_receiver.R`

### What it does

- Iterates over **sender** and **receiver** gene lists from `enrichment_inputs/`
- Converts gene IDs (SYMBOL → ENTREZ, fallback ENSEMBL → ENTREZ)
- Runs **GO (ALL)** and **KEGG** enrichment for each group
- Saves results (CSV) and dotplots (if `enrichplot` available)

### Output structure
```
data/E-GEOD-93593/visualizations/{method}/enrichment_results/
  Day{time}_{cell_type}__sender/
    {method}_Day{time}_{cell_type}_sender_GO_results.csv
    {method}_Day{time}_{cell_type}_sender_KEGG_results.csv
  Day{time}_{cell_type}__receiver/
    {method}_Day{time}_{cell_type}_receiver_GO_results.csv
    {method}_Day{time}_{cell_type}_receiver_KEGG_results.csv
```

### Common runtime notes

- Some KEGG calls can time out (remote KEGG REST API).
- Many groups have small gene sets → no significant results.
- kTotal / kWithin yielded more interpretable terms than Module_Correlation.

---

## 3) Fixing path issues for local DAGs
**Script:** `visualize_celltype_time.py`

### Problem
Saved CSCN objects referenced **absolute server paths** (e.g., `/home/jovyan/...`), which broke local runs.

### Fix

- `load_dags_for_group()` now **prefers local DAG folders**:
  `data/E-GEOD-93593/DAG/{method}/Day{time}_{cell_type}/result_*.pkl`
- Only falls back to `_cscn` if DAGs are missing.

---

## 4) Enrichment summary + visualization (Python)
**Script:** `summarize_enrichment.py`

### Purpose

- Summarize most enriched terms for **top groups** (by filtered specific edges)
- Generate sender/receiver **heatmaps** and **bar charts**

### Outputs (per method)
```
data/E-GEOD-93593/visualizations/{method}/enrichment_summary/
  explanation_table.csv
  GO_sender_heatmap.png
  GO_receiver_heatmap.png
  KEGG_sender_heatmap.png
  KEGG_receiver_heatmap.png
  GO_sender_receiver_counts.png
  KEGG_sender_receiver_counts.png
```

### Explanation table
`explanation_table.csv` columns:

- `method`
- `group`
- `role` (sender/receiver)
- `go_top_terms`
- `kegg_top_terms`

### Fix applied

- Group names used **underscores** in folders but **spaces** in earlier summary.
- Added `group_to_safe()` to normalize group names when reading files.

---

## 5) Interpretation rules applied

- Only **kTotal** and **kWithin** used for biological interpretation.
- Module_Correlation often yields noisy or metabolic/low-specificity terms.
- Sender terms interpreted as **upstream signaling**, receiver as **downstream effect/structure**.

---

## 6) Files added/updated
**New files:**

- `run_enrichment_sender_receiver.R`
- `summarize_enrichment.py`
- `docs/visualization_plan.md` (updated to reflect enrichment outputs)
- `docs/enrichment_work_log.md` (this file)

**Updated scripts:**

- `visualize_celltype_time.py` (filtered edges + enrichment inputs + DAG loading fix)

---

## 7) Commands used
**Generate filtered specific edges + enrichment inputs**
```
python visualize_celltype_time.py \
  --methods Module_Correlation,kWithin,kTotal \
  --min-edge-count 2 \
  --edge-frac 0.05 \
  --max-dags 4
```

**Run GO/KEGG enrichment (sender & receiver)**
```
Rscript run_enrichment_sender_receiver.R
```

**Summarize enrichment and plot**
```
python summarize_enrichment.py --methods kTotal,kWithin --top-groups 5 --top-terms 10
```

---

## 8) Next recommended checks

- Verify `explanation_table.csv` has non-empty GO/KEGG for top groups.
- If many are empty:
  - Increase `--top-groups` to include more candidates,
  - Lower filtering stringency for pseudogenes,
  - Aggregate sender/receiver together for larger gene sets.
