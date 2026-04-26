# Biomarker Analysis Workflow

This document summarizes the current standard biomarker workflow used for the dataset-specific runs already stored under `data/`.

Important distinction:

- The repository now has a newer config-driven `cscn` CLI workflow.
- The biomarker results currently present under `data/` were mainly produced with the dataset-specific scripts under `scripts/prep/`, `scripts/biomarker/`, and `scripts/enrichment/`.

This document describes that practical workflow, because it is the one that matches the outputs already on disk.

## Goal

For a given two-group comparison:

1. prepare a clean per-group expression input
2. run DESeq2 to obtain a top gene list
3. build group-specific CSCN graphs
4. merge those graphs into a global causal graph
5. estimate gene-level causal effects on the group label
6. export biomarkers
7. run GO/KEGG enrichment on the biomarker set

## Standard Directory Pattern

For one dataset `data/<DATASET>/`, the common outputs are:

```text
data/<DATASET>/
  output_deseq/
    ... count/metadata tables used by DESeq2
    ... full DESeq2 results
    ... top gene list
    ... sampled group matrices (.npy)
    ... sampled cell ids
    ... top genes actually used by CSCN
  DAG/
    <run_name>/
      <group_a>/
      <group_b>/
  Biomarkers_<run_name>.csv
  enrichment_inputs/
  enrichment_results/
```

Typical key files:

- `output_deseq/*_all_genes.csv`
- `output_deseq/*_top150_genes.csv` or `*_top10_genes.csv`
- `output_deseq/<run_name>_<group>.npy`
- `output_deseq/<run_name>_<group>_sampled_cells.csv`
- `output_deseq/<run_name>_top150_genes_used.csv`
- `Biomarkers_<run_name>.csv`
- `enrichment_results/Biomarkers/*_GO_results.csv`
- `enrichment_results/Biomarkers/*_KEGG_results.csv`

## Step 1. Prepare Dataset-Specific Inputs

Purpose:

- filter cells into the two groups being compared
- aggregate raw counts for DESeq2
- produce count matrix and metadata tables for the contrast

Main scripts:

- `scripts/prep/prepare_GSE115978.py`
- `scripts/prep/prepare_GSE121893.py`
- `scripts/prep/prepare_GSE131907.py`
- `scripts/prep/prepare_GSE132465.py`
- `scripts/prep/prepare_GSE159115.R`

What this step usually produces in `output_deseq/`:

- count matrix for DESeq2
- metadata table for DESeq2
- sometimes excluded sample table
- sometimes alternative aggregation outputs depending on dataset

Examples:

- `count_matrix_by_sample.csv`
- `metadata_by_sample.csv`
- `count_matrix_<run_slug>_by_sample.csv`
- `metadata_<run_slug>_by_sample.csv`

Dataset-specific filtering is encoded inside each script. Examples:

- `GSE121893`: disease group, region, cell compartment
- `GSE131907`: normal lung epithelial vs tumor malignant epithelial
- `GSE132465`: normal epithelial vs tumor epithelial
- `GSE159115`: paired `ccRCC` tumor vs `PT-B/PT-C` normal

## Step 2. Run DESeq2 To Define Candidate Genes

Purpose:

- rank genes for the requested case-vs-control comparison
- save a top gene list used to seed CSCN

Main scripts:

- `scripts/prep/run_deseq2_GSE115978.R`
- `scripts/prep/run_deseq2_GSE121893.R`
- `scripts/prep/run_deseq2_GSE131907.R`
- `scripts/prep/run_deseq2_GSE132465.R`
- `scripts/prep/run_deseq2_GSE159115.R`

Common behavior:

- reads prepared count matrix and metadata
- runs `DESeq2`
- filters low-count genes
- writes:
  - full DESeq2 result table
  - top-N gene list used by the biomarker step

Common outputs:

- `deseq2_<contrast>_all_genes.csv`
- `deseq2_<contrast>_top150_genes.csv`

Notes:

- Most datasets use `top150`.
- Some existing runs use smaller effective lists, for example `top10`.
- The biomarker stage later records the final genes actually found in both groups as `*_top150_genes_used.csv` or similar.

## Step 3. Build Group-Specific CSCN Inputs

Purpose:

- sample cells from each group
- extract expression only for the DESeq2 top genes
- save per-group matrices used by CSCN

Main scripts:

- `scripts/biomarker/Biomarker_GSE115978.py`
- `scripts/biomarker/Biomarker_GSE121893.py`
- `scripts/biomarker/Biomarker_GSE131907.py`
- `scripts/biomarker/Biomarker_GSE132465.py`
- `scripts/biomarker/Biomarker_GSE138852.py`
- `scripts/biomarker/Biomarker_GSE159115.py`

Shared helper code:

- `src/biomarker/datasets.py`

What happens here:

1. load the top gene list from DESeq2
2. select eligible cells for each group
3. randomly sample a fixed number of cells per group
4. extract expression for the requested genes
5. save:
   - `.npy` matrix per group
   - sampled cell ids per group
   - top genes actually used

Shared output behavior from `save_prepared_inputs()`:

- `output_deseq/<dataset_or_run_name>_<group>.npy`
- `output_deseq/<dataset_or_run_name>_<group>_sampled_cells.csv`
- `output_deseq/<dataset_or_run_name>_top150_genes_used.csv`

Normalization:

- count-based datasets are typically normalized with `log1p(CPM)` inside `normalize_log1p()`
- some datasets already use normalized TPM-like matrices and skip raw-count normalization at this stage

## Step 4. Run CSCN Separately For Each Group

Purpose:

- infer group-specific causal structure over the selected gene set

What the biomarker scripts do:

- create one CSCN object per group
- run `run_core()`
- run `run_pc_concurrently()`
- save incremental DAG results and final CSCN object

Common outputs:

- `DAG/<run_name>/<group>/result_*.pkl`
- `<run_name>_<group>_cscn`

These outputs are later loaded back for biomarker identification.

## Step 5. Merge Group Graphs And Estimate Biomarkers

Purpose:

- convert graph node ids back to gene names
- compose a global graph from both groups
- evaluate each candidate gene as a causal treatment variable for the disease label

Shared helper code:

- `src/biomarker/graph_utils.py`

Core logic:

1. load all saved group DAGs
2. map node ids to gene symbols
3. compose group graphs into one global directed graph
4. add a synthetic sink node named `DISEASE`
5. for each candidate gene:
   - find confounders with `confounder_method="classic"`
   - run causal analysis against outcome `DISEASE`
   - keep the gene if `adjustment_formula != 0`

Current biomarker inclusion rule:

- a gene is exported as a biomarker if causal analysis succeeds and the estimated `adjustment_formula` is nonzero

Current output columns:

- `gene`
- `ACE`
- `n_confounders` when enabled

Current sorting:

- biomarkers are usually sorted by `abs(ACE)` descending

Main output:

- `Biomarkers_<run_name>.csv`

## Step 6. Run Functional Enrichment

Purpose:

- compare biomarker genes with the broader DESeq2-selected gene pool
- summarize whether the biomarker set captures coherent biology

Main script:

- `scripts/enrichment/GO_KEGG.R`

What it does:

1. read biomarker CSV
2. read the CSCN-used gene list
3. create two sets:
   - `Biomarkers`
   - `DESeq2_only_genes` = top genes used by CSCN minus biomarkers
4. convert gene ids
5. run:
   - GO enrichment
   - KEGG enrichment
6. save CSV tables and plots

Main outputs:

- `enrichment_inputs/Biomarkers.csv`
- `enrichment_inputs/DESeq2_only_genes.csv`
- `enrichment_results/Biomarkers/*_GO_results.csv`
- `enrichment_results/Biomarkers/*_KEGG_results.csv`
- `enrichment_results/DESeq2_only_genes/*_GO_results.csv`
- `enrichment_results/DESeq2_only_genes/*_KEGG_results.csv`
- corresponding barplots and dotplots

## Step 7. Interpret The Results

The practical review order should be:

1. inspect the biomarker CSV
2. check whether the top genes are biologically plausible for the dataset contrast
3. inspect biomarker GO/KEGG enrichment
4. compare biomarker enrichment against `DESeq2_only_genes`
5. decide whether the biomarker set is:
   - strong and coherent
   - directionally plausible but weak
   - noisy or poorly aligned

## Minimum Quality Checks

Before accepting a biomarker run, check:

- both groups have enough eligible cells
- sampled cells were successfully written
- the final `*_top150_genes_used.csv` is not unexpectedly small
- both group CSCN objects were saved
- `Biomarkers_<run_name>.csv` is not empty
- enrichment did not fail because of gene-id conversion or KEGG network errors

Warning signs:

- biomarkers dominated by pseudogenes or `RP11-*` entries
- enrichment driven by only one gene across many terms
- no overlap between biomarker biology and known disease/cell-state biology
- very small used-gene count relative to requested top genes

## Current Standard Run Order

For the legacy dataset-specific workflow, the practical order is:

1. `prepare_<dataset>`
2. `run_deseq2_<dataset>`
3. `Biomarker_<dataset>`
4. `GO_KEGG.R`
5. manual biological review

In compact form:

```text
raw data
-> prepare comparison-specific count/metadata tables
-> DESeq2 top-gene ranking
-> sample cells and extract group matrices
-> CSCN per group
-> merge group DAGs and run causal effect estimation
-> Biomarkers_<run_name>.csv
-> GO/KEGG enrichment
-> biological interpretation
```

## Recommended Next Improvement

The current enrichment workflow uses the full biomarker set without splitting by direction. For several tumor-vs-normal comparisons this mixes opposite biological signals.

The next refinement should be:

1. split biomarkers into `ACE > 0` and `ACE < 0`
2. run enrichment separately for both directions
3. optionally add ranked GSEA using the full biomarker or DESeq2 score vector

## Relation To The New CLI

For future new datasets, the preferred repository direction is:

- `cscn prepare`
- `cscn run`
- `cscn aggregate`
- `cscn biomarker`

But when documenting or reproducing the biomarker outputs currently stored in `data/`, the dataset-specific workflow above is the correct reference process.
