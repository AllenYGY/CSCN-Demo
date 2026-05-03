# Cross-dataset Biomarker Enrichment Summary

This note records the current cross-dataset summary figure set built from biomarker enrichment results.

Included datasets, in fixed left-to-right order:

1. `BreastTumer`
2. `SCP259`
3. `GSE159115`
4. `GSE138852`

The summary figures use:

- biomarker GO enrichment only
- biomarker KEGG enrichment only
- no `DESeq2-only` signal in the cross-dataset main figures

## Output Files

- [GO summary PDF](/Users/allenygy/Research/CSCN/results/paper_figures/cross_dataset/cross_dataset_go_summary.pdf)
- [GO summary PNG](/Users/allenygy/Research/CSCN/results/paper_figures/cross_dataset/cross_dataset_go_summary.png)
- [KEGG summary PDF](/Users/allenygy/Research/CSCN/results/paper_figures/cross_dataset/cross_dataset_kegg_summary.pdf)
- [KEGG summary PNG](/Users/allenygy/Research/CSCN/results/paper_figures/cross_dataset/cross_dataset_kegg_summary.png)
- [GO mapped term table](/Users/allenygy/Research/CSCN/results/paper_figures/cross_dataset/cross_dataset_biomarker_go_terms.csv)
- [KEGG mapped term table](/Users/allenygy/Research/CSCN/results/paper_figures/cross_dataset/cross_dataset_biomarker_kegg_terms.csv)
- [Process mapping notes](/Users/allenygy/Research/CSCN/results/paper_figures/cross_dataset/cross_dataset_biomarker_process_mapping.md)

Generator script:

- [make_cross_dataset_enrichment_summary.py](/Users/allenygy/Research/CSCN/scripts/figures/make_cross_dataset_enrichment_summary.py)

## Figure Design

Each dataset occupies one column.

Each row is a manually harmonized biological-process theme rather than a raw GO/KEGG term.

Cell encoding:

- color intensity: strongest `-log10(adjusted p)` for that dataset-theme pair
- cell text: number of mapped raw terms contributing to that theme
- `–`: no mapped signal for that theme

## GO Themes

- `Secretory / barrier defense`
- `Immune / antigen presentation`
- `Sulfur / redox regulation`
- `Epithelial / adhesion remodeling`
- `ATP / mitochondrial energy`
- `Nucleotide biosynthesis`
- `Glycolysis / carbon metabolism`
- `Neuronal / synapse organization`

## KEGG Themes

- `Mitochondrial respiration / OXPHOS`
- `Neurodegeneration-associated pathways`
- `Immune / infection`
- `Antigen presentation`
- `Proteostasis / proteasome`
- `Sulfur metabolism`
- `Glycolysis / carbon metabolism`
- `Hypoxia / metabolic signaling`

## Interpretation Intent

This figure set is designed to answer a single question:

- what biological axis is emphasized by the final biomarker set in each dataset?

The intended high-level reading is:

- `BreastTumer`: ATP / nucleotide biosynthesis plus immune-antigen and proteostasis signals
- `SCP259`: secretory-barrier, immune-interface, sulfur/redox
- `GSE159115`: glycolysis / carbon metabolism and hypoxia-associated signaling
- `GSE138852`: neuronal / synapse organization

## Notes

- `GSE138852` currently has no significant biomarker KEGG result, so its KEGG column is expected to be sparse or empty.
- The cross-dataset summary is complementary to the existing single-dataset figures under:
  - `results/paper_figures/SCP259/`
  - `results/paper_figures/GSE159115/`
  - `results/paper_figures/GSE138852/`
