# Biomarker Summary

This note summarizes the biomarker outputs currently present under `data/` and evaluates whether the biomarker signals and enrichment results are biologically consistent with each dataset's comparison design.

Scope:

- `SCP259`
- `GSE115978`
- `GSE121893`
- `GSE131907`
- `GSE132465`
- `GSE138852`
- `GSE159115`

Evaluation basis:

- biomarker CSVs under each dataset directory
- current `enrichment_results/` outputs already generated under `data/`
- comparison definitions encoded in the dataset-specific biomarker scripts

## Overall Ranking

From strongest biological consistency to weakest:

1. `GSE159115`
2. `GSE138852`
3. `SCP259`
4. `GSE132465`
5. `GSE131907`
6. `GSE121893`
7. `GSE115978`

## Dataset Summary Table

| Dataset | Disease / cancer | Comparison | Representative biomarkers | Biomarker count | Main enrichment signal | Assessment | Sample Size (Cells) | Gene Size |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `GSE159115` | clear cell renal cell carcinoma | `ccRCC tumor` vs `PT-B/PT-C normal` | `CD68`, `TM4SF18`, `DCDC2`, `SLC25A25`, `CRACR2B`, `PFKL` | `6` | `glycolysis`, `pyruvate metabolism`, `HIF-1`, `carbon metabolism`, `PPP` | Strong match | `100 / 100` | `150` |
| `GSE138852` | Alzheimer's disease | `AD` vs `ct` | `GPM6A`, `NTRK2`, `CTNNA2`, `SLC1A2`, `GPC5`, `NRXN1`, `SPP1` | `7` | `synapse organization`, `synapse structure/activity`, `axonogenesis`, `neuron migration` | Strong match | `3000 / 3000` | `150` |
| `SCP259` | ulcerative colitis | `Inflamed crypt/proliferative epithelial` vs `Healthy crypt/proliferative epithelial` | `TFF1`, `S100P`, `PLA2G2A`, `PI3`, `LYZ`, `REG4`, `PDIA3`, `CEACAM5` | `24` | biomarker: `peptidase inhibitor activity`, `secretory granule lumen`, `MHC assembly`, `sulfur metabolism`; DESeq2-only: `fatty acid oxidation`, `acute inflammatory response` | Strong match | `4000 / 4000` | `150` |
| `GSE132465` | colorectal cancer | `tumor epithelial` vs `normal epithelial` | `CLCA4`, `ZG16`, `PHGR1`, `GUCA2A`, `GUCA2B`, `CA7`, `ITLN1`, `MS4A12`, `LYPD8` | `20` | `chloride transport`, `nitrogen metabolism`, `mineral absorption`, metal-ion/homeostasis terms | Good match, but mainly normal epithelial program | `1000 / 1000` | `150` |
| `GSE131907` | lung cancer | `normal lung epithelial` vs `tumor malignant epithelial` | `ABCC3`, `LINC00152`, `CRABP2`, `CEACAM5`, `PMEPA1`, `ALDH1A1`, `PVRL4`, `FAM83A` | `70` | weak GO; `PPAR signaling` in DESeq2-only set | Moderate match | `2000 / 2000` | `150` |
| `GSE121893` | dilated heart failure | `dHF` vs `N` | `RGS5`, `VWF`, `IFIT3`, `ITLN1`, `MTRNR2L13` | `22`  | `apoptosis`, `complement/coagulation`, `platelet activation`, `ECM/focal adhesion` | Partial match, noisy | `500 / 500` | `150` |
| `GSE115978` | melanoma | `malignant treatment-naive` vs `post-treatment` | `ANO1`, `GZMM`, `SOX11`, `NEFL`, `NEFM`, `CHRNA1` | `46` | only `neurofilament` and `postsynaptic cytoskeleton` | Weak match | `800 / 800` | `150` |

## Per-Dataset Notes

### `GSE159115`

Comparison:

- `ccRCC tumor` vs `PT-B/PT-C normal`

Dataset context:

- Clear cell renal cell carcinoma, compared against paired normal proximal tubule-related kidney epithelium.

Interpretation:

- This is the cleanest result in the current set.
- The enrichment pattern is highly consistent with classic `ccRCC` biology:
  - glycolytic shift
  - hypoxia/HIF-axis activation
  - central carbon metabolism rewiring
- The `DESeq2-only` enrichment is especially convincing and contains the strongest disease-level signal in the whole benchmark set.

Verdict:

- `Strong match`

### `GSE138852`

Comparison:

- `AD` vs `ct`

Dataset context:

- Alzheimer's disease versus control brain samples.

Interpretation:

- The biomarker set is tightly centered on neuronal and synaptic genes.
- GO terms such as synapse regulation, synapse assembly, axonogenesis, and neuron migration are internally coherent.
- This dataset shows one of the most structured biomarker-to-enrichment mappings in the repository.

Verdict:

- `Strong match`

### `SCP259`

Comparison:

- `Inflamed crypt/proliferative epithelial` vs `Healthy crypt/proliferative epithelial`

Dataset context:

- Ulcerative colitis, focused on inflamed versus healthy proliferative crypt epithelium.

Interpretation:

- This result is strongly consistent with inflammatory epithelial remodeling in ulcerative colitis.
- The biomarker enrichment points to a coherent secretory / barrier-defense / immune-interface program:
  - peptidase and endopeptidase inhibitor activity
  - vesicle and secretory granule lumen terms
  - MHC protein complex assembly and peptide antigen assembly
  - sulfurtransferase / sulfur metabolism
- Representative biomarker genes such as `S100P`, `TIMP1`, `PI3`, `LYZ`, `CEACAM5`, `PDIA3`, and `HLA-DMA` support an inflamed secretory epithelial state rather than a generic proliferation-only signature.
- The `DESeq2-only` set complements this well:
  - `fatty acid oxidation`
  - `lipid oxidation`
  - `fatty acid beta-oxidation`
  - `acute inflammatory response`
- Taken together, the comparison reads as a shift away from epithelial metabolic homeostasis and toward inflammatory, secretory, and immune-interfacing functions.
- The `KEGG` layer is also coherent:
  - biomarker set: `sulfur metabolism`
  - `DESeq2-only` set: `fatty acid degradation`

Verdict:

- `Strong match`

### `GSE132465`

Comparison:

- `tumor epithelial` vs `normal epithelial`

Dataset context:

- Colorectal cancer epithelial cells versus normal epithelial cells.

Interpretation:

- The dominant signal is not a classic tumor-proliferation program.
- Instead, the biomarkers strongly resemble mature normal colorectal epithelial identity:
  - ion transport
  - secretory/absorptive epithelial function
  - mucosal barrier/homeostasis features
- This still makes biological sense for a tumor-vs-normal epithelial contrast, because loss of differentiated epithelial function is itself informative.

Verdict:

- `Good match`
- Best interpreted as `normal epithelial program loss`, not as a direct tumor-mechanism panel.

### `GSE131907`

Comparison:

- `normal lung epithelial` vs `tumor malignant epithelial`

Dataset context:

- Lung cancer malignant epithelial cells versus normal lung epithelial cells.

Interpretation:

- The biomarker genes themselves are fairly plausible for malignant epithelial lung programs.
- The enrichment results are weaker than the gene list:
  - biomarker GO terms are narrow and not very disease-defining
  - the stronger enrichment appears in the `DESeq2-only` set rather than the biomarker set
- This means the selected biomarkers are not wrong, but the current enrichment does not summarize them particularly well.

Verdict:

- `Moderate match`

### `GSE121893`

Comparison:

- `dHF` vs `N`

Dataset context:

- Dilated heart failure versus normal control, with a mixed heart-cell composition background.

Interpretation:

- There are some reasonable cardiovascular/stress-associated clues:
  - `RGS5`
  - `VWF`
  - apoptosis-related terms
  - complement/coagulation and adhesion-related KEGG terms
- However, the biomarker list contains substantial noise from pseudogenes and `RP11-`/`MTRNR2L` entries.
- The enrichment is therefore directionally plausible, but not clean enough to treat as a strong disease signature.

Verdict:

- `Partial match`
- Usable with caution, but not one of the strongest showcase datasets.

### `GSE115978`

Comparison:

- `malignant treatment-naive` vs `post-treatment`

Dataset context:

- Melanoma malignant cells before versus after treatment.

Interpretation:

- The current biomarker set leans toward a neural-like or dedifferentiation-like pattern:
  - `SOX11`
  - `NEFL`
  - `NEFM`
  - `CHRNA1`
  - protocadherin-family genes
- That is not fully incompatible with melanoma state changes, but the enrichment is too sparse.
- The biomarker enrichment only returns two GO terms, both driven by `NEFL/NEFM`, and does not strongly recover treatment-response or immune-evasion themes.

Verdict:

- `Weak match`

## Practical Recommendation

If the current biomarker results need to be prioritized for reporting or downstream interpretation, use them in this order:

- Primary showcase:
  - `GSE159115`
  - `GSE138852`
  - `SCP259`
  - `GSE132465`
- Secondary / acceptable with qualification:
  - `GSE131907`
- Use cautiously:
  - `GSE121893`
  - `GSE115978`

## Important Caveat

The current enrichment workflow performs over-representation analysis on the full biomarker list without splitting by effect direction (`ACE > 0` vs `ACE < 0`). For several datasets, especially tumor-vs-normal comparisons, this likely weakens interpretability by mixing biologically opposite signals into one enrichment run.

If these summaries are going to be used in a formal result section, the next improvement should be:

1. split biomarkers by direction
2. run enrichment separately for positive and negative biomarkers
3. optionally add ranked GSEA on the full biomarker or DESeq2 score vector
