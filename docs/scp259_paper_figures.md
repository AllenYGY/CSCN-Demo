# SCP259 Paper Figures

This note defines the current paper-oriented figure set for the `SCP259` biomarker result.

Main-figure strategy:

- use `SCP259` as the single showcase biomarker dataset
- do not place `GSE132465`, `GSE131907`, `GSE121893`, or `GSE115978` in the main figure set
- keep those datasets in summary tables or supplementary material only

## Output Files

Generated figure outputs:

- [Figure 1 PDF](/results/paper_figures/SCP259/scp259_figure_1_biomarker_panel.pdf)
- ![Figure 1 PNG](/results/paper_figures/SCP259/scp259_figure_1_biomarker_panel.png)
- [Figure 2 PDF](/results/paper_figures/SCP259/scp259_figure_2_biomarker_enrichment.pdf)
- ![Figure 2 PNG](/results/paper_figures/SCP259/scp259_figure_2_biomarker_enrichment.png)
- [Figure 3 PDF](/results/paper_figures/SCP259/scp259_figure_3_deseq_only_enrichment.pdf)
- ![Figure 3 PNG](/results/paper_figures/SCP259/scp259_figure_3_deseq_only_enrichment.png)
- [Blueprint Notes](/results/paper_figures/SCP259/scp259_figure_blueprint.md)

Reproducible script:

- [make_biomarker_paper_figures.py](/Users/allenygy/Research/CSCN/scripts/figures/make_biomarker_paper_figures.py)

## Figure 1

Title:

- `SCP259 biomarker panel`

Panels:

- `A`: ranked biomarker effect sizes (`ACE`) with positive and negative directions separated by color
- `B`: functional category annotation of the final biomarker panel

Main message:

- CSCN yields a compact biomarker panel that remains biologically interpretable rather than diffuse.

## Figure 2

Title:

- `SCP259 biomarker enrichment`

Panels:

- `A`: biomarker GO enrichment
- `B`: biomarker KEGG enrichment

Main message:

- the CSCN biomarker set emphasizes an inflammatory secretory / barrier-defense / antigen-interface axis

## Figure 3

Title:

- `SCP259 DESeq2-only enrichment`

Panels:

- `A`: `DESeq2-only` GO enrichment
- `B`: `DESeq2-only` KEGG enrichment

Main message:

- the `DESeq2-only` background is broader and leans more toward fatty-acid oxidation and inflammatory context

## Main-Text Positioning

Recommended use in the manuscript:

- `Figure 1` should introduce the final biomarker panel
- `Figure 2` should support the biomarker-level mechanistic interpretation
- `Figure 3` should provide the DESeq2-only background contrast
- the cross-dataset comparison remains outside the main figure set and should stay in summary tables or supplementary material

## Excluded Cases

Do not include these as separate main-text figure panels:

- `GSE132465`
- `GSE131907`
- `GSE121893`
- `GSE115978`

Reason:

- they weaken the single-dataset mechanism story centered on `SCP259`
- they are better used as support for generalizability than as co-equal visual showcases
