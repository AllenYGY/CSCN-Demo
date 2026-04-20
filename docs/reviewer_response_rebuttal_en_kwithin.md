We thank the reviewer for this important point and agree that a single randomly selected cell-level CSCN may not be representative given single-cell variability. In the revised analysis, we replaced the single-cell illustration with **consensus CSCNs** computed for each *cell type × time point* by aggregating multiple cell-level DAGs and retaining edges that recur above a consensus threshold (thus emphasizing stable, reproducible structure). This yields more robust networks and directly addresses representativeness. To demonstrate systematic differences across distinct cell types, we now compare consensus networks side-by-side at the same time point (e.g., Day 54 across multiple cell types) and quantify inter–cell type similarity using Jaccard indices of edge sets. We also preserve the original time-course question by presenting **time-grid consensus CSCNs** for each cell type, enabling direct visualization of temporal evolution within the same lineage. Quantitative summaries of node/edge counts and density accompany these panels to support the visual comparisons. Together, these additions show that differences are systematic across cell types rather than artifacts of a single cell, and they provide a stronger, more representative demonstration than the original Figure 3C.

Referenced figures/tables (example renderings, kWithin):

- Consensus networks per time × cell type (Day 54 examples)
![](../data/E-GEOD-93593/visualizations/kWithin/combined/Day54_cortical_interneuron.png)
![](../data/E-GEOD-93593/visualizations/kWithin/combined/Day54_medial_ganglionic_eminence_neuron.png)
![](../data/E-GEOD-93593/visualizations/kWithin/combined/Day54_telencephalic_progenitor_cell.png)

- Time-course grids per cell type
![](../data/E-GEOD-93593/visualizations/kWithin/time_grids/Dayall_cortical_interneuron.png)

- Cell-type similarity (Jaccard) per time point
![](../data/E-GEOD-93593/visualizations/kWithin/similarity/time_54.png)

- GO/KEGG dotplots from consensus nodes (example: Day54 cortical interneuron)
![](../data/E-GEOD-93593/visualizations/kWithin/enrichment_results_consensus/Day54_cortical_interneuron/kWithin_Day54_cortical_interneuron_consensus_GO_dotplot.png)
![](../data/E-GEOD-93593/visualizations/kWithin/enrichment_results_consensus/Day54_cortical_interneuron/kWithin_Day54_cortical_interneuron_consensus_KEGG_dotplot.png)

- Quantitative metrics summary (table)
`data/E-GEOD-93593/visualizations/kWithin/metrics_summary.csv`
