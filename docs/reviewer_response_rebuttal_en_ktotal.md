We thank the reviewer for this important point and agree that a single randomly selected cell-level CSCN may not be representative given single-cell variability. To address this, we now aggregate multiple cell-level DAGs within each *cell type × time point* and retain only consensus edges that recur above a threshold, producing **consensus CSCNs** that emphasize stable structure (Figure R1). This directly improves representativeness relative to the original single‑cell illustration.

**Figure R1. Consensus CSCNs (Day 54 examples).**  
![](../data/E-GEOD-93593/visualizations/kTotal/combined/Day54_cortical_interneuron.png)  
![](../data/E-GEOD-93593/visualizations/kTotal/combined/Day54_medial_ganglionic_eminence_neuron.png)  
![](../data/E-GEOD-93593/visualizations/kTotal/combined/Day54_telencephalic_progenitor_cell.png)

To demonstrate systematic differences across distinct cell types (rather than random variability), we compare consensus networks side‑by‑side at the same time point and quantify their structure. At Day 54, the **telencephalic progenitor** network shows higher edge count and density (47 edges, density 0.0419) than **cortical interneuron** (43 edges, density 0.0361) or **MGE neuron** (42 edges, density 0.0333), indicating reproducible structural differences across cell types. This is also reflected in the low-overlap pattern in the Jaccard similarity heatmap (Figure R2). Quantitative values are from `data/E-GEOD-93593/visualizations/kTotal/metrics_summary.csv`.
Corresponding consensus-edge CSVs:  
`data/E-GEOD-93593/visualizations/kTotal/consensus_edges/Day54_cortical_interneuron.csv`  
`data/E-GEOD-93593/visualizations/kTotal/consensus_edges/Day54_medial_ganglionic_eminence_neuron.csv`  
`data/E-GEOD-93593/visualizations/kTotal/consensus_edges/Day54_telencephalic_progenitor_cell.csv`

**Figure R2. Cell‑type similarity (Jaccard) at Day 54.**  
![](../data/E-GEOD-93593/visualizations/kTotal/similarity/time_54.png)
Underlying edge sets (per cell type, Day 54):  
`data/E-GEOD-93593/visualizations/kTotal/consensus_edges/Day54_{cell_type}.csv`  
Full merged edge table:  
`data/E-GEOD-93593/visualizations/kTotal/consensus_edges/consensus_edges_summary.csv`

We also preserve the original time‑course question by presenting **time‑grid consensus CSCNs** for each cell type. For example, in MGE neuron, the consensus network contracts over time (edges: Day26=29 → Day54=42 → Day100=17 → Day125=11), while MGE progenitor shows a late expansion (edges: Day26=24 → Day54=40 → Day100=37 → Day125=57), demonstrating lineage‑specific temporal remodeling (Figure R3). These trends are supported by `data/E-GEOD-93593/visualizations/kTotal/metrics_summary.csv`.

**Figure R3. Time‑course grids (cell types with all time points).**  
![](../data/E-GEOD-93593/visualizations/kTotal/time_grids/Dayall_medial_ganglionic_eminence_neuron.png)  
![](../data/E-GEOD-93593/visualizations/kTotal/time_grids/Dayall_medial_ganglionic_eminence_progenitor.png)
Underlying consensus-edge CSVs (MGE neuron):  
`data/E-GEOD-93593/visualizations/kTotal/consensus_edges/Day26_medial_ganglionic_eminence_neuron.csv`  
`data/E-GEOD-93593/visualizations/kTotal/consensus_edges/Day54_medial_ganglionic_eminence_neuron.csv`  
`data/E-GEOD-93593/visualizations/kTotal/consensus_edges/Day100_medial_ganglionic_eminence_neuron.csv`  
`data/E-GEOD-93593/visualizations/kTotal/consensus_edges/Day125_medial_ganglionic_eminence_neuron.csv`  
Underlying consensus-edge CSVs (MGE progenitor):  
`data/E-GEOD-93593/visualizations/kTotal/consensus_edges/Day26_medial_ganglionic_eminence_progenitor.csv`  
`data/E-GEOD-93593/visualizations/kTotal/consensus_edges/Day54_medial_ganglionic_eminence_progenitor.csv`  
`data/E-GEOD-93593/visualizations/kTotal/consensus_edges/Day100_medial_ganglionic_eminence_progenitor.csv`  
`data/E-GEOD-93593/visualizations/kTotal/consensus_edges/Day125_medial_ganglionic_eminence_progenitor.csv`

Finally, the consensus‑node GO/KEGG enrichment provides biological support for these structural differences (Figure R4; Tables R5–R6). The GO summary indicates that **membrane raft/microdomain** and **PLC/G‑protein signaling** are enriched at intermediate stages (e.g., cortical interneuron Day54; MGE neuron Day26/100), while later MGE neuron time points shift toward **neurotransmitter exocytosis** (Day125: dense‑core granule exocytosis; calcium‑regulated neurotransmitter release). In MGE progenitors, the late stage is enriched for **regulation of synaptic vesicle fusion** (Day125), whereas telencephalic progenitors show **early synaptic vesicle fusion regulation** (Day26) followed by **cell‑cell adhesion** (Day54). KEGG summaries reinforce these trends, with **synaptic vesicle cycle** and **glutamatergic synapse** emerging in MGE neuron/progenitor at later time points (Day100–125), consistent with maturation of synaptic machinery. Full per‑group GO/KEGG summaries are provided in Table R5 (`data/E-GEOD-93593/visualizations/kTotal/enrichment_results_consensus/consensus_GO_summary.csv`) and Table R6 (`data/E-GEOD-93593/visualizations/kTotal/enrichment_results_consensus/consensus_KEGG_summary.csv`).

**Figure R4. GO/KEGG dotplots from consensus nodes (Day54 cortical interneuron).**  
![](../data/E-GEOD-93593/visualizations/kTotal/enrichment_results_consensus/Day54_cortical_interneuron/kTotal_Day54_cortical_interneuron_consensus_GO_dotplot.png)  
![](../data/E-GEOD-93593/visualizations/kTotal/enrichment_results_consensus/Day54_cortical_interneuron/kTotal_Day54_cortical_interneuron_consensus_KEGG_dotplot.png)
Corresponding enrichment CSVs:  
`data/E-GEOD-93593/visualizations/kTotal/enrichment_results_consensus/Day54_cortical_interneuron/kTotal_Day54_cortical_interneuron_consensus_GO_results.csv`  
`data/E-GEOD-93593/visualizations/kTotal/enrichment_results_consensus/Day54_cortical_interneuron/kTotal_Day54_cortical_interneuron_consensus_KEGG_results.csv`
