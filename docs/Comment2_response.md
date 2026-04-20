
### Response to Reviewer’s Comment on Figure 3C

**Reviewer Comment:** *In Figure 3C, the depicted single-cell causal network may not be representative due to inherent data variability. A more compelling demonstration would involve comparing networks across multiple distinct cell types to illustrate systematic differences.*

**Response:**
We thank the reviewer for this insightful point. We agree that a single randomly selected cell-level CSCN may be subject to stochastic noise and may not fully represent the collective regulatory logic of a functional cell population. To address this, we have implemented a **Consensus CSCN framework**.
The dataset includes multiple neural lineages (e.g., cortical interneuron, MGE neuron/progenitor, telencephalic progenitor, etc.), enabling cross–cell-type comparisons. Some lineages are absent at specific time points; our comparisons are therefore made within the available groups.

Instead of visualizing individual cells, we now aggregate multiple cell-level DAGs within each *cell type × time point* and retain only high-confidence edges that recur above a predefined (empirical) consensus threshold. This produces **consensus CSCNs** that emphasize the stable, reproducible regulatory backbone of the population, directly improving representativeness (Figure R1).

**Figure R1. Consensus CSCNs at Day 54.** The following figures illustrate the aggregated causal structures for three major lineages. Each network is a consensus of hundreds of individual cell DAGs (e.g., $N=214$ for interneurons, $N=166$ for MGE neurons), ensuring that the displayed topology is robust to individual-cell variability.
![](../data/E-GEOD-93593/visualizations/kTotal/combined/Day54_cortical_interneuron.png)  
![](../data/E-GEOD-93593/visualizations/kTotal/combined/Day54_medial_ganglionic_eminence_neuron.png)  
![](../data/E-GEOD-93593/visualizations/kTotal/combined/Day54_telencephalic_progenitor_cell.png)

To directly address the request for "compelling demonstration of systematic differences," we performed a side-by-side comparison of these consensus networks. At Day 54, the **telencephalic progenitor** network shows significantly higher edge density (47 edges, density 0.0419) compared to **cortical interneurons** (43 edges, density 0.0361) or **MGE neurons** (42 edges, density 0.0333).

This structural divergence is quantitatively validated by our **Jaccard similarity analysis** across all four time points (Figures R2a–R2d). Across Day 26, Day 54, Day 100, and Day 125, the heatmaps consistently show low overlap between distinct cell types, indicating that CSCN captures lineage‑specific causal hierarchies rather than random noise.
At **Day 26**, cross‑type similarity is uniformly low, with most pairwise overlaps near zero.  
At **Day 54**, the spread widens slightly but remains low overall, indicating persistent cell‑type separation.  
At **Day 100**, low similarity is again dominant, with only a few pairs showing moderate overlap.  
At **Day 125**, overlaps remain minimal across the panel, reinforcing that the directed networks are cell‑type specific rather than stochastic.

**Figure R2a. Cell-type similarity (Jaccard) at Day 26.**  
![](../data/E-GEOD-93593/visualizations/kTotal/similarity/time_26.png)

**Figure R2b. Cell-type similarity (Jaccard) at Day 54.**  
![](../data/E-GEOD-93593/visualizations/kTotal/similarity/time_54.png)

**Figure R2c. Cell-type similarity (Jaccard) at Day 100.**  
![](../data/E-GEOD-93593/visualizations/kTotal/similarity/time_100.png)

**Figure R2d. Cell-type similarity (Jaccard) at Day 125.**  
![](../data/E-GEOD-93593/visualizations/kTotal/similarity/time_125.png)

We also preserved the temporal dimension by presenting **time-grid consensus CSCNs** (Figure R3). This longitudinal view shows that networks do not fluctuate randomly but evolve in a lineage-specific manner. For example, the **MGE neuron** consensus network exhibits a "contraction" pattern over time (decreasing from 42 edges at Day 54 to 11 edges by Day 125), while the **MGE progenitor** network shows a late-stage expansion (increasing from 24 edges at Day 26 to 57 edges by Day 125), reflecting systematic temporal remodeling during maturation.

**Figure R3. Time-course grids for MGE lineages.** ![](../data/E-GEOD-93593/visualizations/kTotal/time_grids/Dayall_medial_ganglionic_eminence_neuron.png)  
![](../data/E-GEOD-93593/visualizations/kTotal/time_grids/Dayall_medial_ganglionic_eminence_progenitor.png)

Finally, the biological relevance of these systematic structural differences is supported by **GO/KEGG enrichment analysis** of the consensus nodes (Figure R4). For instance, in **cortical interneurons (Day 54)**, we observed a concentration of causal edges related to **membrane raft/microdomain** and **G-protein signaling**. In contrast, later stages of **MGE neurons (Day 125)** shift toward **neurotransmitter exocytosis** and **synaptic vesicle cycle** pathways, consistent with the functional maturation of synaptic machinery. These findings demonstrate that the differences identified by CSCN are biologically meaningful and systemically organized across cell types and time.

**Figure R4. GO/KEGG dotplots from consensus nodes (Day 54 cortical interneuron).**  
![](../data/E-GEOD-93593/visualizations/kTotal/enrichment_results_consensus/Day54_cortical_interneuron/kTotal_Day54_cortical_interneuron_consensus_GO_dotplot.png)  
![](../data/E-GEOD-93593/visualizations/kTotal/enrichment_results_consensus/Day54_cortical_interneuron/kTotal_Day54_cortical_interneuron_consensus_KEGG_dotplot.png)

**Figure R5. GO/KEGG dotplots from consensus nodes (Day 125 MGE neuron).**  
![](../data/E-GEOD-93593/visualizations/kTotal/enrichment_results_consensus/Day125_medial_ganglionic_eminence_neuron/kTotal_Day125_medial_ganglionic_eminence_neuron_consensus_GO_dotplot.png)  
![](../data/E-GEOD-93593/visualizations/kTotal/enrichment_results_consensus/Day125_medial_ganglionic_eminence_neuron/kTotal_Day125_medial_ganglionic_eminence_neuron_consensus_KEGG_dotplot.png)

In summary, by shifting to a consensus-based approach supported by quantitative similarity metrics and biological enrichment, we provide a more robust and representative demonstration of the systematic causal differences across cell types, directly addressing the reviewer's concerns.
