### Response to Reviewer’s Comment on Edge Directionality

**Reviewer Comment:** *The necessity of edge directionality in the inferred networks lacks experimental support or substantive discussion. Additional analyses are required to validate the importance of these directed relationships.*

**Response:**
We appreciate this concern and agree that directionality should be supported by evidence beyond a purely algorithmic claim. To address this, we performed a **functional asymmetry analysis** based on edge direction: for each cell type × time point, we split genes into **senders** (sources of outgoing edges) and **receivers** (targets of incoming edges), and then performed GO/KEGG enrichment separately for the two sets. The resulting enrichments consistently show **non‑overlapping, biologically coherent roles** for senders versus receivers, supporting the functional relevance of edge directionality (see Figures R5–R9).

At the network level, sender and receiver nodes occupy distinct topological roles, with clear separation between source‑dominated and target‑dominated genes, which is visible in the directed sender/receiver graphs (Figure R5).

For example, in **Day26 lateral ganglionic eminence neurons**, sender genes are enriched for **L‑glutamate import, neurotransmitter transport, and exocytosis**, while receiver genes are enriched for **G‑protein/PLC signaling complexes**, indicating a directional flow from neurotransmitter‑related processes (sender) to downstream signaling/response modules (receiver) (Figure R6).

Similarly, in **Day125 MGE neurons**, sender genes are enriched for **neurotransmitter transport and calcium‑dependent exocytosis**, whereas receiver genes show enrichment for **G‑protein beta‑subunit binding and phospholipase C activity**, again indicating a directional split between upstream signaling and downstream response components (Figure R7).

KEGG pathway analysis mirrors these sender/receiver splits (Figure R8), including paired examples from Day26 LGE neuron and Day125 MGE neuron, reinforcing that directionality aligns with biologically interpretable upstream‑to‑downstream relationships rather than stochastic structure.

Finally, across multiple cell types and time points, the sender vs receiver enrichment heatmaps show a consistent asymmetry pattern, supporting that directionality captures systematic functional roles (Figure R9).

**Figure R5. Sender/receiver network topology (directed, color‑coded).**  
![](../data/E-GEOD-93593/visualizations/kTotal/sender_receiver_graphs/Day26_lateral_ganglionic_eminence_neuron.png)  
![](../data/E-GEOD-93593/visualizations/kTotal/sender_receiver_graphs/Day125_medial_ganglionic_eminence_neuron.png)

**Figure R6. Sender vs Receiver GO enrichments (Day26 LGE neuron).**  
![](../data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day26_lateral_ganglionic_eminence_neuron__sender/kTotal_Day26_lateral_ganglionic_eminence_neuron_sender_GO_dotplot.png)  
![](../data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day26_lateral_ganglionic_eminence_neuron__receiver/kTotal_Day26_lateral_ganglionic_eminence_neuron_receiver_GO_dotplot.png)

**Figure R7. Sender vs Receiver GO enrichments (Day125 MGE neuron).**  
![](../data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day125_medial_ganglionic_eminence_neuron__sender/kTotal_Day125_medial_ganglionic_eminence_neuron_sender_GO_dotplot.png)  
![](../data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day125_medial_ganglionic_eminence_neuron__receiver/kTotal_Day125_medial_ganglionic_eminence_neuron_receiver_GO_dotplot.png)

**Figure R8. Sender vs Receiver KEGG enrichments (Day26 LGE neuron and Day125 MGE neuron).**  
![](../data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day26_lateral_ganglionic_eminence_neuron__sender/kTotal_Day26_lateral_ganglionic_eminence_neuron_sender_KEGG_dotplot.png)  
![](../data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day26_lateral_ganglionic_eminence_neuron__receiver/kTotal_Day26_lateral_ganglionic_eminence_neuron_receiver_KEGG_dotplot.png)  
![](../data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day125_medial_ganglionic_eminence_neuron__sender/kTotal_Day125_medial_ganglionic_eminence_neuron_sender_KEGG_dotplot.png)  
![](../data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day125_medial_ganglionic_eminence_neuron__receiver/kTotal_Day125_medial_ganglionic_eminence_neuron_receiver_KEGG_dotplot.png)

**Figure R9. Global sender/receiver enrichment patterns.**  
![](../data/E-GEOD-93593/visualizations/kTotal/enrichment_summary/GO_sender_heatmap.png)  
![](../data/E-GEOD-93593/visualizations/kTotal/enrichment_summary/GO_receiver_heatmap.png)

In summary, while we do not claim experimental validation of directionality, the **sender/receiver enrichment asymmetry** provides functional evidence that directed edges capture biologically meaningful roles (signal dispatch vs signal reception/response), addressing the reviewer’s concern.
