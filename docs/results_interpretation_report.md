# Results Interpretation Report (GO/KEGG)

> Scope: Interpretation is based on the **filtered specific-edge** enrichment results and uses **kTotal / kWithin only**. Each claim below cites the exact CSV file(s) that contain the supporting GO/KEGG terms.

---

## 1) Global biological theme
**Observation:** Enriched terms repeatedly align with a neurodevelopmental signaling axis:
**GPCR / PLC / PI signaling → membrane raft / microdomain organization → filopodium / morphology → synaptic & neurotransmitter transport**.

**Evidence CSVs:**

- **kTotal summary table** (fast index):
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_summary/explanation_table.csv`
- Representative GO/KEGG CSVs used below (see sections 2–5).

---

## 2) Early LGE neuron (Day26) — synaptic signaling appears early
**Key interpretation:** Day26 LGE neurons show early synaptic signaling and glutamatergic features.

**Evidence (GO / KEGG):**

- **GO (sender):** L‑glutamate import, exocytosis, filopodium
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day26_lateral_ganglionic_eminence_neuron__sender/kTotal_Day26_lateral_ganglionic_eminence_neuron_sender_GO_results.csv`
- **KEGG (sender):** Glutamatergic synapse
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day26_lateral_ganglionic_eminence_neuron__sender/kTotal_Day26_lateral_ganglionic_eminence_neuron_sender_KEGG_results.csv`

---

## 3) Mid-stage cortical interneuron (Day54) — GPCR/PLC signaling + membrane platforms
**Key interpretation:** Day54 cortical interneurons show strong membrane‑signal coupling: upstream GPCR/PLC signaling on the sender side and membrane raft organization on receiver side.

**Evidence (GO / KEGG):**

- **GO (sender):** G‑protein beta‑subunit binding / phospholipase C activity / GTPase complex
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day54_cortical_interneuron__sender/kTotal_Day54_cortical_interneuron_sender_GO_results.csv`
- **GO (receiver):** membrane raft / membrane microdomain / basal plasma membrane
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day54_cortical_interneuron__receiver/kTotal_Day54_cortical_interneuron_receiver_GO_results.csv`
- **KEGG (sender):** Chemokine signaling / Ras / Inositol phosphate metabolism (shared signaling modules)
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day54_cortical_interneuron__sender/kTotal_Day54_cortical_interneuron_sender_KEGG_results.csv`

---

## 4) Late MGE progenitor (Day125) — neurotransmitter transport & signaling transition
**Key interpretation:** Late MGE progenitors show neurotransmitter transport and signaling readiness, consistent with maturation trajectory.

**Evidence (GO / KEGG):**

- **GO (sender):** L‑glutamate import / neurotransmitter transport / amino acid transport
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day125_medial_ganglionic_eminence_progenitor__sender/kTotal_Day125_medial_ganglionic_eminence_progenitor_sender_GO_results.csv`
- **GO (receiver):** G‑protein beta‑subunit binding / phospholipase C activity
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day125_medial_ganglionic_eminence_progenitor__receiver/kTotal_Day125_medial_ganglionic_eminence_progenitor_receiver_GO_results.csv`
- **KEGG (receiver):** Ras / Chemokine signaling (shared signal modules)
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day125_medial_ganglionic_eminence_progenitor__receiver/kTotal_Day125_medial_ganglionic_eminence_progenitor_receiver_KEGG_results.csv`

---

## 5) Late neural progenitor (Day125) — morphology & membrane domains
**Key interpretation:** Late neural progenitors retain migratory/structural features (filopodium, membrane domains) while also showing signal‑complex enrichment.

**Evidence (GO / KEGG):**

- **GO (sender):** filopodium / membrane raft / membrane microdomain / heterotrimeric G‑protein complex
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day125_neural_progenitor_cell__sender/kTotal_Day125_neural_progenitor_cell_sender_GO_results.csv`

---

## 6) Sender vs Receiver directionality (functional split)
**Key interpretation:** Sender gene sets emphasize upstream signaling; receiver sets emphasize membrane structure and response.

**Evidence:**

- Sender GO terms often include **G‑protein / PLC / neurotransmitter transport / exocytosis**.
  - Example: `Day54_cortical_interneuron__sender` GO CSV
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day54_cortical_interneuron__sender/kTotal_Day54_cortical_interneuron_sender_GO_results.csv`
- Receiver GO terms often include **membrane raft / microdomain / basal membrane**.
  - Example: `Day54_cortical_interneuron__receiver` GO CSV
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day54_cortical_interneuron__receiver/kTotal_Day54_cortical_interneuron_receiver_GO_results.csv`

---

## 7) Why KEGG “disease” terms appear (not pathology)
Terms like **Alcoholism / HIV / KSHV** reflect shared signaling modules (GPCR/PLC/Ras) and do **not** indicate disease in this context.

**Evidence:**

- These KEGG terms co‑occur with GPCR/PLC signals in the same groups:
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day54_cortical_interneuron__sender/kTotal_Day54_cortical_interneuron_sender_KEGG_results.csv`
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day26_lateral_ganglionic_eminence_neuron__sender/kTotal_Day26_lateral_ganglionic_eminence_neuron_sender_KEGG_results.csv`

---

## 8) Summary statement
The enriched pathways and GO terms from the **filtered specific-edge** networks consistently map to neurodevelopmental signaling and structural modules. The sender/receiver split further supports biological directionality: **signal initiation → membrane organization → structural/synaptic outcomes**.

---

### Notes

- For kWithin, enrichment is dominated by **structural / ribosomal / adhesion** terms; this is consistent with proliferative or migratory programs but less directly tied to neuronal signaling.
- The kTotal enrichment provides the strongest evidence for signaling validity in this dataset.
