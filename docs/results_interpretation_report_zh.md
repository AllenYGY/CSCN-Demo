# 结果解读报告（GO/KEGG）

> 范围：基于**过滤伪基因后的特异边**富集结果，仅使用 **kTotal / kWithin**。每条结论都给出对应的 CSV 路径作为证据。

---

## 1) 总体生物学主线
**观察：** 富集条目反复指向同一条神经发育通路链条：
**GPCR / PLC / PI 信号 → 膜微区组织 → 突起/形态结构（filopodium）→ 突触与递质运输**。

**证据（总表）：**

- `data/E-GEOD-93593/visualizations/kTotal/enrichment_summary/explanation_table.csv`

---

## 2) 早期 LGE 神经元（Day26）— 突触信号早期出现
**结论：** Day26 的 LGE 神经元出现谷氨酸突触相关信号，符合神经元早期连通建立。

**证据（GO / KEGG）：**

- **GO（sender）：** L‑glutamate import / exocytosis / filopodium
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day26_lateral_ganglionic_eminence_neuron__sender/kTotal_Day26_lateral_ganglionic_eminence_neuron_sender_GO_results.csv`
- **KEGG（sender）：** Glutamatergic synapse
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day26_lateral_ganglionic_eminence_neuron__sender/kTotal_Day26_lateral_ganglionic_eminence_neuron_sender_KEGG_results.csv`

---

## 3) 中期皮层抑制性神经元（Day54）— 信号 + 膜平台联动
**结论：** Day54 cortical interneuron 表现为“上游信号 + 膜微区结构响应”。

**证据（GO / KEGG）：**

- **GO（sender）：** G‑protein beta‑subunit binding / phospholipase C activity / GTPase complex
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day54_cortical_interneuron__sender/kTotal_Day54_cortical_interneuron_sender_GO_results.csv`
- **GO（receiver）：** membrane raft / membrane microdomain / basal plasma membrane
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day54_cortical_interneuron__receiver/kTotal_Day54_cortical_interneuron_receiver_GO_results.csv`
- **KEGG（sender）：** Chemokine signaling / Ras / Inositol phosphate metabolism
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day54_cortical_interneuron__sender/kTotal_Day54_cortical_interneuron_sender_KEGG_results.csv`

---

## 4) 晚期 MGE 前体（Day125）— 递质运输与功能化
**结论：** 晚期 MGE progenitor 进入“神经递质运输 + 信号转导准备”阶段。

**证据（GO / KEGG）：**

- **GO（sender）：** L‑glutamate import / neurotransmitter transport / amino acid transport
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day125_medial_ganglionic_eminence_progenitor__sender/kTotal_Day125_medial_ganglionic_eminence_progenitor_sender_GO_results.csv`
- **GO（receiver）：** G‑protein beta‑subunit binding / phospholipase C activity
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day125_medial_ganglionic_eminence_progenitor__receiver/kTotal_Day125_medial_ganglionic_eminence_progenitor_receiver_GO_results.csv`
- **KEGG（receiver）：** Ras / Chemokine signaling
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day125_medial_ganglionic_eminence_progenitor__receiver/kTotal_Day125_medial_ganglionic_eminence_progenitor_receiver_KEGG_results.csv`

---

## 5) 晚期 neural progenitor（Day125）— 结构与膜平台仍占主导
**结论：** late neural progenitor 仍保留迁移/突起结构特征。

**证据（GO）：**

- **GO（sender）：** filopodium / membrane raft / membrane microdomain / heterotrimeric G‑protein complex
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day125_neural_progenitor_cell__sender/kTotal_Day125_neural_progenitor_cell_sender_GO_results.csv`

---

## 6) sender vs receiver 的方向性一致性
**结论：** sender 富集信号转导类通路，receiver 富集膜结构/形态响应通路，符合“信号源 → 结构执行端”。

**证据：**

- sender 例子：
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day54_cortical_interneuron__sender/kTotal_Day54_cortical_interneuron_sender_GO_results.csv`
- receiver 例子：
  - `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day54_cortical_interneuron__receiver/kTotal_Day54_cortical_interneuron_receiver_GO_results.csv`

---

## 7) KEGG “疾病通路”解释
KEGG 中的 **Alcoholism / HIV / KSHV** 等并非疾病结论，而是通用信号模块（GPCR / PLC / Ras）的再利用。

**证据：**

- `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day54_cortical_interneuron__sender/kTotal_Day54_cortical_interneuron_sender_KEGG_results.csv`
- `data/E-GEOD-93593/visualizations/kTotal/enrichment_results/Day26_lateral_ganglionic_eminence_neuron__sender/kTotal_Day26_lateral_ganglionic_eminence_neuron_sender_KEGG_results.csv`

---

## 8) 总结性表述
基于过滤后的特异边富集结果，**神经发育信号通路（GPCR/PLC/PI）→ 膜微区组织 → 突触/递质运输**的链条在多个关键组反复出现，且 sender/receiver 的功能分工方向一致。因此可以认为最终通路结果具有明确的生物学合理性和可信度。

---

### 备注

- kWithin 的富集更偏结构/增殖通路（黏附、核糖体、小亚基等），解释价值略弱于 kTotal。
- 若需要更强证据，可增加置换检验或合并 sender+receiver 进行富集。
