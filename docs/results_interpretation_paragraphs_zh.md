# 结果解读（论文段落版，精炼）

基于过滤伪基因后的特异边富集结果，我们观察到一条稳定且符合神经发育规律的功能链：**GPCR/PLC/PI 信号转导 → 膜微区/膜筏组织 → 神经突起与突触结构 → 递质运输/突触功能**。这一主线在多个关键组反复出现，且在 sender 与 receiver 之间体现出明确的功能分工，支持通路合理性与方向性。证据可见于 kTotal 的汇总解释表：`data/E-GEOD-93593/visualizations/kTotal/enrichment_summary/explanation_table.csv`。

在早期阶段（Day26），LGE neuron 的 sender 基因集显著富集 **L‑glutamate import、exocytosis、filopodium**，并在 KEGG 中出现 **Glutamatergic synapse**，提示神经元早期即开始构建突触相关信号与递质外排能力（`.../Day26_lateral_ganglionic_eminence_neuron__sender/..._GO_results.csv` 与 `..._KEGG_results.csv`）。同一时期的 telencephalic progenitor 也出现 **PLC/G‑protein 相关信号** 与 **synaptic vesicle 组织**相关条目，表明前体细胞已具备突触组织的早期分子基础（`.../Day26_telencephalic_progenitor_cell__sender/..._GO_results.csv` 与 `...__receiver/..._GO_results.csv`）。

中期阶段（Day54），cortical interneuron 的 sender 富集 **G‑protein beta‑subunit binding / phospholipase C activity / GTPase complex**，而 receiver 富集 **membrane raft / membrane microdomain / basal plasma membrane**，呈现“信号进入—膜平台组织—结构响应”的典型分工，这与中期抑制性神经元的信号整合和突起塑形一致（`.../Day54_cortical_interneuron__sender/..._GO_results.csv` 与 `...__receiver/..._GO_results.csv`）。

晚期阶段（Day125），MGE progenitor 的 sender 显著富集 **neurotransmitter transport / L‑glutamate import**，receiver 则富集 **G‑protein/PLC 信号**，提示其网络正向成熟神经元功能过渡（`.../Day125_medial_ganglionic_eminence_progenitor__sender/..._GO_results.csv` 与 `...__receiver/..._GO_results.csv`）。同时 neural progenitor 的 sender 富集 **filopodium / membrane raft / G‑protein complex**，说明其仍保留前体迁移与形态动态特征（`.../Day125_neural_progenitor_cell__sender/..._GO_results.csv`）。

需要说明的是，KEGG 中出现的“Alcoholism / HIV / KSHV”等条目不代表疾病结论，而是由于这些通路包含 GPCR/PLC/Ras 等通用信号基因，在神经发育背景下可视为共享信号模块（例如 `.../Day54_cortical_interneuron__sender/..._KEGG_results.csv`）。综合上述，最终通路结果具有明确的神经发育生物学合理性，并且 sender/receiver 的方向性分工支持网络结构的可解释性。
