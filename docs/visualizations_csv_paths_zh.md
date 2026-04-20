**可视化图表对应的 CSV / 数据源路径**

说明：下列路径中的 `{method}` 为 `kTotal` / `kWithin` / `Module_Correlation`；  
`{time}` 为 `26/54/100/125`；`{cell_type}` 为对应细胞类型（空格已替换为 `_`）。

---

**1) 共识网络图（Combined graphs）**  
图：  

- `data/E-GEOD-93593/visualizations/{method}/combined/Day{time}_{cell_type}.png`  
对应 CSV（共识边）：  
- `data/E-GEOD-93593/visualizations/{method}/consensus_edges/Day{time}_{cell_type}.csv`  
全量汇总：  
- `data/E-GEOD-93593/visualizations/{method}/consensus_edges/consensus_edges_summary.csv`

**2) Edge frequency 直方图（Edge hist）**  
图：  

- `data/E-GEOD-93593/visualizations/{method}/edge_hist/Day{time}_{cell_type}.png`  
对应 CSV：  
- 与共识边 CSV 相同：`.../consensus_edges/Day{time}_{cell_type}.csv`  
说明：`count` 列即边出现次数。

**3) Time‑course grids（同一细胞类型多时间点网格）**  
图：  

- `data/E-GEOD-93593/visualizations/{method}/time_grids/Dayall_{cell_type}.png`  
对应 CSV：  
- 该细胞类型在各时间点的共识边 CSV：  
  - `.../consensus_edges/Day26_{cell_type}.csv`  
  - `.../consensus_edges/Day54_{cell_type}.csv`  
  - `.../consensus_edges/Day100_{cell_type}.csv`  
  - `.../consensus_edges/Day125_{cell_type}.csv`

**4) 相似度热图（Jaccard similarity）**  
图：  

- `data/E-GEOD-93593/visualizations/{method}/similarity/time_{time}.png`  
对应 CSV / 数据源：  
- 计算基于同一时间点各 cell type 的共识边集合  
- 可用 `.../consensus_edges/consensus_edges_summary.csv` 过滤 `time={time}` 后按 `cell_type` 分组计算 Jaccard
汇总指标 CSV：  
- `data/E-GEOD-93593/visualizations/{method}/similarity/similarity_summary.csv`

**5) 趋势图（Edges/Density time trend）**  
图：  

- `data/E-GEOD-93593/visualizations/{method}/trends/Daytrend_{cell_type}.png`  
对应 CSV：  
- `data/E-GEOD-93593/visualizations/{method}/metrics_summary.csv`

**6) Sender/Receiver 全局图（3 色节点）**  
图：  

- `data/E-GEOD-93593/visualizations/{method}/sender_receiver_graphs/Day{time}_{cell_type}.png`  
对应 CSV：  
- `data/E-GEOD-93593/visualizations/{method}/specific_edges_filtered/Day{time}_{cell_type}.csv`  
说明：  
- sender-only = 只出现在 `from`  
- receiver-only = 只出现在 `to`  
- both = 同时出现在 `from` 与 `to`

**7) 特异边列表（非图）**  
CSV：  

- 原始：`data/E-GEOD-93593/visualizations/{method}/specific_edges/Day{time}_{cell_type}.csv`  
- 过滤伪基因：`data/E-GEOD-93593/visualizations/{method}/specific_edges_filtered/Day{time}_{cell_type}.csv`

**8) 单个 DAG 图（Individual DAGs）**  
图：  

- `data/E-GEOD-93593/visualizations/{method}/individual/Day{time}_{cell_type}/DAG_{id}.png`  
数据源：  
- `data/E-GEOD-93593/DAG/{method}/Day{time}_{cell_type}/result_{id}.pkl`

**9) Sender/Receiver GO/KEGG dotplot**  
图：  

- `data/E-GEOD-93593/visualizations/{method}/enrichment_results/Day{time}_{cell_type}__sender/{method}_Day{time}_{cell_type}_sender_GO_dotplot.png`  
- `data/E-GEOD-93593/visualizations/{method}/enrichment_results/Day{time}_{cell_type}__receiver/{method}_Day{time}_{cell_type}_receiver_KEGG_dotplot.png`  
对应 CSV：  
- GO：`.../enrichment_results/.../{method}_Day{time}_{cell_type}_{role}_GO_results.csv`  
- KEGG：`.../enrichment_results/.../{method}_Day{time}_{cell_type}_{role}_KEGG_results.csv`

**10) 共识节点 GO/KEGG dotplot（time+celltype）**  
图：  

- `data/E-GEOD-93593/visualizations/{method}/enrichment_results_consensus/Day{time}_{cell_type}/{method}_Day{time}_{cell_type}_consensus_GO_dotplot.png`  
对应 CSV：  
- GO：`.../enrichment_results_consensus/Day{time}_{cell_type}/{method}_Day{time}_{cell_type}_consensus_GO_results.csv`  
- KEGG：`.../enrichment_results_consensus/Day{time}_{cell_type}/{method}_Day{time}_{cell_type}_consensus_KEGG_results.csv`

**10.1) 共识节点 GO/KEGG 汇总表（按 cell type × time）**  
CSV：  

- GO 汇总：`data/E-GEOD-93593/visualizations/{method}/enrichment_results_consensus/consensus_GO_summary.csv`  
- KEGG 汇总：`data/E-GEOD-93593/visualizations/{method}/enrichment_results_consensus/consensus_KEGG_summary.csv`  
说明：每行包含 `time` 和 `cell_type`，并给出该组的 top terms。

**11) 富集总结热图 / 解释表**  
图：  

- `data/E-GEOD-93593/visualizations/{method}/enrichment_summary/GO_sender_heatmap.png`  
- `data/E-GEOD-93593/visualizations/{method}/enrichment_summary/KEGG_receiver_heatmap.png`  
对应 CSV：  
- `data/E-GEOD-93593/visualizations/{method}/enrichment_summary/explanation_table.csv`  
- 原始富集结果见 `.../enrichment_results/.../*_GO_results.csv` 与 `*_KEGG_results.csv`
