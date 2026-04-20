# 图注版（中文）

**图1 细胞类型×时间的共识网络（time grids）**
每个面板为某一细胞类型在 Day26/54/100/125 的共识网络（由该组所有 DAG 合并得到）。缺失时间点标注 “No cells at Day X”。节点为基因，边为方向性关系，环形布局用于对比时间演化。

**图2 sender/receiver 全局网络（sender_receiver_graphs）**
基于过滤后的特异边构建的方向性网络。蓝色节点为 sender-only，绿色为 receiver-only，金色为兼具 sender/receiver，节点大小与度数成正比。用于直观展示方向性分工。

**图3 特异边频率分布（edge_hist）**
显示每个组中边出现频次的分布，横轴为边 (u→v)，纵轴为该边在组内 DAG 中出现的次数，突出最稳定的调控关系。

**图4 细胞类型差异热图（similarity）**
同一时间点的细胞类型相似度矩阵，基于共识边集合的 Jaccard 指数。颜色越深表示网络越相似，越浅表示差异越大。

**图5 sender/receiver 富集热图（enrichment_summary）**
展示特异边最多的 top groups 在 GO/KEGG 中的显著通路（-log10 p）。sender 与 receiver 分别绘制，用于比较信号源与结构响应的功能分工。

**图6 sender/receiver 富集数量对比（enrichment_summary）**
条形图展示每个 top group 的 GO/KEGG 显著条目数量（sender vs receiver），反映方向性通路富集的强弱差异。

**图7 代表性通路证据（解释表）**
`explanation_table.csv` 汇总 top groups 的 GO/KEGG 关键词（含 p 值），用于论文正文中的通路合理性说明。
