**审稿意见回应大纲（基于当前结果）**

1) **认可问题**
- 单个随机细胞的 DAG 可能不代表整体（单细胞噪声与变异性）。

2) **方法更新：由“单细胞网络”转为“多细胞共识网络”**
- 对每个 *cell type × time* 聚合多细胞 DAG，保留“共识边阈值后的网络”，提高代表性与稳定性。
- 结果图（示例路径）：
  - `data/E-GEOD-93593/visualizations/kTotal/combined/Day{time}_{cell_type}.png`
  - `data/E-GEOD-93593/visualizations/kWithin/combined/Day{time}_{cell_type}.png`

3) **新增：多细胞类型对比，展示系统性差异**
- 在同一时间点并列不同 cell type 的共识网络，强调结构差异。
- 结果图（示例路径）：
  - `data/E-GEOD-93593/visualizations/kTotal/combined/Day54_*`
  - `data/E-GEOD-93593/visualizations/kWithin/combined/Day54_*`

4) **新增：细胞类型相似度矩阵（定量对比）**
- 使用 Jaccard 相似度展示 cell type 间网络差异，避免仅凭视觉判断。
- 结果图（示例路径）：
  - `data/E-GEOD-93593/visualizations/kTotal/similarity/time_54.png`
  - `data/E-GEOD-93593/visualizations/kWithin/similarity/time_54.png`

5) **保留时间维度，但改为“同一细胞类型的 time‑course 共识网络”**
- 同一 cell type 的多时间点网络变化，体现时间演化与发育轨迹。
- 结果图（示例路径）：
  - `data/E-GEOD-93593/visualizations/kTotal/time_grids/Dayall_{cell_type}.png`
  - `data/E-GEOD-93593/visualizations/kWithin/time_grids/Dayall_{cell_type}.png`

6) **补充定量证据**
- 节点数/边数/密度等统计指标，证明差异具有可重复性与统计支撑。
- 汇总表（示例路径）：
  - `data/E-GEOD-93593/visualizations/kTotal/metrics_summary.csv`
  - `data/E-GEOD-93593/visualizations/kWithin/metrics_summary.csv`

7) **结论性回应（可以写进 rebuttal）**
- 通过“共识网络 + 多细胞类型对比 + 时间序列 + 定量指标”，替代单细胞展示，系统性差异更清晰，直接回应审稿人对代表性的担忧。
