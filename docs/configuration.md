# CSCN Configuration

`cscn` 使用单个 YAML 文件驱动整个工作流。配置文件既决定输入数据如何解析，也决定预处理、CSCN 运行参数、共识聚合、biomarker 后处理和 viewer 启动方式。

这份文档关注两件事：

- 配置字段的真实语义和约束
- 每个字段最终如何影响 `runs/<run_name>/` 输出

如果你需要看更完整的模块设计、目录布局和二次开发约定，请同时阅读 [Design And Development](design_and_development.md)。

## Design Principles

当前配置体系遵循三个原则：

- 单一入口：所有 CLI 子命令都从同一个 YAML 读取参数，避免分散配置。
- 相对路径友好：配置中的相对路径都相对于配置文件所在目录解析，而不是相对于当前 shell 工作目录。
- staged workflow：`prepare`、`run`、`aggregate`、`biomarker`、`viewer` 共用一份配置，但每个阶段只消费自己需要的字段。

## Minimal Examples

### Minimal h5ad Example

```yaml
run_name: pbmc_demo

input:
  format: h5ad
  path: ../data/pbmc/pbmc.h5ad
  obs_group_key: cell_type

preprocess:
  normalize: true
  log1p: true
  sample_per_group: 200
  random_seed: 42
  gene_selection:
    top_n: 150

run:
  output_dir: ../runs
  groups: ["B cell", "T cell"]
  max_workers: 4

aggregate:
  consensus: true
  consensus_threshold_mode: auto

biomarker:
  enabled: false
```

### Minimal Tables Example

```yaml
run_name: generic_tables_demo

input:
  format: tables
  expr_path: ../data/demo/expression.csv
  metadata_path: ../data/demo/metadata.csv
  expr_orientation: cells_by_genes
  expr_cell_id_column: cell_id
  metadata_cell_id_column: cell_id
  obs_group_key: condition

preprocess:
  gene_selection:
    top_n: 150

run:
  output_dir: ../runs
```

## Top-Level Fields

### `run_name`

- 类型：`string`
- 默认行为：如果未提供，自动使用配置文件文件名去掉扩展名后的值
- 约束：不能为空
- 影响：
  - 最终 run 目录为 `run.output_dir/run_name`
  - 路径组件中的 `/` 和 `\` 会在运行时被替换为 `_`

示例：

```yaml
run_name: gse121893_generic
```

### `input`

定义原始数据如何加载成统一的 `cells x genes` 表达矩阵，以及如何把 metadata 与表达矩阵对齐。

#### `input.format`

- 可选值：`h5ad`、`tables`
- 必填

#### `input.path`

- 仅 `h5ad` 模式使用
- 指向 `.h5ad` 文件

#### `input.expr_path`

- 仅 `tables` 模式使用
- 指向表达矩阵文件

#### `input.metadata_path`

- 仅 `tables` 模式使用
- 指向 metadata 文件

#### `input.layer`

- 仅 `h5ad` 模式可用
- 省略时读取 `adata.X`
- 设置后读取 `adata.layers[layer]`

#### `input.gene_key`

- `h5ad` 模式：如果该列存在于 `adata.var.columns`，则用该列作为基因名
- `tables` 模式且 `expr_orientation: genes_by_cells`：作为“基因名所在列”的列名
- 省略时：
  - `h5ad` 读取 `var_names`
  - `genes_by_cells` 读取表达矩阵第一列

#### `input.obs_group_key`

- 用于定义分组运行 CSCN 的 metadata 列
- 省略时，系统只构造单个分组 `all`
- 如果配置了 `run.groups`，这里必须存在且能匹配到这些 group

#### `input.expr_orientation`

- 仅 `tables` 模式使用
- 可选值：`cells_by_genes`、`genes_by_cells`
- 默认：`cells_by_genes`

#### `input.expr_cell_id_column`

- `cells_by_genes` 模式下可选
- 如果指定且列存在，表达矩阵会用该列作为 cell id 索引
- 如果不指定，系统会尝试自动判断第一列是否像 cell id

#### `input.metadata_cell_id_column`

- `tables` 模式下用于 metadata 对齐
- 默认：`cell_id`
- 特殊行为：如果该值显式设置为空字符串，代码会尝试寻找 `Unnamed:*` 列作为 cell id 列，这主要用于兼容某些历史 CSV 导出

#### `input.expr_delimiter`

- 表达矩阵分隔符，可选
- 省略时自动推断：
  - `.tsv`、`.txt`、`.tsv.gz` 使用制表符
  - 其他默认使用逗号

#### `input.metadata_delimiter`

- metadata 分隔符，可选
- 推断规则与 `expr_delimiter` 相同

#### `input.obs_cell_id_key`

- 仅 `h5ad` 模式使用
- 如果该列存在于 `adata.obs.columns`，则用它作为 cell id
- 否则退回 `obs_names`

#### `input` 示例

```yaml
input:
  format: tables
  expr_path: ../../data/GSE121893/GSE121893_human_heart_sc_umi.csv.gz
  metadata_path: ../../data/GSE121893/GSE121893_covariates.csv.gz
  expr_orientation: genes_by_cells
  metadata_cell_id_column: cell_id
  obs_group_key: disease
```

#### `input` 阶段内实际发生的事情

无论输入格式是什么，系统最终都会构造一个统一的 `LoadedDataset`：

- `expression`：`pandas.DataFrame`，形状为 `cells x genes`
- `metadata`：索引与 `expression.index` 对齐
- `cell_ids`：按表达矩阵顺序排列的 cell id
- `gene_names`：去重后的基因名

额外规则：

- 重名基因会被重命名为 `GENE`, `GENE__1`, `GENE__2` 这种形式
- 表达矩阵中的非数值内容会被转成数值，无法解析的值会被填为 `0.0`
- 只有同时出现在表达矩阵和 metadata 中的 cell id 会被保留

### `preprocess`

定义 `prepare` 阶段对原始矩阵的处理方式。

#### `preprocess.normalize`

- 类型：`bool`
- 默认：`true`
- 行为：按细胞做 library-size normalization，目标总量固定为 `1e6`

公式：

```text
normalized = matrix / row_sum * 1e6
```

#### `preprocess.log1p`

- 类型：`bool`
- 默认：`true`
- 行为：在 normalization 之后做 `np.log1p`

#### `preprocess.sample_per_group`

- 类型：`int | null`
- 默认：`null`
- 行为：对每个 group 随机无放回采样固定数量细胞
- 约束：必须为正整数；如果某个 group 细胞数不足，会直接报错

#### `preprocess.random_seed`

- 类型：`int`
- 默认：`42`
- 用于 group 内采样复现

#### `preprocess.gene_selection.strategy`

- 类型：`string`
- 默认：`variance`
- 现状：当前代码会保留该字段，但实际只实现了“方差排序”这一路径

#### `preprocess.gene_selection.top_n`

- 类型：`int`
- 默认：`150`
- 行为：按全局方差排序后取前 `N` 个基因

#### `preprocess.gene_selection.gene_list_path`

- 类型：`path | null`
- 行为：如果提供，优先使用外部 gene list，而不是 `top_n`
- 匹配规则：只保留同时存在于载入数据中的基因
- 失败条件：如果一个都匹配不上，会直接报错

#### `preprocess` 示例

```yaml
preprocess:
  normalize: true
  log1p: true
  sample_per_group: 500
  random_seed: 42
  gene_selection:
    top_n: 150
```

#### `prepare` 阶段输出

执行 `cscn prepare` 后，会产出：

- `inputs/genes.csv`：节点索引到基因名的映射
- `inputs/groups.csv`：分组清单及每组细胞数
- `inputs/cell_metadata.csv`：采样后保留的 metadata
- `inputs/run_summary.json`：简要统计
- `matrices/<group>.npy`：每组用于 CSCN 的浮点矩阵
- `matrices/<group>_cells.csv`：矩阵行顺序对应的 cell id
- `config.snapshot.yaml`：配置快照

### `run`

定义 CSCN 主算法运行参数。

#### `run.output_dir`

- 类型：`path`
- 默认：`runs`
- 所有 run 输出的根目录

#### `run.groups`

- 类型：`list[str] | null`
- 省略时表示运行所有分组
- 设置后只运行指定分组
- 如果配置的分组在 metadata 中不存在，会报错

#### `run.max_workers`

- 类型：`int | null`
- 省略时自动取 `min(8, os.cpu_count())`
- 用于 `run_pc_concurrently`

#### `run.sigmoid_score`

- 类型：`float`
- 默认：`0.1`
- 直接透传给底层 `CSCN` 对象

#### `run.significance_level`

- 类型：`float`
- 默认：`0.01`

#### `run.max_cond_vars`

- 类型：`int`
- 默认：`20`

#### `run.use_bitmap`

- 类型：`bool`
- 默认：`true`

#### `run.using_nmf`

- 类型：`bool`
- 默认：`false`
- 调用 `cscn.run_core(..., usingNMF=...)` 时传入

#### `run.show_progress`

- 类型：`bool`
- 默认：`false`
- 控制底层 CSCN 是否输出内部进度

#### `run.progress_interval`

- 类型：`int`
- 默认：`100`
- 同时传给 `run_pc_concurrently` 和底层 CSCN 进度显示

#### `run` 示例

```yaml
run:
  output_dir: ../../runs
  groups: ["N", "dHF"]
  max_workers: 8
  sigmoid_score: 0.1
  significance_level: 0.05
  max_cond_vars: 20
  use_bitmap: true
  using_nmf: false
  show_progress: false
  progress_interval: 100
```

#### `run` 阶段输出

执行 `cscn run` 后，会为每个 group 产出：

- `dags/<group>/result_<i>.pkl`
- `objects/<group>_cscn.pkl`

这里的 `result_<i>.pkl` 与 `matrices/<group>_cells.csv` 的第 `i` 行一一对应。

### `aggregate`

定义单细胞 DAG 如何汇总成 group-level 共识图。

#### `aggregate.consensus`

- 类型：`bool`
- 默认：`true`
- 如果为 `false`，`aggregate` 阶段只会写一个空的 `summary.json`，不会生成 group 共识边文件

#### `aggregate.consensus_threshold_mode`

- 类型：`"auto" | int`
- 默认：`auto`

`auto` 的当前定义：

```text
max(2, ceil(0.05 * num_dags))
```

解释：

- 至少要求某条边在 2 个 DAG 中出现
- DAG 数很多时，用 5% 作为阈值

#### `aggregate` 输出

执行 `cscn aggregate` 后，会产出：

- `consensus/<group>.csv`
- `consensus/summary.json`

`consensus/<group>.csv` 字段为：

- `from`
- `to`
- `count`
- `threshold`
- `num_dags`

### `biomarker`

定义可选的 biomarker/causal 后处理。

#### `biomarker.enabled`

- 类型：`bool`
- 默认：`false`
- 仅 `run-all` 会根据该开关决定是否自动执行 biomarker
- 单独执行 `cscn biomarker` 时，不会检查这个开关

#### `biomarker.case_group`

- 类型：`string | null`
- 含义：病例组

#### `biomarker.control_group`

- 类型：`string | null`
- 含义：对照组

#### `biomarker.outcome_column`

- 类型：`string`
- 默认：`DISEASE`
- 也会作为 sink node 名称传给 biomarker 逻辑

#### `biomarker.confounder_method`

- 类型：`string`
- 默认：`classic`

#### `biomarker` 限制

当前实现只支持明确的二组比较：

- `case_group` 和 `control_group` 都必须提供
- 两者不能相同
- 两者都必须存在于 prepared run 的 group 集合中

#### `biomarker` 输出

- `biomarker/biomarkers.csv`
- `biomarker/causal_summary.json`

### `viewer`

定义 `cscn viewer` 启动 FastAPI 服务时使用的参数。

#### `viewer.enabled`

- 当前仅作为配置声明字段保留
- `cscn viewer` 命令本身不会检查它

#### `viewer.host`

- 默认：`127.0.0.1`

#### `viewer.port`

- 默认：`8000`

#### `viewer.reload`

- 默认：`false`
- 传给 `uvicorn.run(..., reload=...)`

## Stage Consumption Matrix

不同命令只使用配置的一部分：

| Command | Uses |
| --- | --- |
| `cscn prepare` | `run_name`, `input.*`, `preprocess.*`, `run.output_dir`, `run.groups` |
| `cscn run` | `run_name`, `run.output_dir`, `run.max_workers`, `run.sigmoid_score`, `run.significance_level`, `run.max_cond_vars`, `run.use_bitmap`, `run.using_nmf`, `run.show_progress`, `run.progress_interval` |
| `cscn aggregate` | `run_name`, `run.output_dir`, `aggregate.*` |
| `cscn biomarker` | `run_name`, `run.output_dir`, `biomarker.*` |
| `cscn run-all` | 上述全部 |
| `cscn viewer` | `run.output_dir`, `viewer.*`, `run_name` |

## Path Resolution Rules

配置中的相对路径都基于“配置文件所在目录”解析。

例如配置位于：

```text
examples/GSE121893/config.yaml
```

那么：

```yaml
run:
  output_dir: ../../runs
```

最终会解析到仓库根目录下的 `runs/`，而不是当前 shell 所在目录。

## Validation and Failure Modes

常见报错来源：

- `input.format` 不是 `h5ad` 或 `tables`
- `tables` 模式缺少 `expr_path` 或 `metadata_path`
- `metadata` 里找不到 `metadata_cell_id_column`
- 表达矩阵与 metadata 没有共享 cell id
- `obs_group_key` 不存在
- `run.groups` 中指定的 group 不存在
- `sample_per_group` 大于某个 group 的细胞数
- `gene_list_path` 中请求的基因一个都不在数据里
- `aggregate.consensus_threshold_mode` 既不是 `auto` 也不是整数

## Recommended Templates

### 小型 `h5ad` 数据集

```yaml
run_name: small_h5ad_demo

input:
  format: h5ad
  path: ../data/demo/demo.h5ad
  obs_group_key: cell_type

preprocess:
  normalize: true
  log1p: true
  sample_per_group: 200
  gene_selection:
    top_n: 150

run:
  output_dir: ../runs
  max_workers: 4

aggregate:
  consensus: true
  consensus_threshold_mode: auto

biomarker:
  enabled: false

viewer:
  enabled: true
  host: 127.0.0.1
  port: 8000
  reload: false
```

### 历史 CSV/GZ 数据集

```yaml
run_name: legacy_table_demo

input:
  format: tables
  expr_path: ../data/legacy/expression.csv.gz
  metadata_path: ../data/legacy/metadata.csv.gz
  expr_orientation: genes_by_cells
  metadata_cell_id_column: ""
  obs_group_key: disease

preprocess:
  normalize: true
  log1p: true
  sample_per_group: 500
  random_seed: 42
  gene_selection:
    top_n: 150

run:
  output_dir: ../runs
  groups: ["control", "case"]
  max_workers: 8

aggregate:
  consensus: true
  consensus_threshold_mode: auto

biomarker:
  enabled: true
  case_group: case
  control_group: control
  outcome_column: DISEASE
  confounder_method: classic
```
