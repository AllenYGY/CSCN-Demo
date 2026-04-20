# CSCN Configuration

`cscn` 使用 YAML 配置文件。推荐先从 `examples/` 里的模板复制，再按数据集修改。

## Minimal h5ad Example

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

## Minimal tables Example

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

## Key Fields

### `input`

- `format`: `h5ad` or `tables`
- `path`: `.h5ad` 文件路径
- `expr_path`: 表达矩阵文件路径
- `metadata_path`: metadata 文件路径
- `layer`: `.h5ad` 时可选，默认读 `X`
- `gene_key`: `.h5ad` 时可选，默认读 `var_names`
- `obs_group_key`: 用于按 group 跑 CSCN 的 metadata 列
- `expr_orientation`: `cells_by_genes` or `genes_by_cells`
- `expr_cell_id_column`: 表达矩阵中的 cell id 列名
- `metadata_cell_id_column`: metadata 中的 cell id 列名

### `preprocess`

- `normalize`: 是否做 library-size normalization
- `log1p`: 是否做 `log1p`
- `sample_per_group`: 每个 group 采样多少细胞；省略则保留全部
- `random_seed`: 随机种子
- `gene_selection.top_n`: 默认内置筛选保留多少基因
- `gene_selection.gene_list_path`: 外部 gene list，存在时优先生效

### `run`

- `output_dir`: 标准 run 输出根目录
- `groups`: 可选，只跑指定 groups
- `max_workers`: 并行跑 PC 的 worker 数
- `sigmoid_score`
- `significance_level`
- `max_cond_vars`
- `use_bitmap`
- `using_nmf`

### `aggregate`

- `consensus`: 是否生成共识边
- `consensus_threshold_mode`: `auto` 或整数

`auto` 当前定义为：

```text
max(2, ceil(0.05 * num_dags))
```

### `biomarker`

- `enabled`: 是否在 `run-all` 时执行 biomarker
- `case_group`
- `control_group`
- `outcome_column`
- `confounder_method`

`biomarker` 当前只支持明确的二组比较。

### `viewer`

- `enabled`
- `host`
- `port`
- `reload`

`cscn viewer` 会从 `run.output_dir` 自动发现可浏览的 run。
