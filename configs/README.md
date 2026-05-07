# Configs

这些配置是新 `cscn` CLI 的配置模板，不再依赖旧的 dataset-specific 主脚本。

## Available Configs

- `GSE121893/config.yaml`
- `GSE138852/config.yaml`
- `GSE128639/rna_only.yaml`
- `GSE128639/adt_only.yaml`
- `GSE128639/rna_adt_joint.yaml`
- `GSE128639/shared3000_rna_only.yaml`
- `GSE128639/shared3000_adt_only.yaml`
- `GSE128639/shared3000_rna_adt_joint.yaml`
- `GSE164378/rna_only.yaml`
- `GSE164378/adt_only.yaml`
- `GSE164378/rna_adt_joint.yaml`
- `SCP2046/sham1_full_adaptive_block_prior.yaml`

## GSE128639 Preparation

`GSE128639` 需要先把 GEO `MNC` 原始矩阵整理成 CSCN 可直接读取的精选输入表：

```bash
python3 scripts/prep/prepare_GSE128639.py
```

默认会在 `data/GSE128639/cscn_inputs/` 生成：

- `gse128639_mnc_metadata.csv.gz`
- `gse128639_mnc_rna_only_expression.tsv.gz`
- `gse128639_mnc_adt_only_expression.tsv.gz`
- `gse128639_mnc_rna_adt_joint_expression.tsv.gz`

同时会额外生成一套固定 `shared3000` 子集：

- 总细胞数：`3000`
- 分层字段：`hto_best_label`
- 随机种子：`42`
- 代表文件：
  - `gse128639_mnc_shared3000_metadata.csv.gz`
  - `gse128639_mnc_shared3000_rna_only_expression.tsv.gz`
  - `gse128639_mnc_shared3000_adt_only_expression.tsv.gz`
  - `gse128639_mnc_shared3000_rna_adt_joint_expression.tsv.gz`

当前 `GSE128639` 输出使用统一标准化后的 cell id（将 RNA/ADT 中的 `.1` 后缀转换为 `-1`），并附带基于 HTO 原始计数的最小 metadata 摘要字段。

## GSE164378 Preparation

`GSE164378` 需要先把 GEO `3P` 原始矩阵整理成 CSCN 可直接读取的 `cells x features` 表格：

```bash
python3 scripts/prep/prepare_GSE164378.py
```

默认会在 `data/GSE164378/cscn_inputs/` 生成：

- `gse164378_3p_rna_only_expression.csv.gz`
- `gse164378_3p_adt_only_expression.csv.gz`
- `gse164378_3p_rna_adt_joint_expression.csv.gz`
- `gse164378_3p_metadata.csv.gz`

当前 `GSE164378` 预处理会额外生成一套固定 shared-cell 子集：

- 总细胞数：`2000`
- 分层字段：`celltype.l1`
- 当前默认分配：8 个 `celltype.l1` 大类各 `250`

当前三份 `GSE164378` config 默认直接使用这套 fixed shared-cell 输入：

- 单个 CSCN 运行组：`all`
- 不再在 `prepare` 阶段二次随机采样

## Run

```bash
cscn run-all --config configs/GSE121893/config.yaml
```

分步运行：

```bash
cscn prepare --config configs/GSE121893/config.yaml
cscn run --config configs/GSE121893/config.yaml
cscn aggregate --config configs/GSE121893/config.yaml
```

启动 viewer：

```bash
cscn viewer --config configs/GSE121893/config.yaml
```

这些 examples 的目标是演示：

- 如何把真实 scRNA 数据接到统一 config
- 如何从标准 run 目录继续做 aggregate / biomarker / viewer
- 如何逐步替换旧的 `scripts/biomarker/*.py`
