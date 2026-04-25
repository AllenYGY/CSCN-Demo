# Configs

这些配置是新 `cscn` CLI 的配置模板，不再依赖旧的 dataset-specific 主脚本。

## Available Configs

- `GSE121893/config.yaml`
- `GSE138852/config.yaml`
- `SCP2046/sham1_full_adaptive_block_prior.yaml`

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
