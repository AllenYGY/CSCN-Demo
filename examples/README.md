# Examples

这些配置是新 `cscn` CLI 的示例入口，不再依赖旧的 dataset-specific 主脚本。

## Available Examples

- `GSE121893/config.yaml`
- `GSE138852/config.yaml`

## Run

```bash
cscn run-all --config examples/GSE121893/config.yaml
```

分步运行：

```bash
cscn prepare --config examples/GSE121893/config.yaml
cscn run --config examples/GSE121893/config.yaml
cscn aggregate --config examples/GSE121893/config.yaml
```

启动 viewer：

```bash
cscn viewer --config examples/GSE121893/config.yaml
```

这些 examples 的目标是演示：

- 如何把真实 scRNA 数据接到统一 config
- 如何从标准 run 目录继续做 aggregate / biomarker / viewer
- 如何逐步替换旧的 `scripts/biomarker/*.py`
