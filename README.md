# CSCN

`CSCN` 现在以“通用 scRNA CSCN 仓库”为主，而不是只围绕特定 biomarker 数据集脚本组织。

当前命名约定：

- GitHub 仓库名：`CSCN-Demo`
- Python 包名：`cscn`
- CLI 名：`cscn`

当前正式入口是：

- Python 包：`src/cscn/`
- 统一 CLI：`cscn`
- 配置文件：YAML
- 标准输出目录：`runs/<run_name>/`

`biomarker/causal` 保留为可选后处理，DAG viewer 保留为官方可选能力。

## Install

完整安装步骤见：

- [Installation](docs/installation.md)

最短路径是：

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .[viewer,dev]
cd apps/dag-viewer
npm install
npm run build
cd ../..
```

## Quickstart

准备并运行一个配置：

```bash
cscn run-all --config examples/GSE121893/config.yaml
```

分步运行：

```bash
cscn prepare --config examples/GSE121893/config.yaml
cscn run --config examples/GSE121893/config.yaml
cscn aggregate --config examples/GSE121893/config.yaml
cscn biomarker --config examples/GSE121893/config.yaml
```

启动 viewer：

```bash
cscn viewer --config examples/GSE121893/config.yaml
```

## Supported Inputs

第一版支持两类输入：

- `.h5ad`
- 表达矩阵文件 + cell metadata 文件

最低输入契约是：

- 表达矩阵能解析成 `cells x genes`
- metadata 能提供 cell id
- 如需按 group 跑 CSCN，metadata 里要有一个 group 列

默认基因筛选是内置方差 Top-N；也可以通过 `gene_list_path` 提供外部 gene list 覆盖。

## Workflow

标准工作流采用 `staged + run-all`：

- `cscn prepare`
- `cscn run`
- `cscn aggregate`
- `cscn biomarker`
- `cscn viewer`
- `cscn run-all`

标准 run 目录结构：

```text
runs/<run_name>/
  config.snapshot.yaml
  inputs/
  matrices/
  dags/
  objects/
  consensus/
  biomarker/
  logs/
```

## Repo Layout

- `src/cscn/`: 正式主实现
- `src/biomarker/`: 现有 biomarker 与旧 viewer 兼容实现
- `apps/dag-viewer/`: DAG viewer 前端
- `examples/`: 新 CLI/config 的示例配置
- `docs/`: 使用与设计文档
- `scripts/biomarker/`: 旧的、数据集定制化脚本，保留作参考
- `legacy/`: 历史实验与旧实现

## Legacy Notes

仓库里仍然保留 `scripts/biomarker/Biomarker_GSE121893.py`、`scripts/biomarker/Biomarker_GSE138852.py` 这类旧脚本，但它们不再是正式主入口。根目录下同名脚本现在只是兼容 shim，会转发到 `cscn run-all --config examples/.../config.yaml`。

新的开发和新数据集接入，应该优先使用：

- `src/cscn/`
- `cscn` CLI
- `examples/` 里的配置模板

## Docs

- [Installation](docs/installation.md)
- [Configuration](docs/configuration.md)
- [Design And Development](docs/design_and_development.md)
- [DAG Viewer](docs/dag_viewer.md)
- [Examples](examples/README.md)
