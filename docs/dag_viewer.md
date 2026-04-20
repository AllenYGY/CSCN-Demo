# DAG Viewer

`cscn viewer` 是新的官方可选入口，用来浏览标准 run 目录里的单细胞 DAG 和共识 DAG。

前端工程仍然在：

- `apps/dag-viewer`

## Viewer Scope

viewer 现在默认面向 `runs/<run_name>/` 结构，而不是只面向 biomarker 数据集目录。

它会自动发现：

- `runs/<run_name>/dags/<group>/result_*.pkl`
- `runs/<run_name>/consensus/<group>.csv`
- `runs/<run_name>/inputs/genes.csv`
- `runs/<run_name>/matrices/<group>_cells.csv`

## 界面行为

- 根目录下拉列表来自 `run.output_dir`
- 每个 run 默认显示成单一“当前视图 / 默认视图”
- `Group` 列表只有同时存在时间和 cell type 时才显示 `Day`
- 默认勾选“仅保留有边相连的节点”
- 可以取消勾选查看孤立节点

## 共识图说明

- 如果 `consensus/<group>.csv` 已存在，viewer 直接读取它
- 否则会从 `dags/<group>/result_*.pkl` 动态聚合
- 动态聚合阈值默认是：

```text
max(2, ceil(0.05 * num_dags))
```

- 共识图只保留被共识边连接到的节点，所以节点数通常少于单细胞 DAG

## 安装

Python：

```bash
pip install -e .[viewer]
```

前端：

```bash
cd apps/dag-viewer
npm install
npm run build
```

## 运行

```bash
cscn viewer --config examples/GSE121893/config.yaml
```

浏览器打开：

- [http://127.0.0.1:8000](http://127.0.0.1:8000)

## API

- `GET /api/health`
- `GET /api/roots`
- `POST /api/scan`
- `POST /api/graph`

`/api/scan` 请求体示例：

```json
{
  "rootPath": "runs/gse121893_generic"
}
```

`/api/graph` 请求体示例：

```json
{
  "rootPath": "runs/gse121893_generic",
  "method": "gse121893_generic",
  "groupKey": "dHF",
  "mode": "single",
  "resultId": 0
}
```
