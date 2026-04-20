# DAG Viewer

`cscn viewer` 是新的官方可选入口，用来浏览标准 `runs/<run_name>/` 目录中的单细胞 DAG 和 group 共识 DAG。

和旧版 viewer 相比，当前实现的核心变化不是“画图方式”而是“数据约定”：

- 输入不再要求是某个固定 biomarker 数据集目录
- viewer 默认面向 `cscn prepare/run/aggregate` 产出的标准 run 目录
- 后端优先读取预计算好的共识边 CSV，缺失时才动态回退到 pickle 聚合

如果你要看整体架构、目录设计和扩展点，请同时阅读 [Design And Development](design_and_development.md)。

## Scope

前端工程仍然位于：

- `apps/dag-viewer`

Python 后端入口位于：

- `src/cscn/viewer.py`

其中 `src/cscn/viewer.py` 负责：

- 暴露新的 FastAPI API
- 面向标准 run 目录做 root/group/graph 发现
- 托管前端静态资源

而图对象归一化、pickle 兼容读取、共识 CSV 反序列化等基础逻辑，仍复用：

- `src/biomarker/dag_viewer/service.py`

这意味着当前 viewer 架构是“新 run-root 适配层 + 旧图序列化服务”的组合，而不是完全重写。

## Compatible Run Layout

viewer 认为一个目录是兼容 run root，当且仅当它至少满足：

- 目录存在
- `dags/` 子目录存在
- `inputs/genes.csv` 文件存在

标准目录结构如下：

```text
runs/<run_name>/
  config.snapshot.yaml
  inputs/
    genes.csv
    groups.csv
    cell_metadata.csv
    run_summary.json
  matrices/
    <group>.npy
    <group>_cells.csv
  dags/
    <group>/
      result_0.pkl
      result_1.pkl
      ...
  objects/
    <group>_cscn.pkl
  consensus/
    <group>.csv
    summary.json
  biomarker/
  logs/
```

viewer 重点依赖的文件：

- `inputs/genes.csv`
- `inputs/groups.csv`
- `matrices/<group>_cells.csv`
- `dags/<group>/result_*.pkl`
- `consensus/<group>.csv`

## Data Discovery Model

### Root Discovery

`GET /api/roots` 会扫描 `run.output_dir`，只收集兼容 run root 的目录。

每个 root 会返回：

- `path`
- `label`
- `compatible`
- `rootKind`
- `collectionKind`
- `collectionLabel`
- `collectionCount`

当前实现中：

- `rootKind` 固定为 `run`
- `collectionKind` 固定为 `single`
- `collectionLabel` 固定为 `当前视图`

也就是说，新版 viewer 的一个 root 就对应一个 `cscn` 运行结果目录，而不是旧模型中的“dataset root + 多个方法集合”。

### Group Discovery

`POST /api/scan` 会扫描单个 run root 下的所有 group：

- 枚举 `dags/*`
- 仅保留内部存在 `result_*.pkl` 的 group 目录
- 从 `inputs/groups.csv` 读取更友好的 group label

每个 group 返回：

- `method`
- `groupKey`
- `time`
- `cellType`
- `numDags`
- `hasNodeMap`
- `hasConsensusCsv`
- `resultIds`

当前 run 模式下这些字段的语义是：

- `method`：固定取 run 名的安全版本，用于兼容前端已有状态模型
- `groupKey`：真实 group 目录名
- `time`：当前总是空字符串
- `cellType`：优先取 `groups.csv` 里的 `group_label`
- `numDags`：该 group 下 `result_*.pkl` 数量
- `hasNodeMap`：`inputs/genes.csv` 是否存在
- `hasConsensusCsv`：`consensus/<group>.csv` 是否存在
- `resultIds`：可用的单细胞 DAG 索引

## Graph Loading Semantics

`POST /api/graph` 支持两种模式：

- `single`
- `consensus`

### Single Mode

请求参数：

- `rootPath`
- `method`
- `groupKey`
- `mode: "single"`
- `resultId`

后端行为：

1. 校验 `rootPath` 是否是兼容 run root
2. 校验 `method` 是否等于 run 名
3. 定位 `dags/<group>/result_<resultId>.pkl`
4. 反序列化 pickle 图对象
5. 用 `inputs/genes.csv` 中的节点映射把整数节点转成基因名
6. 根据 `matrices/<group>_cells.csv` 把 `resultId` 反查到对应 cell id
7. 输出标准化后的节点、边和元数据

重要约定：

- `result_<i>.pkl` 的 `i` 被视作矩阵行索引，而不是任意主键
- `matrices/<group>_cells.csv` 的第 `i` 行就是该 DAG 对应的 cell id

### Consensus Mode

请求参数：

- `rootPath`
- `method`
- `groupKey`
- `mode: "consensus"`

后端行为分两种路径：

#### 路径 1：读取预计算 CSV

如果 `consensus/<group>.csv` 存在：

- 直接从 CSV 读取 `from/to/count/threshold/num_dags`
- 仅构造被共识边引用到的节点

#### 路径 2：动态聚合 pickle

如果 `consensus/<group>.csv` 不存在：

- 枚举该 group 下所有 `result_*.pkl`
- 逐个归一化图对象
- 统计边出现次数
- 使用 `auto` 阈值生成共识边

默认动态阈值：

```text
max(2, ceil(0.05 * num_dags))
```

因此共识图通常具有两个特点：

- 边数明显少于所有单细胞 DAG 边的并集
- 节点数只覆盖参与共识边的节点，不包含纯孤立节点

## API

### `GET /api/health`

返回依赖状态和前端构建状态。

示例响应：

```json
{
  "dependencies": {
    "fastapi": true,
    "uvicorn": true,
    "pgmpy": true,
    "networkx": true,
    "pandas": true
  },
  "supportedRootKinds": ["dataset", "dag"],
  "frontendBuilt": true,
  "frontendDist": "/abs/path/apps/dag-viewer/dist"
}
```

说明：

- `supportedRootKinds` 来自旧 viewer service 的兼容状态，当前 `cscn viewer` 实际暴露的是 run-root 模型
- `frontendBuilt` 用于快速判断前端静态资源是否已经构建

### `GET /api/roots`

扫描 `run.output_dir` 下的 run root。

示例：

```json
{
  "roots": [
    {
      "path": "/abs/path/runs/gse121893_generic",
      "label": "gse121893_generic",
      "compatible": true,
      "rootKind": "run",
      "collectionKind": "single",
      "collectionLabel": "当前视图",
      "collectionCount": 1
    }
  ],
  "defaultRootPath": "/abs/path/runs/gse121893_generic",
  "dataRoot": "/abs/path/runs"
}
```

### `POST /api/scan`

请求体：

```json
{
  "rootPath": "/abs/path/runs/gse121893_generic"
}
```

响应会描述这个 run 中有哪些 group 可以被浏览。

### `POST /api/graph`

单细胞 DAG 请求体：

```json
{
  "rootPath": "/abs/path/runs/gse121893_generic",
  "method": "gse121893_generic",
  "groupKey": "dHF",
  "mode": "single",
  "resultId": 0
}
```

共识 DAG 请求体：

```json
{
  "rootPath": "/abs/path/runs/gse121893_generic",
  "method": "gse121893_generic",
  "groupKey": "dHF",
  "mode": "consensus"
}
```

### Graph Response Shape

无论 `single` 还是 `consensus`，后端都会返回：

```json
{
  "mode": "single",
  "nodes": [
    {
      "id": "0",
      "label": "TNNT2",
      "rawId": "0",
      "mapped": true,
      "degree": 3,
      "inDegree": 1,
      "outDegree": 2
    }
  ],
  "edges": [
    {
      "id": "0->1",
      "source": "0",
      "target": "1",
      "count": 1,
      "weight": 1
    }
  ],
  "metadata": {
    "method": "gse121893_generic",
    "groupKey": "dHF",
    "numNodes": 42,
    "numEdges": 51,
    "resultId": 0,
    "cellRunId": "AAAC..."
  }
}
```

元数据差异：

- `single` 常见字段：`resultId`、`cellRunId`
- `consensus` 常见字段：`threshold`、`numDags`

## Frontend Serving

`cscn viewer` 会同时承担 API 服务和前端静态文件托管：

- `/assets/*` 挂载到 `apps/dag-viewer/dist/assets`
- `/` 返回 `dist/index.html`
- 非 API 路由统一走 SPA fallback

如果前端未构建，根路径会返回一个内嵌的“缺少前端构建产物”页面，而不是直接 404。

## Installation

Python 依赖：

```bash
pip install -e .[viewer]
```

前端构建：

```bash
cd apps/dag-viewer
npm install
npm run build
```

## Run

```bash
cscn viewer --config examples/GSE121893/config.yaml
```

默认访问地址：

- [http://127.0.0.1:8000](http://127.0.0.1:8000)

注意：

- `viewer.enabled` 目前不会阻止命令启动
- viewer 会把 `run.output_dir` 作为可扫描根目录
- 当前配置对应的 `run_name` 会作为默认 root 优先打开

## Error Model

后端对业务错误统一抛出 `DagViewerError`，再在 FastAPI 层转换成 `400 Bad Request`。

常见错误包括：

- `rootPath` 为空
- 目录不是兼容 run root
- `method` 与当前 run 名不匹配
- `groupKey` 不存在
- `single` 模式缺少 `resultId`
- `resultId` 不是整数
- 指定的 `result_<id>.pkl` 文件不存在
- pickle 依赖的类或第三方库缺失

## Development Notes

### 为什么还依赖 `src/biomarker/dag_viewer/service.py`

原因很实际：

- 旧 service 已经包含了 pickle 图对象兼容读取逻辑
- 旧 service 已经把多种图对象归一化为统一的 nodes/edges 结构
- 新 viewer 只需要重做“如何发现数据根目录”和“如何对接 run 布局”

所以当前分层更像：

```text
src/cscn/viewer.py
  -> 负责 run root 发现、API 暴露、前端托管
  -> 调用 biomarker.dag_viewer.service 中的图归一化能力
```

### 如果要扩展 viewer

常见改动点：

- 新增 root 扫描规则：修改 `_is_run_root` 与 `discover_run_roots`
- 新增 group 元数据：修改 `_load_group_labels`、`scan_run_root`
- 新增图模式：修改 `load_run_graph` 和前端状态模型
- 替换共识逻辑：优先修改 `src/cscn/aggregate.py`，viewer 只负责消费产物

### 建议的演进方向

- 把 run viewer 与 legacy viewer 的共用序列化逻辑抽到 `src/cscn/viewer_service.py`
- 让 `GET /api/health` 反映真实的 run-root 能力，而不是旧模型里的 root kind
- 为 graph payload 增加版本号，降低前后端协同改动成本
