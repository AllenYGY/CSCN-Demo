# 安装说明

这份文档用于说明如何在本地安装和启动 `CSCN` 项目。

## 环境要求

- `Python >= 3.10`
- `Node.js >= 18`
- `npm >= 9`

如果只运行 CSCN 后端流程，`Node.js` 不是必须的；如果要使用网页 DAG viewer，则需要额外安装前端依赖并构建静态页面。

## 获取代码

```bash
git clone git@github.com:AllenYGY/CSCN-Demo.git
cd CSCN-Demo
```

如果你是通过其他 URL 克隆，进入你本地实际的仓库目录即可。

## 推荐安装方式

推荐使用项目内虚拟环境，避免系统 Python 的包管理限制。

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

安装基础依赖：

```bash
pip install -e .
```

如果要运行 viewer 和测试：

```bash
pip install -e .[viewer,dev]
```

## 前端安装

只有在需要网页 DAG viewer 时才需要这一步。

```bash
cd apps/dag-viewer
npm install
npm run build
cd ../..
```

构建完成后，viewer 静态文件会出现在 `apps/dag-viewer/dist/`。

## 安装后验证

先确认 CLI 已可用：

```bash
cscn --help
```

如果你的 shell 里还没有 `cscn` 命令，也可以直接这样验证：

```bash
python -m cscn.cli --help
```

如果要验证 viewer 依赖是否完整：

```bash
cscn viewer --config configs/GSE121893/config.yaml
```

然后在浏览器打开：

- [http://127.0.0.1:8000](http://127.0.0.1:8000)

## 数据目录约定

当前示例配置默认假设数据位于 `data/` 目录下，例如：

```text
data/
  GSE121893/
  GSE138852/
```

仓库不会跟踪这些大数据文件；如果你换了数据位置，需要同步修改 `configs/` 或你自己的配置文件。

## 常用安装组合

只跑 CSCN：

```bash
pip install -e .
```

跑 CSCN + viewer：

```bash
pip install -e .[viewer]
cd apps/dag-viewer
npm install
npm run build
cd ../..
```

跑 CSCN + viewer + 测试：

```bash
pip install -e .[viewer,dev]
cd apps/dag-viewer
npm install
npm run build
cd ../..
```

## 常见问题

### 1. `externally-managed-environment`

这是系统 Python 拒绝直接安装包的提示。不要继续往系统环境里装，改用项目虚拟环境：

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .[viewer,dev]
```

### 2. `cscn: command not found`

通常是因为：

- 当前没有激活虚拟环境
- 你用的不是安装过本项目的 Python

先执行：

```bash
source .venv/bin/activate
```

如果仍然有问题，用：

```bash
python -m cscn.cli --help
```

### 3. viewer 启动后提示前端构建产物不存在

说明还没执行前端构建：

```bash
cd apps/dag-viewer
npm install
npm run build
cd ../..
```

### 4. 读取 DAG pickle 时报 `pgmpy` 相关错误

原始 CSCN DAG pickle 依赖 `pgmpy`。请确认你安装的是带 viewer 的完整依赖：

```bash
pip install -e .[viewer]
```

如果你只装了基础依赖，某些 pickle 读取功能可能无法使用。

### 5. 示例配置可以解析，但运行时报找不到数据文件

这通常不是安装问题，而是本地缺少示例数据。请检查：

- `data/GSE121893/`
- `data/GSE138852/`

或者直接改用你自己的配置文件，把 `input.path`、`expr_path`、`metadata_path` 指向真实数据位置。

## 下一步

安装完成后，建议按这个顺序继续：

1. 阅读 [Configuration](configuration.md)
2. 查看 [Configs](../configs/README.md)
3. 如需网页可视化，阅读 [DAG Viewer](dag_viewer.md)
