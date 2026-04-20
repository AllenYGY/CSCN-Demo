import { startTransition, useEffect, useRef, useState } from "react";
import GraphCanvas from "./components/GraphCanvas";
import { fetchGraph, fetchHealth, fetchRoots, scanDagRoot } from "./api";
import type {
  CollectionOption,
  GraphMode,
  GraphResponse,
  HealthResponse,
  RootOption,
  RootsResponse,
  ScanGroup,
  ScanResponse,
} from "./types";

export default function App() {
  const scanRequestRef = useRef(0);
  const graphRequestRef = useRef(0);
  const [rootsData, setRootsData] = useState<RootsResponse | null>(null);
  const [selectedRootPath, setSelectedRootPath] = useState("");
  const [health, setHealth] = useState<HealthResponse | null>(null);
  const [scanData, setScanData] = useState<ScanResponse | null>(null);
  const [selectedCollectionKey, setSelectedCollectionKey] = useState("");
  const [selectedGroupKey, setSelectedGroupKey] = useState("");
  const [mode, setMode] = useState<GraphMode>("single");
  const [selectedResultId, setSelectedResultId] = useState<number | undefined>(undefined);
  const [graph, setGraph] = useState<GraphResponse | null>(null);
  const [searchText, setSearchText] = useState("");
  const [showAllLabels, setShowAllLabels] = useState(false);
  const [onlyConnectedNodes, setOnlyConnectedNodes] = useState(true);
  const [selectedNodeId, setSelectedNodeId] = useState<string | null>(null);
  const [selectedEdgeId, setSelectedEdgeId] = useState<string | null>(null);
  const [scanError, setScanError] = useState("");
  const [graphError, setGraphError] = useState("");
  const [rootsError, setRootsError] = useState("");
  const [isScanning, setIsScanning] = useState(false);
  const [isLoadingGraph, setIsLoadingGraph] = useState(false);

  useEffect(() => {
    void fetchHealth().then(setHealth).catch(() => null);
    void fetchRoots()
      .then((nextRoots) => {
        setRootsData(nextRoots);
        setRootsError("");
        const compatibleRoots = nextRoots.roots.filter((root) => root.compatible);
        const fallbackRoot = nextRoots.defaultRootPath ?? compatibleRoots[0]?.path ?? "";
        setSelectedRootPath(fallbackRoot);
      })
      .catch((error) => {
        setRootsError(error instanceof Error ? error.message : "加载数据集列表失败。");
      });
  }, []);

  useEffect(() => {
    if (!selectedRootPath) {
      return;
    }
    void handleScan(selectedRootPath);
  }, [selectedRootPath]);

  async function handleScan(rootPath: string = selectedRootPath) {
    const requestId = ++scanRequestRef.current;
    setIsScanning(true);
    setScanError("");
    setGraphError("");
    setScanData(null);
    setGraph(null);
    startTransition(() => {
      setSelectedCollectionKey("");
      setSelectedGroupKey("");
      setSelectedResultId(undefined);
      setMode("single");
      setSelectedNodeId(null);
      setSelectedEdgeId(null);
      setSearchText("");
      setShowAllLabels(false);
    });
    try {
      const nextScan = await scanDagRoot(rootPath);
      if (requestId !== scanRequestRef.current) {
        return;
      }
      setScanData(nextScan);

      const firstCollection = nextScan.collections[0];
      const firstCollectionKey = firstCollection?.key ?? "";
      const firstGroups = nextScan.groups.filter((group) => group.method === firstCollectionKey);
      const firstGroup = firstGroups[0];

      startTransition(() => {
        setSelectedCollectionKey(firstCollectionKey);
        setSelectedGroupKey(firstGroup?.groupKey ?? "");
        setSelectedResultId(firstGroup?.resultIds[0]);
        setMode("single");
        setGraph(null);
        setSelectedNodeId(null);
        setSelectedEdgeId(null);
        setSearchText("");
        setShowAllLabels(false);
      });
    } catch (error) {
      if (requestId !== scanRequestRef.current) {
        return;
      }
      setScanData(null);
      setGraph(null);
      setScanError(error instanceof Error ? error.message : "扫描目录失败。");
    } finally {
      if (requestId === scanRequestRef.current) {
        setIsScanning(false);
      }
    }
  }

  const compatibleRoots = rootsData?.roots.filter((root) => root.compatible) ?? [];
  const selectedRoot = rootsData?.roots.find((root) => root.path === selectedRootPath) ?? null;
  const collections = scanData?.collections ?? [];
  const displayGraph = filterGraphByConnectivity(graph, onlyConnectedNodes);
  const graphInsights = graph ? summarizeGraph(graph, displayGraph ?? graph) : null;
  const selectedCollection =
    collections.find((collection) => collection.key === selectedCollectionKey) ?? null;
  const groupsForCollection = scanData
    ? scanData.groups.filter((group) => group.method === selectedCollectionKey)
    : [];
  const selectedGroup =
    groupsForCollection.find((group) => group.groupKey === selectedGroupKey) ?? null;
  const selectedNode = displayGraph?.nodes.find((node) => node.id === selectedNodeId) ?? null;
  const selectedEdge = displayGraph?.edges.find((edge) => edge.id === selectedEdgeId) ?? null;
  const selectedNodeUpstream = selectedNode
    ? displayGraph?.edges
        .filter((edge) => edge.target === selectedNode.id)
        .map(
          (edge) =>
            displayGraph.nodes.find((node) => node.id === edge.source)?.label ?? edge.source,
        )
        .sort((left, right) => left.localeCompare(right))
    : [];
  const selectedNodeDownstream = selectedNode
    ? displayGraph?.edges
        .filter((edge) => edge.source === selectedNode.id)
        .map(
          (edge) =>
            displayGraph.nodes.find((node) => node.id === edge.target)?.label ?? edge.target,
        )
        .sort((left, right) => left.localeCompare(right))
    : [];

  useEffect(() => {
    if (!selectedGroup && groupsForCollection.length > 0) {
      const fallbackGroup = groupsForCollection[0];
      setSelectedGroupKey(fallbackGroup.groupKey);
      setSelectedResultId(fallbackGroup.resultIds[0]);
    }
  }, [groupsForCollection, selectedGroup]);

  useEffect(() => {
    if (!selectedGroup) {
      return;
    }
    if (mode === "single") {
      if (selectedResultId === undefined || !selectedGroup.resultIds.includes(selectedResultId)) {
        setSelectedResultId(selectedGroup.resultIds[0]);
      }
    }
  }, [mode, selectedGroup, selectedResultId]);

  useEffect(() => {
    if (!scanData || scanData.rootPath !== selectedRootPath || !selectedCollectionKey || !selectedGroupKey) {
      return;
    }
    if (mode === "single" && selectedResultId === undefined) {
      return;
    }

    const requestId = ++graphRequestRef.current;
    setIsLoadingGraph(true);
    setGraphError("");
    void fetchGraph({
      rootPath: selectedRootPath,
      method: selectedCollectionKey,
      groupKey: selectedGroupKey,
      mode,
      resultId: mode === "single" ? selectedResultId : undefined,
    })
      .then((nextGraph) => {
        if (requestId !== graphRequestRef.current) {
          return;
        }
        setGraph(nextGraph);
        setSelectedNodeId(null);
        setSelectedEdgeId(null);
      })
      .catch((error) => {
        if (requestId !== graphRequestRef.current) {
          return;
        }
        setGraph(null);
        setGraphError(error instanceof Error ? error.message : "加载图失败。");
      })
      .finally(() => {
        if (requestId === graphRequestRef.current) {
          setIsLoadingGraph(false);
        }
      });
  }, [selectedRootPath, scanData, selectedCollectionKey, selectedGroupKey, mode, selectedResultId]);

  useEffect(() => {
    if (!displayGraph) {
      return;
    }
    if (selectedNodeId && !displayGraph.nodes.some((node) => node.id === selectedNodeId)) {
      setSelectedNodeId(null);
    }
    if (selectedEdgeId && !displayGraph.edges.some((edge) => edge.id === selectedEdgeId)) {
      setSelectedEdgeId(null);
    }
  }, [displayGraph, selectedEdgeId, selectedNodeId]);

  function selectCollection(nextCollectionKey: string) {
    const nextGroups = scanData?.groups.filter((group) => group.method === nextCollectionKey) ?? [];
    const firstGroup = nextGroups[0];
    startTransition(() => {
      setSelectedCollectionKey(nextCollectionKey);
      setSelectedGroupKey(firstGroup?.groupKey ?? "");
      setSelectedResultId(firstGroup?.resultIds[0]);
      setSelectedNodeId(null);
      setSelectedEdgeId(null);
    });
  }

  function selectGroup(nextGroupKey: string) {
    const nextGroup = groupsForCollection.find((group) => group.groupKey === nextGroupKey);
    setSelectedGroupKey(nextGroupKey);
    setSelectedResultId(nextGroup?.resultIds[0]);
    setSelectedNodeId(null);
    setSelectedEdgeId(null);
  }

  function selectMode(nextMode: GraphMode) {
    setMode(nextMode);
    setSelectedNodeId(null);
    setSelectedEdgeId(null);
  }

  return (
    <div className="app-shell">
      <aside className="workbench-panel controls-panel">
        <div className="panel-section">
          <p className="eyebrow">CSCN DAG Viewer</p>
          <h1>CSCN 可视化</h1>
          <p className="lede">
            按数据集自动识别可用 DAG 视图。不同目录结构会用更合适的层级呈现，而不是一律套成筛法。
          </p>
        </div>

        <div className="panel-section">
          <label className="field-label" htmlFor="rootSelect">
            数据集
          </label>
          <select
            id="rootSelect"
            className="select-input"
            value={selectedRootPath}
            onChange={(event) => setSelectedRootPath(event.target.value)}
            disabled={!compatibleRoots.length}
          >
            {!compatibleRoots.length ? <option value="">没有可用数据集</option> : null}
            {(rootsData?.roots ?? []).map((root) => (
              <option
                key={root.path}
                value={root.path}
                disabled={!root.compatible}
              >
                {root.compatible
                  ? formatRootLabel(root)
                  : `${root.label} · 不可用`}
              </option>
            ))}
          </select>
          <button
            className="primary-button"
            onClick={() => void handleScan()}
            disabled={isScanning || !selectedRootPath}
          >
            {isScanning ? "扫描中…" : "扫描目录"}
          </button>
          {selectedRoot ? (
            <p className="lede small">
              {selectedRoot.path}
            </p>
          ) : null}
          {rootsError ? <p className="error-banner">{rootsError}</p> : null}
          {scanError ? <p className="error-banner">{scanError}</p> : null}
        </div>

        {scanData ? (
          <div className="panel-section">
            <span className="field-label">{scanData.collectionLabel}</span>
            {scanData.showCollectionSelector ? (
              <div className="chip-row">
                {collections.map((collection) => (
                  <button
                    key={collection.key}
                    className={collection.key === selectedCollectionKey ? "chip active" : "chip"}
                    onClick={() => selectCollection(collection.key)}
                  >
                    {collection.displayLabel}
                  </button>
                ))}
              </div>
            ) : selectedCollection ? (
              <div className="detail-card">
                <h4>{selectedCollection.displayLabel}</h4>
                <p>{selectedCollection.rawLabel}</p>
              </div>
            ) : null}
          </div>
        ) : null}

        <div className="panel-section">
          <label className="field-label" htmlFor="groupSelect">
            Group
          </label>
          <select
            id="groupSelect"
            className="select-input"
            value={selectedGroupKey}
            onChange={(event) => selectGroup(event.target.value)}
          >
            {groupsForCollection.map((group) => (
              <option key={`${group.method}:${group.groupKey}`} value={group.groupKey}>
                {formatGroupLabel(group)}
              </option>
            ))}
          </select>
          {selectedGroup ? (
            <div className="meta-grid">
              <div>
                <span className="meta-label">cell type</span>
                <strong>{selectedGroup.cellType || "全部"}</strong>
              </div>
              <div>
                <span className="meta-label">time</span>
                <strong>{selectedGroup.time || "未知"}</strong>
              </div>
              <div>
                <span className="meta-label">DAG 数</span>
                <strong>{selectedGroup.numDags}</strong>
              </div>
              <div>
                <span className="meta-label">共识 CSV</span>
                <strong>{selectedGroup.hasConsensusCsv ? "有" : "回退聚合"}</strong>
              </div>
            </div>
          ) : null}
        </div>

        <div className="panel-section">
          <span className="field-label">模式</span>
          <div className="chip-row">
            <button
              className={mode === "single" ? "chip active" : "chip"}
              onClick={() => selectMode("single")}
            >
              单细胞
            </button>
            <button
              className={mode === "consensus" ? "chip active" : "chip"}
              onClick={() => selectMode("consensus")}
            >
              共识
            </button>
          </div>

          {mode === "single" && selectedGroup ? (
            <>
              <label className="field-label" htmlFor="resultSelect">
                result_id
              </label>
              <select
                id="resultSelect"
                className="select-input"
                value={selectedResultId ?? ""}
                onChange={(event) => setSelectedResultId(Number.parseInt(event.target.value, 10))}
              >
                {selectedGroup.resultIds.map((resultId) => (
                  <option key={resultId} value={resultId}>
                    result_{resultId}
                  </option>
                ))}
              </select>
            </>
          ) : null}
        </div>

        <div className="panel-section">
          <label className="field-label" htmlFor="searchInput">
            基因搜索
          </label>
          <input
            id="searchInput"
            className="select-input"
            value={searchText}
            onChange={(event) => setSearchText(event.target.value)}
            placeholder="输入 gene / raw id"
          />
          <label className="checkbox-row">
            <input
              type="checkbox"
              checked={onlyConnectedNodes}
              onChange={(event) => setOnlyConnectedNodes(event.target.checked)}
            />
            <span>仅保留有边相连的节点</span>
          </label>
          <label className="checkbox-row">
            <input
              type="checkbox"
              checked={showAllLabels}
              onChange={(event) => setShowAllLabels(event.target.checked)}
            />
            <span>强制显示全部标签</span>
          </label>
        </div>

        {health ? (
          <div className="panel-section">
            <span className="field-label">运行依赖</span>
            <div className="dependency-list">
              {Object.entries(health.dependencies).map(([name, installed]) => (
                <span key={name} className={installed ? "status-pill ok" : "status-pill warn"}>
                  {name}: {installed ? "OK" : "缺失"}
                </span>
              ))}
            </div>
          </div>
        ) : null}
      </aside>

      <main className="workspace-panel">
        <section className="workspace-header workbench-panel">
          <div>
            <p className="eyebrow">当前视图</p>
            <h2>
              {selectedGroup
                ? formatGroupHeading(selectedGroup, mode)
                : rootsError
                  ? "等待数据集列表"
                  : "等待目录扫描"}
            </h2>
            <p className="lede small">
              {displayGraph?.metadata.cellRunId
                ? `对应 cell/run: ${displayGraph.metadata.cellRunId}`
                : selectedCollection
                  ? `${scanData?.collectionLabel ?? "当前视图"}: ${selectedCollection.rawLabel}`
                  : "单细胞模式会在元数据可用时显示原始 run id。"}
            </p>
          </div>
          <div className="status-cluster">
            {displayGraph ? <span className="status-pill ok">{displayGraph.metadata.numNodes} nodes</span> : null}
            {displayGraph ? <span className="status-pill ok">{displayGraph.metadata.numEdges} edges</span> : null}
            {displayGraph?.metadata.numDags ? (
              <span className="status-pill ok">num_dags {displayGraph.metadata.numDags}</span>
            ) : null}
            {displayGraph?.metadata.threshold ? (
              <span className="status-pill ok">threshold {displayGraph.metadata.threshold}</span>
            ) : null}
            {onlyConnectedNodes ? <span className="status-pill warn">已隐藏孤立节点</span> : null}
            {isLoadingGraph ? <span className="status-pill warn">载入中</span> : null}
          </div>
        </section>

        <section className="graph-panel workbench-panel">
          {graphError ? <p className="error-banner">{graphError}</p> : null}
          {!displayGraph && !graphError ? (
            <div className="empty-state">
              <p>扫描目录后会在这里显示 DAG 画布。</p>
            </div>
          ) : (
            <GraphCanvas
              graph={displayGraph}
              searchText={searchText}
              showAllLabels={showAllLabels}
              selectedNodeId={selectedNodeId}
              selectedEdgeId={selectedEdgeId}
              onSelectNode={(nodeId) => {
                setSelectedNodeId(nodeId);
                setSelectedEdgeId(null);
              }}
              onSelectEdge={(edgeId) => {
                setSelectedEdgeId(edgeId);
                setSelectedNodeId(null);
              }}
              onClearSelection={() => {
                setSelectedNodeId(null);
                setSelectedEdgeId(null);
              }}
            />
          )}
        </section>
      </main>

      <aside className="workbench-panel detail-panel">
        <div className="panel-section">
          <p className="eyebrow">细节面板</p>
          <h3>图元详情</h3>
          {!selectedNode && !selectedEdge ? (
            <p className="lede small">点击节点查看上下游，点击边查看方向与权重。</p>
          ) : null}
        </div>

        {selectedNode ? (
          <div className="panel-section">
            <span className="field-label">节点</span>
            <div className="detail-card">
              <h4>{selectedNode.label}</h4>
              <p>raw id: {selectedNode.rawId}</p>
              <p>mapped: {selectedNode.mapped ? "是" : "否"}</p>
              <p>
                degree / in / out: {selectedNode.degree} / {selectedNode.inDegree} /{" "}
                {selectedNode.outDegree}
              </p>
            </div>

            <span className="field-label">上游</span>
            <ul className="detail-list">
              {selectedNodeUpstream?.length ? (
                selectedNodeUpstream.map((label) => <li key={`up:${label}`}>{label}</li>)
              ) : (
                <li>无</li>
              )}
            </ul>

            <span className="field-label">下游</span>
            <ul className="detail-list">
              {selectedNodeDownstream?.length ? (
                selectedNodeDownstream.map((label) => <li key={`down:${label}`}>{label}</li>)
              ) : (
                <li>无</li>
              )}
            </ul>
          </div>
        ) : null}

        {selectedEdge ? (
          <div className="panel-section">
            <span className="field-label">边</span>
            <div className="detail-card">
              <h4>{selectedEdge.id}</h4>
              <p>
                {displayGraph?.nodes.find((node) => node.id === selectedEdge.source)?.label ?? selectedEdge.source}
                {" → "}
                {displayGraph?.nodes.find((node) => node.id === selectedEdge.target)?.label ?? selectedEdge.target}
              </p>
              <p>count: {selectedEdge.count}</p>
              <p>weight: {selectedEdge.weight}</p>
            </div>
          </div>
        ) : null}

        {displayGraph ? (
          <div className="panel-section">
            <span className="field-label">图摘要</span>
            <div className="meta-grid">
              <div>
                <span className="meta-label">method</span>
                <strong>{selectedCollection?.displayLabel ?? displayGraph.metadata.method}</strong>
              </div>
              <div>
                <span className="meta-label">group</span>
                <strong>{displayGraph.metadata.groupKey}</strong>
              </div>
              <div>
                <span className="meta-label">mode</span>
                <strong>{displayGraph.mode}</strong>
              </div>
              <div>
                <span className="meta-label">result_id</span>
                <strong>{displayGraph.metadata.resultId ?? "—"}</strong>
              </div>
              <div>
                <span className="meta-label">映射节点</span>
                <strong>
                  {graphInsights?.mappedNodes ?? 0} / {graphInsights?.totalNodes ?? 0}
                </strong>
              </div>
              <div>
                <span className="meta-label">孤立节点</span>
                <strong>{graphInsights?.isolatedNodes ?? 0}</strong>
              </div>
              <div>
                <span className="meta-label">source</span>
                <strong>{graphInsights?.sourceNodes ?? 0}</strong>
              </div>
              <div>
                <span className="meta-label">sink</span>
                <strong>{graphInsights?.sinkNodes ?? 0}</strong>
              </div>
            </div>
          </div>
        ) : null}

        {graphInsights ? (
          <div className="panel-section">
            <span className="field-label">结构观察</span>
            <div className="detail-card">
              <p>
                当前画布显示 {graphInsights.displayNodes} / {graphInsights.totalNodes} 个节点，
                {graphInsights.displayEdges} / {graphInsights.totalEdges} 条边。
              </p>
              <p>
                {onlyConnectedNodes
                  ? `已隐藏 ${graphInsights.hiddenIsolatedNodes} 个孤立节点。`
                  : "当前保留全部节点，包括孤立节点。"}
              </p>
            </div>
            <span className="field-label">高连接节点</span>
            <ul className="detail-list">
              {graphInsights.topHubs.length ? (
                graphInsights.topHubs.map((node) => (
                  <li key={node.id}>
                    {node.label} · degree {node.degree}
                  </li>
                ))
              ) : (
                <li>无</li>
              )}
            </ul>
          </div>
        ) : null}
      </aside>
    </div>
  );
}

function formatGroupLabel(group: ScanGroup) {
  if (group.time && group.cellType) {
    return `Day ${group.time} · ${group.cellType} · ${group.numDags} DAGs`;
  }
  if (group.cellType) {
    return `${group.cellType} · ${group.numDags} DAGs`;
  }
  return `Group ${group.groupKey} · ${group.numDags} DAGs`;
}

function formatGroupHeading(group: ScanGroup, mode: GraphMode) {
  const base = group.time && group.cellType
    ? `Day ${group.time} · ${group.cellType}`
    : group.cellType || group.groupKey;
  return mode === "single" ? `${base} / 单细胞 DAG` : `${base} / 共识 DAG`;
}

function formatRootLabel(root: RootOption) {
  if (root.collectionCount && root.collectionLabel) {
    return `${root.label} · ${root.collectionCount} 个${root.collectionLabel}`;
  }
  return root.label;
}

function filterGraphByConnectivity(
  graph: GraphResponse | null,
  onlyConnectedNodes: boolean,
): GraphResponse | null {
  if (!graph || !onlyConnectedNodes) {
    return graph;
  }

  const connectedNodeIds = new Set<string>();
  for (const edge of graph.edges) {
    connectedNodeIds.add(edge.source);
    connectedNodeIds.add(edge.target);
  }

  const nodes = graph.nodes.filter((node) => connectedNodeIds.has(node.id));
  const edges = graph.edges.filter(
    (edge) => connectedNodeIds.has(edge.source) && connectedNodeIds.has(edge.target),
  );

  return {
    ...graph,
    nodes,
    edges,
    metadata: {
      ...graph.metadata,
      numNodes: nodes.length,
      numEdges: edges.length,
    },
  };
}

function summarizeGraph(graph: GraphResponse, displayGraph: GraphResponse) {
  const mappedNodes = graph.nodes.filter((node) => node.mapped).length;
  const isolatedNodes = graph.nodes.filter((node) => node.degree === 0).length;
  const sourceNodes = displayGraph.nodes.filter(
    (node) => node.inDegree === 0 && node.outDegree > 0,
  ).length;
  const sinkNodes = displayGraph.nodes.filter(
    (node) => node.outDegree === 0 && node.inDegree > 0,
  ).length;
  const topHubs = [...displayGraph.nodes]
    .sort((left, right) => right.degree - left.degree || left.label.localeCompare(right.label))
    .slice(0, 5);

  return {
    totalNodes: graph.nodes.length,
    totalEdges: graph.edges.length,
    displayNodes: displayGraph.nodes.length,
    displayEdges: displayGraph.edges.length,
    mappedNodes,
    isolatedNodes,
    hiddenIsolatedNodes: Math.max(0, graph.nodes.length - displayGraph.nodes.length),
    sourceNodes,
    sinkNodes,
    topHubs,
  };
}
