import { useDeferredValue, useEffect, useRef } from "react";
import cytoscape, {
  type Core,
  type EdgeSingular,
  type LayoutOptions,
  type NodeSingular,
} from "cytoscape";
import dagre from "cytoscape-dagre";
import type { GraphResponse } from "../types";

cytoscape.use(dagre);

export interface GraphCanvasProps {
  graph: GraphResponse | null;
  searchText: string;
  showAllLabels: boolean;
  selectedNodeId: string | null;
  selectedEdgeId: string | null;
  onSelectNode: (nodeId: string) => void;
  onSelectEdge: (edgeId: string) => void;
  onClearSelection: () => void;
}

export function getVisibleLabelIds(args: {
  graph: GraphResponse | null;
  selectedNodeId: string | null;
  searchText: string;
  showAllLabels: boolean;
}): Set<string> {
  const { graph, selectedNodeId, searchText, showAllLabels } = args;
  if (!graph) {
    return new Set<string>();
  }

  if (graph.nodes.length <= 80 || showAllLabels) {
    return new Set(graph.nodes.map((node) => node.id));
  }

  const keep = new Set<string>();
  const query = searchText.trim().toLowerCase();
  if (query) {
    for (const node of graph.nodes) {
      if (
        node.label.toLowerCase().includes(query) ||
        node.rawId.toLowerCase().includes(query)
      ) {
        keep.add(node.id);
      }
    }
  }

  if (selectedNodeId) {
    keep.add(selectedNodeId);
    for (const edge of graph.edges) {
      if (edge.source === selectedNodeId) {
        keep.add(edge.target);
      }
      if (edge.target === selectedNodeId) {
        keep.add(edge.source);
      }
    }
  }

  return keep;
}

export default function GraphCanvas({
  graph,
  searchText,
  showAllLabels,
  selectedNodeId,
  selectedEdgeId,
  onSelectNode,
  onSelectEdge,
  onClearSelection,
}: GraphCanvasProps) {
  const containerRef = useRef<HTMLDivElement | null>(null);
  const cyRef = useRef<Core | null>(null);
  const selectNodeRef = useRef(onSelectNode);
  const selectEdgeRef = useRef(onSelectEdge);
  const clearSelectionRef = useRef(onClearSelection);
  const deferredSearch = useDeferredValue(searchText);

  useEffect(() => {
    selectNodeRef.current = onSelectNode;
    selectEdgeRef.current = onSelectEdge;
    clearSelectionRef.current = onClearSelection;
  }, [onSelectNode, onSelectEdge, onClearSelection]);

  useEffect(() => {
    if (!containerRef.current) {
      return;
    }

    const cy = cytoscape({
      container: containerRef.current,
      minZoom: 0.25,
      maxZoom: 3,
      wheelSensitivity: 0.18,
      style: [
        {
          selector: "node",
          style: {
            "background-color": "data(fillColor)",
            "border-color": "#14333a",
            "border-width": 1.5,
            width: "data(size)",
            height: "data(size)",
            label: "data(displayLabel)",
            "font-size": 12,
            "font-family": "Avenir Next, Segoe UI, sans-serif",
            color: "#18313a",
            "text-wrap": "wrap",
            "text-max-width": "90px",
            "text-valign": "center",
            "text-halign": "center",
            "overlay-opacity": 0,
          },
        },
        {
          selector: "edge",
          style: {
            width: "data(lineWidth)",
            "line-color": "data(lineColor)",
            "target-arrow-shape": "triangle",
            "target-arrow-color": "data(lineColor)",
            "curve-style": "bezier",
            opacity: 0.9,
            "arrow-scale": 1.05,
          },
        },
        {
          selector: ".muted",
          style: {
            opacity: 0.14,
          },
        },
        {
          selector: ".selected-node",
          style: {
            "border-width": 3,
            "border-color": "#7b2f19",
          },
        },
        {
          selector: ".upstream-node",
          style: {
            "background-color": "#c8744d",
          },
        },
        {
          selector: ".downstream-node",
          style: {
            "background-color": "#4fa38d",
          },
        },
        {
          selector: ".search-match",
          style: {
            "border-width": 3,
            "border-color": "#d8a543",
          },
        },
        {
          selector: ".selected-edge",
          style: {
            width: 6,
            "line-color": "#7b2f19",
            "target-arrow-color": "#7b2f19",
          },
        },
        {
          selector: ".focus-edge-upstream",
          style: {
            width: 5,
            "line-color": "#c8744d",
            "target-arrow-color": "#c8744d",
          },
        },
        {
          selector: ".focus-edge-downstream",
          style: {
            width: 5,
            "line-color": "#4fa38d",
            "target-arrow-color": "#4fa38d",
          },
        },
      ],
    });

    cy.on("tap", "node", (event) => {
      const node = event.target as NodeSingular;
      selectNodeRef.current(node.id());
    });

    cy.on("tap", "edge", (event) => {
      const edge = event.target as EdgeSingular;
      selectEdgeRef.current(edge.id());
    });

    cy.on("tap", (event) => {
      if (event.target === cy) {
        clearSelectionRef.current();
      }
    });

    cyRef.current = cy;
    return () => {
      cy.destroy();
      cyRef.current = null;
    };
  }, []);

  useEffect(() => {
    const cy = cyRef.current;
    if (!cy) {
      return;
    }

    cy.elements().remove();
    if (!graph) {
      return;
    }

    cy.add(buildElements(graph));
    const layout = cy.layout({
      name: "dagre",
      padding: 36,
      rankDir: "LR",
      nodeSep: 30,
      rankSep: 100,
      edgeSep: 18,
      animate: false,
    } as LayoutOptions);
    layout.run();
    cy.fit(undefined, 36);
  }, [graph]);

  useEffect(() => {
    const cy = cyRef.current;
    if (!cy || !graph) {
      return;
    }

    applyInteractionState({
      cy,
      graph,
      searchText: deferredSearch,
      showAllLabels,
      selectedNodeId,
      selectedEdgeId,
    });
  }, [graph, deferredSearch, showAllLabels, selectedNodeId, selectedEdgeId]);

  return <div ref={containerRef} className="graph-canvas" data-testid="graph-canvas" />;
}

function buildElements(graph: GraphResponse) {
  const maxDegree = Math.max(...graph.nodes.map((node) => node.degree), 1);
  const maxWeight = Math.max(...graph.edges.map((edge) => edge.weight), 1);

  return [
    ...graph.nodes.map((node) => {
      const degreeRatio = node.degree / maxDegree;
      return {
        data: {
          ...node,
          displayLabel: node.label,
          size: 24 + degreeRatio * 26,
          fillColor: interpolateColor("#d8e7dd", "#2d7663", degreeRatio),
        },
      };
    }),
    ...graph.edges.map((edge) => {
      const weightRatio = edge.weight / maxWeight;
      return {
        data: {
          ...edge,
          lineWidth: 1.5 + weightRatio * 4.5,
          lineColor: interpolateColor("#c7d4d9", "#7b2f19", weightRatio),
        },
      };
    }),
  ];
}

function applyInteractionState(args: {
  cy: Core;
  graph: GraphResponse;
  searchText: string;
  showAllLabels: boolean;
  selectedNodeId: string | null;
  selectedEdgeId: string | null;
}) {
  const { cy, graph, searchText, showAllLabels, selectedNodeId, selectedEdgeId } = args;
  const visibleLabels = getVisibleLabelIds({
    graph,
    selectedNodeId,
    searchText,
    showAllLabels,
  });

  cy.elements().removeClass(
    "muted selected-node upstream-node downstream-node search-match selected-edge focus-edge-upstream focus-edge-downstream",
  );

  cy.nodes().forEach((node) => {
    const label = String(node.data("label") ?? "");
    node.data("displayLabel", visibleLabels.has(node.id()) ? label : "");
  });

  const query = searchText.trim().toLowerCase();
  if (query) {
    cy.nodes().forEach((node) => {
      const label = String(node.data("label") ?? "").toLowerCase();
      const rawId = String(node.data("rawId") ?? "").toLowerCase();
      if (label.includes(query) || rawId.includes(query)) {
        node.addClass("search-match");
      }
    });
  }

  if (selectedNodeId) {
    const selectedNode = cy.getElementById(selectedNodeId);
    if (selectedNode.nonempty()) {
      const upstreamEdges = selectedNode.incomers("edge");
      const upstreamNodes = upstreamEdges.sources();
      const downstreamEdges = selectedNode.outgoers("edge");
      const downstreamNodes = downstreamEdges.targets();
      const focus = selectedNode
        .union(upstreamEdges)
        .union(upstreamNodes)
        .union(downstreamEdges)
        .union(downstreamNodes);

      cy.elements().difference(focus).addClass("muted");
      selectedNode.addClass("selected-node");
      upstreamNodes.addClass("upstream-node");
      downstreamNodes.addClass("downstream-node");
      upstreamEdges.addClass("focus-edge-upstream");
      downstreamEdges.addClass("focus-edge-downstream");
      return;
    }
  }

  if (selectedEdgeId) {
    const selectedEdge = cy.getElementById(selectedEdgeId);
    if (selectedEdge.nonempty()) {
      const focus = selectedEdge.union(selectedEdge.source()).union(selectedEdge.target());
      cy.elements().difference(focus).addClass("muted");
      selectedEdge.addClass("selected-edge");
      selectedEdge.source().addClass("selected-node");
      selectedEdge.target().addClass("downstream-node");
    }
  }
}

function interpolateColor(start: string, end: string, ratio: number) {
  const clamp = Math.max(0, Math.min(1, ratio));
  const startRgb = hexToRgb(start);
  const endRgb = hexToRgb(end);
  const mixed = startRgb.map((value, index) =>
    Math.round(value + (endRgb[index] - value) * clamp),
  );
  return `rgb(${mixed[0]}, ${mixed[1]}, ${mixed[2]})`;
}

function hexToRgb(hex: string) {
  const normalized = hex.replace("#", "");
  return [
    Number.parseInt(normalized.slice(0, 2), 16),
    Number.parseInt(normalized.slice(2, 4), 16),
    Number.parseInt(normalized.slice(4, 6), 16),
  ];
}
