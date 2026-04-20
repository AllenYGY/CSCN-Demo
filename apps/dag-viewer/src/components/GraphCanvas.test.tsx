import { getVisibleLabelIds } from "./GraphCanvas";
import type { GraphResponse } from "../types";

function buildGraph(nodeCount: number): GraphResponse {
  return {
    mode: "single",
    nodes: Array.from({ length: nodeCount }, (_, index) => ({
      id: `n${index}`,
      label: `Gene-${index}`,
      rawId: `raw-${index}`,
      mapped: true,
      degree: index === 0 ? 2 : 1,
      inDegree: index === 0 ? 1 : 0,
      outDegree: index === 0 ? 1 : 1,
    })),
    edges: [
      {
        id: "n0->n1",
        source: "n0",
        target: "n1",
        count: 1,
        weight: 1,
      },
    ],
    metadata: {
      method: "kTotal",
      groupKey: "Day54_cortical_interneuron",
      numNodes: nodeCount,
      numEdges: 1,
    },
  };
}

describe("getVisibleLabelIds", () => {
  it("shows all labels for compact graphs", () => {
    const graph = buildGraph(8);
    const visible = getVisibleLabelIds({
      graph,
      selectedNodeId: null,
      searchText: "",
      showAllLabels: false,
    });
    expect(visible.size).toBe(8);
  });

  it("limits labels to search matches and selected neighborhoods on larger graphs", () => {
    const graph = buildGraph(120);
    const visible = getVisibleLabelIds({
      graph,
      selectedNodeId: "n0",
      searchText: "Gene-35",
      showAllLabels: false,
    });

    expect(visible.has("n0")).toBe(true);
    expect(visible.has("n1")).toBe(true);
    expect(visible.has("n35")).toBe(true);
    expect(visible.has("n80")).toBe(false);
  });
});
