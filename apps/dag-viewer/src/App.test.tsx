import { render, screen, waitFor } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import App from "./App";

vi.mock("./components/GraphCanvas", () => ({
  default: ({
    graph,
  }: {
    graph: { metadata?: { groupKey?: string }; nodes?: unknown[] } | null;
  }) => (
    <div data-testid="graph-canvas-mock">
      {graph?.metadata?.groupKey ?? "empty"}:{graph?.nodes?.length ?? 0}
    </div>
  ),
}));

describe("App", () => {
  afterEach(() => {
    vi.restoreAllMocks();
  });

  it("scans the root path and loads the default single graph", async () => {
    const fetchMock = vi.fn(async (input: RequestInfo | URL, init?: RequestInit) => {
      const url = String(input);
      if (url === "/api/roots") {
        return mockJsonResponse({
          roots: [
            {
              path: "/tmp/data/E-GEOD-93593",
              label: "E-GEOD-93593",
              compatible: true,
              rootKind: "dataset",
              collectionKind: "method",
              collectionLabel: "方法",
              collectionCount: 1,
            },
          ],
          defaultRootPath: "/tmp/data/E-GEOD-93593",
          dataRoot: "/tmp/data",
        });
      }
      if (url === "/api/health") {
        return mockJsonResponse({
          dependencies: {
            fastapi: true,
            uvicorn: true,
            pgmpy: true,
            networkx: true,
            pandas: true,
          },
          supportedRootKinds: ["dataset", "dag"],
          frontendBuilt: false,
          frontendDist: "/tmp/dist",
        });
      }
      if (url === "/api/scan") {
        return mockJsonResponse({
          rootPath: "/tmp/data/E-GEOD-93593",
          rootKind: "dataset",
          rootLabel: "E-GEOD-93593",
          methods: ["kTotal"],
          collections: [
            {
              key: "kTotal",
              displayLabel: "kTotal",
              rawLabel: "kTotal",
              kind: "method",
            },
          ],
          collectionKind: "method",
          collectionLabel: "方法",
          showCollectionSelector: false,
          groups: [
            {
              method: "kTotal",
              groupKey: "Day54_cortical_interneuron",
              time: "54",
              cellType: "cortical interneuron",
              numDags: 214,
              hasNodeMap: true,
              hasConsensusCsv: true,
              resultIds: [0, 1, 2],
            },
          ],
        });
      }
      if (url === "/api/graph") {
        const body = JSON.parse(String(init?.body));
        expect(body.mode).toBe("single");
        expect(body.resultId).toBe(0);
        return mockJsonResponse({
          mode: "single",
          nodes: [],
          edges: [],
          metadata: {
            method: "kTotal",
            groupKey: "Day54_cortical_interneuron",
            resultId: 0,
            numNodes: 0,
            numEdges: 0,
          },
        });
      }
      throw new Error(`Unexpected fetch call: ${url}`);
    });

    vi.stubGlobal("fetch", fetchMock);

    render(<App />);

    expect(await screen.findByText("CSCN 可视化")).toBeInTheDocument();
    await screen.findByText("Day 54 · cortical interneuron / 单细胞 DAG");
    await screen.findByText("E-GEOD-93593 · 1 个方法");
    await screen.findByTestId("graph-canvas-mock");
    expect(fetchMock).toHaveBeenCalledWith("/api/scan", expect.any(Object));
  });

  it("switches to consensus mode and requests the consensus graph", async () => {
    const fetchMock = vi.fn(async (input: RequestInfo | URL, init?: RequestInit) => {
      const url = String(input);
      if (url === "/api/roots") {
        return mockJsonResponse({
          roots: [
            {
              path: "/tmp/data/E-GEOD-93593",
              label: "E-GEOD-93593",
              compatible: true,
              rootKind: "dataset",
              collectionKind: "method",
              collectionLabel: "方法",
              collectionCount: 1,
            },
          ],
          defaultRootPath: "/tmp/data/E-GEOD-93593",
          dataRoot: "/tmp/data",
        });
      }
      if (url === "/api/health") {
        return mockJsonResponse({
          dependencies: {
            fastapi: true,
            uvicorn: true,
            pgmpy: true,
            networkx: true,
            pandas: true,
          },
          supportedRootKinds: ["dataset", "dag"],
          frontendBuilt: false,
          frontendDist: "/tmp/dist",
        });
      }
      if (url === "/api/scan") {
        return mockJsonResponse({
          rootPath: "/tmp/data/E-GEOD-93593",
          rootKind: "dataset",
          rootLabel: "E-GEOD-93593",
          methods: ["kTotal"],
          collections: [
            {
              key: "kTotal",
              displayLabel: "kTotal",
              rawLabel: "kTotal",
              kind: "method",
            },
          ],
          collectionKind: "method",
          collectionLabel: "方法",
          showCollectionSelector: false,
          groups: [
            {
              method: "kTotal",
              groupKey: "Day54_cortical_interneuron",
              time: "54",
              cellType: "cortical interneuron",
              numDags: 214,
              hasNodeMap: true,
              hasConsensusCsv: true,
              resultIds: [0, 1, 2],
            },
          ],
        });
      }
      if (url === "/api/graph") {
        const body = JSON.parse(String(init?.body));
        return mockJsonResponse({
          mode: body.mode,
          nodes: [],
          edges: [],
          metadata: {
            method: "kTotal",
            groupKey: "Day54_cortical_interneuron",
            numNodes: 0,
            numEdges: 0,
            threshold: body.mode === "consensus" ? 11 : undefined,
            numDags: body.mode === "consensus" ? 214 : undefined,
          },
        });
      }
      throw new Error(`Unexpected fetch call: ${url}`);
    });

    vi.stubGlobal("fetch", fetchMock);

    render(<App />);

    const user = userEvent.setup();
    const consensusButton = await screen.findByRole("button", { name: "共识" });
    await user.click(consensusButton);

    await waitFor(() => {
      const graphCalls = fetchMock.mock.calls
        .filter(([url]) => String(url) === "/api/graph")
        .map(([, init]) => JSON.parse(String(init?.body)));
      expect(graphCalls.some((payload) => payload.mode === "consensus")).toBe(true);
    });
  });

  it("does not reuse the previous dataset method when switching to a single-view dataset", async () => {
    const fetchMock = vi.fn(async (input: RequestInfo | URL, init?: RequestInit) => {
      const url = String(input);
      if (url === "/api/roots") {
        return mockJsonResponse({
          roots: [
            {
              path: "/tmp/data/E-GEOD-93593",
              label: "E-GEOD-93593",
              compatible: true,
              rootKind: "dataset",
              collectionKind: "method",
              collectionLabel: "方法",
              collectionCount: 1,
            },
            {
              path: "/tmp/data/GSE121893",
              label: "GSE121893",
              compatible: true,
              rootKind: "dataset",
              collectionKind: "single",
              collectionLabel: "当前视图",
              collectionCount: 1,
            },
          ],
          defaultRootPath: "/tmp/data/E-GEOD-93593",
          dataRoot: "/tmp/data",
        });
      }
      if (url === "/api/health") {
        return mockJsonResponse({
          dependencies: {
            fastapi: true,
            uvicorn: true,
            pgmpy: true,
            networkx: true,
            pandas: true,
          },
          supportedRootKinds: ["dataset", "dag"],
          frontendBuilt: false,
          frontendDist: "/tmp/dist",
        });
      }
      if (url === "/api/scan") {
        const body = JSON.parse(String(init?.body));
        if (body.rootPath === "/tmp/data/E-GEOD-93593") {
          return mockJsonResponse({
            rootPath: "/tmp/data/E-GEOD-93593",
            rootKind: "dataset",
            rootLabel: "E-GEOD-93593",
            methods: ["kTotal"],
            collections: [
              {
                key: "kTotal",
                displayLabel: "kTotal",
                rawLabel: "kTotal",
                kind: "method",
              },
            ],
            collectionKind: "method",
            collectionLabel: "方法",
            showCollectionSelector: false,
            groups: [
              {
                method: "kTotal",
                groupKey: "Day54_cortical_interneuron",
                time: "54",
                cellType: "cortical interneuron",
                numDags: 214,
                hasNodeMap: true,
                hasConsensusCsv: true,
                resultIds: [0],
              },
            ],
          });
        }
        return mockJsonResponse({
          rootPath: "/tmp/data/GSE121893",
          rootKind: "dataset",
          rootLabel: "GSE121893",
          methods: ["GSE121893_dhf_vs_n_all_all"],
          collections: [
            {
              key: "GSE121893_dhf_vs_n_all_all",
              displayLabel: "默认视图",
              rawLabel: "GSE121893_dhf_vs_n_all_all",
              kind: "single",
            },
          ],
          collectionKind: "single",
          collectionLabel: "当前视图",
          showCollectionSelector: false,
          groups: [
            {
              method: "GSE121893_dhf_vs_n_all_all",
              groupKey: "dHF",
              time: "",
              cellType: "dHF",
              numDags: 10,
              hasNodeMap: true,
              hasConsensusCsv: false,
              resultIds: [0],
            },
          ],
        });
      }
      if (url === "/api/graph") {
        const body = JSON.parse(String(init?.body));
        if (body.rootPath === "/tmp/data/E-GEOD-93593") {
          return mockJsonResponse({
            mode: "single",
            nodes: [],
            edges: [],
            metadata: {
              method: "kTotal",
              groupKey: "Day54_cortical_interneuron",
              resultId: 0,
              numNodes: 0,
              numEdges: 0,
            },
          });
        }
        return mockJsonResponse({
          mode: "single",
          nodes: [],
          edges: [],
          metadata: {
            method: "GSE121893_dhf_vs_n_all_all",
            groupKey: "dHF",
            resultId: 0,
            numNodes: 0,
            numEdges: 0,
          },
        });
      }
      throw new Error(`Unexpected fetch call: ${url}`);
    });

    vi.stubGlobal("fetch", fetchMock);

    render(<App />);

    const user = userEvent.setup();
    await screen.findByText("Day 54 · cortical interneuron / 单细胞 DAG");

    await user.selectOptions(screen.getByLabelText("数据集"), "/tmp/data/GSE121893");

    await screen.findByText("dHF / 单细胞 DAG");
    expect(screen.queryByText(/加载图失败/)).not.toBeInTheDocument();

    const graphCalls = fetchMock.mock.calls
      .filter(([url]) => String(url) === "/api/graph")
      .map(([, options]) => JSON.parse(String(options?.body)));
    const gseCalls = graphCalls.filter(
      (payload) => payload.rootPath === "/tmp/data/GSE121893",
    );

    expect(gseCalls.length).toBeGreaterThan(0);
    expect(
      gseCalls.every((payload) => payload.method === "GSE121893_dhf_vs_n_all_all"),
    ).toBe(true);
  });

  it("hides the collection selector for single-run datasets and shows a dataset list selector", async () => {
    const fetchMock = vi.fn(async (input: RequestInfo | URL) => {
      const url = String(input);
      if (url === "/api/roots") {
        return mockJsonResponse({
          roots: [
            {
              path: "/tmp/data/GSE121893",
              label: "GSE121893",
              compatible: true,
              rootKind: "dataset",
              collectionKind: "single",
              collectionLabel: "当前视图",
              collectionCount: 1,
            },
            {
              path: "/tmp/data/GSE138852",
              label: "GSE138852",
              compatible: false,
              reason: "DAG 根目录为空",
            },
          ],
          defaultRootPath: "/tmp/data/GSE121893",
          dataRoot: "/tmp/data",
        });
      }
      if (url === "/api/health") {
        return mockJsonResponse({
          dependencies: {
            fastapi: true,
            uvicorn: true,
            pgmpy: true,
            networkx: true,
            pandas: true,
          },
          supportedRootKinds: ["dataset", "dag"],
          frontendBuilt: false,
          frontendDist: "/tmp/dist",
        });
      }
      if (url === "/api/scan") {
        return mockJsonResponse({
          rootPath: "/tmp/data/GSE121893",
          rootKind: "dataset",
          rootLabel: "GSE121893",
          methods: ["GSE121893_dhf_vs_n_all_all"],
          collections: [
            {
              key: "GSE121893_dhf_vs_n_all_all",
              displayLabel: "默认视图",
              rawLabel: "GSE121893_dhf_vs_n_all_all",
              kind: "single",
            },
          ],
          collectionKind: "single",
          collectionLabel: "当前视图",
          showCollectionSelector: false,
          groups: [
            {
              method: "GSE121893_dhf_vs_n_all_all",
              groupKey: "dHF",
              time: "",
              cellType: "dHF",
              numDags: 10,
              hasNodeMap: true,
              hasConsensusCsv: false,
              resultIds: [0, 1],
            },
          ],
        });
      }
      if (url === "/api/graph") {
        return mockJsonResponse({
          mode: "single",
          nodes: [],
          edges: [],
          metadata: {
            method: "GSE121893_dhf_vs_n_all_all",
            groupKey: "dHF",
            numNodes: 0,
            numEdges: 0,
          },
        });
      }
      throw new Error(`Unexpected fetch call: ${url}`);
    });

    vi.stubGlobal("fetch", fetchMock);

    render(<App />);

    expect(await screen.findByLabelText("数据集")).toBeInTheDocument();
    expect((await screen.findAllByText("当前视图")).length).toBeGreaterThan(0);
    expect(await screen.findByRole("heading", { name: "默认视图", level: 4 })).toBeInTheDocument();
    expect(screen.getByText("GSE121893_dhf_vs_n_all_all")).toBeInTheDocument();
    expect(screen.getByRole("option", { name: "dHF · 10 DAGs" })).toBeInTheDocument();
  });

  it("filters isolated nodes when only connected nodes are requested", async () => {
    const fetchMock = vi.fn(async (input: RequestInfo | URL, init?: RequestInit) => {
      const url = String(input);
      if (url === "/api/roots") {
        return mockJsonResponse({
          roots: [
            {
              path: "/tmp/data/E-GEOD-93593",
              label: "E-GEOD-93593",
              compatible: true,
              rootKind: "dataset",
              collectionKind: "method",
              collectionLabel: "方法",
              collectionCount: 1,
            },
          ],
          defaultRootPath: "/tmp/data/E-GEOD-93593",
          dataRoot: "/tmp/data",
        });
      }
      if (url === "/api/health") {
        return mockJsonResponse({
          dependencies: {
            fastapi: true,
            uvicorn: true,
            pgmpy: true,
            networkx: true,
            pandas: true,
          },
          supportedRootKinds: ["dataset", "dag"],
          frontendBuilt: false,
          frontendDist: "/tmp/dist",
        });
      }
      if (url === "/api/scan") {
        return mockJsonResponse({
          rootPath: "/tmp/data/E-GEOD-93593",
          rootKind: "dataset",
          rootLabel: "E-GEOD-93593",
          methods: ["kTotal"],
          collections: [
            {
              key: "kTotal",
              displayLabel: "kTotal",
              rawLabel: "kTotal",
              kind: "method",
            },
          ],
          collectionKind: "method",
          collectionLabel: "方法",
          showCollectionSelector: false,
          groups: [
            {
              method: "kTotal",
              groupKey: "Day54_cortical_interneuron",
              time: "54",
              cellType: "cortical interneuron",
              numDags: 214,
              hasNodeMap: true,
              hasConsensusCsv: true,
              resultIds: [0],
            },
          ],
        });
      }
      if (url === "/api/graph") {
        const body = JSON.parse(String(init?.body));
        expect(body.mode).toBe("single");
        return mockJsonResponse({
          mode: "single",
          nodes: [
            {
              id: "A",
              label: "A",
              rawId: "A",
              mapped: true,
              degree: 1,
              inDegree: 0,
              outDegree: 1,
            },
            {
              id: "B",
              label: "B",
              rawId: "B",
              mapped: true,
              degree: 1,
              inDegree: 1,
              outDegree: 0,
            },
            {
              id: "C",
              label: "C",
              rawId: "C",
              mapped: true,
              degree: 0,
              inDegree: 0,
              outDegree: 0,
            },
          ],
          edges: [
            {
              id: "A->B",
              source: "A",
              target: "B",
              count: 1,
              weight: 1,
            },
          ],
          metadata: {
            method: "kTotal",
            groupKey: "Day54_cortical_interneuron",
            resultId: 0,
            numNodes: 3,
            numEdges: 1,
          },
        });
      }
      throw new Error(`Unexpected fetch call: ${url}`);
    });

    vi.stubGlobal("fetch", fetchMock);

    render(<App />);

    const user = userEvent.setup();
    const connectedOnlyCheckbox = await screen.findByRole("checkbox", {
      name: "仅保留有边相连的节点",
    });
    expect(connectedOnlyCheckbox).toBeChecked();
    await screen.findByText("2 nodes");
    expect(screen.getByText("已隐藏孤立节点")).toBeInTheDocument();
    expect(screen.getByTestId("graph-canvas-mock")).toHaveTextContent(
      "Day54_cortical_interneuron:2",
    );

    await user.click(connectedOnlyCheckbox);

    await screen.findByText("3 nodes");
    expect(screen.getByTestId("graph-canvas-mock")).toHaveTextContent(
      "Day54_cortical_interneuron:3",
    );
  });
});

function mockJsonResponse(payload: unknown) {
  return Promise.resolve({
    ok: true,
    json: async () => payload,
  } as Response);
}
