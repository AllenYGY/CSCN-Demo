export type GraphMode = "single" | "consensus";

export interface DependencyStatus {
  fastapi: boolean;
  uvicorn: boolean;
  pgmpy: boolean;
  networkx: boolean;
  pandas: boolean;
}

export interface HealthResponse {
  dependencies: DependencyStatus;
  supportedRootKinds: string[];
  frontendBuilt: boolean;
  frontendDist: string;
}

export interface RootOption {
  path: string;
  label: string;
  compatible: boolean;
  reason?: string;
  rootKind?: string;
  collectionKind?: string;
  collectionLabel?: string;
  collectionCount?: number;
}

export interface RootsResponse {
  roots: RootOption[];
  defaultRootPath: string | null;
  dataRoot: string;
}

export interface CollectionOption {
  key: string;
  displayLabel: string;
  rawLabel: string;
  kind: string;
}

export interface ScanGroup {
  method: string;
  groupKey: string;
  time: string;
  cellType: string;
  numDags: number;
  hasNodeMap: boolean;
  hasConsensusCsv: boolean;
  resultIds: number[];
}

export interface ScanResponse {
  rootPath: string;
  rootKind: string;
  rootLabel: string;
  methods: string[];
  collections: CollectionOption[];
  collectionKind: string;
  collectionLabel: string;
  showCollectionSelector: boolean;
  groups: ScanGroup[];
}

export interface GraphNode {
  id: string;
  label: string;
  rawId: string;
  mapped: boolean;
  degree: number;
  inDegree: number;
  outDegree: number;
}

export interface GraphEdge {
  id: string;
  source: string;
  target: string;
  count: number;
  weight: number;
}

export interface GraphMetadata {
  method: string;
  groupKey: string;
  resultId?: number;
  cellRunId?: string;
  numNodes: number;
  numEdges: number;
  threshold?: number;
  numDags?: number;
}

export interface GraphResponse {
  mode: GraphMode;
  nodes: GraphNode[];
  edges: GraphEdge[];
  metadata: GraphMetadata;
}
