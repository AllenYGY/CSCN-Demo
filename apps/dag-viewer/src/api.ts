import type {
  GraphMode,
  GraphResponse,
  HealthResponse,
  RootsResponse,
  ScanResponse,
} from "./types";

async function fetchJson<T>(path: string, init?: RequestInit): Promise<T> {
  const response = await fetch(path, {
    headers: {
      "Content-Type": "application/json",
      ...(init?.headers ?? {}),
    },
    ...init,
  });

  if (!response.ok) {
    let detail = `${response.status} ${response.statusText}`;
    try {
      const data = (await response.json()) as { detail?: string };
      if (data.detail) {
        detail = data.detail;
      }
    } catch {
      // Ignore JSON parse failures and keep the original status text.
    }
    throw new Error(detail);
  }

  return (await response.json()) as T;
}

export function fetchHealth(): Promise<HealthResponse> {
  return fetchJson<HealthResponse>("/api/health");
}

export function fetchRoots(): Promise<RootsResponse> {
  return fetchJson<RootsResponse>("/api/roots");
}

export function scanDagRoot(rootPath: string): Promise<ScanResponse> {
  return fetchJson<ScanResponse>("/api/scan", {
    method: "POST",
    body: JSON.stringify({ rootPath }),
  });
}

export function fetchGraph(args: {
  rootPath: string;
  method: string;
  groupKey: string;
  mode: GraphMode;
  resultId?: number;
}): Promise<GraphResponse> {
  return fetchJson<GraphResponse>("/api/graph", {
    method: "POST",
    body: JSON.stringify(args),
  });
}
