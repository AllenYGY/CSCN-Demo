from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.preprocessing import StandardScaler

from cscn.core import CSCN
from cscn.layout import safe_component

try:
    import umap
except ModuleNotFoundError:  # pragma: no cover - optional dependency
    umap = None

try:
    import seaborn as sns
except ModuleNotFoundError:  # pragma: no cover - optional dependency
    sns = None


def _log(message: str) -> None:
    print(f"[seqfish-compare] {message}")


def _load_gene_names(run_dir: Path) -> list[str]:
    genes_path = run_dir / "inputs" / "genes.csv"
    frame = pd.read_csv(genes_path)
    return frame["gene_name"].astype(str).tolist()


def _load_group_keys(run_dir: Path) -> list[str]:
    groups_path = run_dir / "inputs" / "groups.csv"
    frame = pd.read_csv(groups_path)
    return frame["group_key"].astype(str).tolist()


def _ckm_path(run_dir: Path, group_key: str) -> Path:
    return run_dir / "ckm" / f"{safe_component(group_key)}_ckm.npy"


def _group_cells_path(run_dir: Path, group_key: str) -> Path:
    return run_dir / "matrices" / f"{safe_component(group_key)}_cells.csv"


def _group_object_path(run_dir: Path, group_key: str) -> Path:
    return run_dir / "objects" / f"{safe_component(group_key)}_cscn.pkl"


def ensure_group_ckm(
    run_dir: Path,
    group_key: str,
    *,
    alpha: float = 0.05,
    beta_transform: str = "log1p",
    strict: bool = True,
) -> np.ndarray:
    ckm_path = _ckm_path(run_dir, group_key)
    if ckm_path.is_file():
        return np.load(ckm_path)

    object_path = _group_object_path(run_dir, group_key)
    if not object_path.is_file():
        raise FileNotFoundError(
            f"Missing CKM and CSCN object for group {group_key}: {ckm_path} / {object_path}"
        )

    ckm_path.parent.mkdir(parents=True, exist_ok=True)
    cscn = CSCN.load_from_file(object_path)
    return cscn.compute_ckm(
        alpha=alpha,
        beta_transform=beta_transform,
        save_path=str(ckm_path),
        strict=strict,
    )


def build_pooled_ckm_df(
    run_dir: Path,
    metadata: pd.DataFrame,
    *,
    alpha: float = 0.05,
    beta_transform: str = "log1p",
    strict: bool = True,
) -> pd.DataFrame:
    gene_names = _load_gene_names(run_dir)
    frames: list[pd.DataFrame] = []

    for group_key in _load_group_keys(run_dir):
        cells_path = _group_cells_path(run_dir, group_key)
        if not cells_path.is_file():
            raise FileNotFoundError(f"Missing group cell list for {group_key}: {cells_path}")
        cells = pd.read_csv(cells_path)["cell_id"].astype(str).tolist()
        ckm = ensure_group_ckm(
            run_dir,
            group_key,
            alpha=alpha,
            beta_transform=beta_transform,
            strict=strict,
        )
        if ckm.shape[0] != len(cells):
            raise ValueError(
                f"CKM row count does not match cells for group {group_key}: {ckm.shape[0]} vs {len(cells)}"
            )
        if ckm.shape[1] != len(gene_names):
            raise ValueError(
                f"CKM column count does not match genes for group {group_key}: {ckm.shape[1]} vs {len(gene_names)}"
            )
        frames.append(pd.DataFrame(ckm, index=cells, columns=gene_names))

    pooled = pd.concat(frames, axis=0)
    if pooled.index.has_duplicates:
        duplicates = pooled.index[pooled.index.duplicated()].unique().tolist()
        raise ValueError(f"Duplicated cell ids across CKM blocks: {duplicates[:10]}")

    missing = metadata.index.difference(pooled.index).tolist()
    extra = pooled.index.difference(metadata.index).tolist()
    if missing or extra:
        raise ValueError(
            f"CKM matrix could not align exactly to metadata. Missing={missing[:10]} Extra={extra[:10]}"
        )
    return pooled.loc[metadata.index].copy()


def build_expression_df(expr_path: Path, metadata: pd.DataFrame, gene_names: list[str]) -> pd.DataFrame:
    expr = pd.read_csv(expr_path)
    if "cell_id" not in expr.columns:
        raise ValueError("Expression baseline requires a `cell_id` column.")
    expr["cell_id"] = expr["cell_id"].astype(str)
    expr = expr.set_index("cell_id")
    missing_genes = [gene for gene in gene_names if gene not in expr.columns]
    if missing_genes:
        raise ValueError(f"Expression baseline is missing run genes: {missing_genes[:10]}")
    missing = metadata.index.difference(expr.index).tolist()
    extra = expr.index.difference(metadata.index).tolist()
    if missing or extra:
        raise ValueError(
            f"Expression matrix could not align exactly to metadata. Missing={missing[:10]} Extra={extra[:10]}"
        )
    return expr.loc[metadata.index, gene_names].astype(float).copy()


def _plot_embedding(
    embedding: np.ndarray,
    labels: pd.Series,
    title: str,
    output_path: Path,
    axis_prefix: str,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(
        {
            f"{axis_prefix}1": embedding[:, 0],
            f"{axis_prefix}2": embedding[:, 1],
            "label": labels.astype(str).tolist(),
        }
    )
    plt.figure(figsize=(8, 6))
    if sns is not None:
        sns.scatterplot(
            data=frame,
            x=f"{axis_prefix}1",
            y=f"{axis_prefix}2",
            hue="label",
            s=24,
            linewidth=0,
            alpha=0.85,
        )
    else:
        for label, group in frame.groupby("label", sort=False):
            plt.scatter(
                group[f"{axis_prefix}1"],
                group[f"{axis_prefix}2"],
                s=24,
                alpha=0.85,
                label=str(label),
            )
        plt.legend()
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def analyze_representation(
    name: str,
    features: pd.DataFrame,
    truth_labels: pd.Series,
    output_dir: Path,
    *,
    random_seed: int = 42,
    n_pcs: int = 2,
    enable_umap: bool = True,
) -> tuple[dict[str, float | int | str], pd.Series]:
    output_dir.mkdir(parents=True, exist_ok=True)
    features.to_csv(output_dir / f"{name}_features.csv")

    scaler = StandardScaler()
    scaled = scaler.fit_transform(features.to_numpy(dtype=float))

    n_clusters = int(truth_labels.nunique())
    kmeans = KMeans(n_clusters=n_clusters, random_state=random_seed, n_init=20)
    clusters = pd.Series(kmeans.fit_predict(scaled), index=features.index, name=f"{name}_cluster")

    pca = PCA(n_components=max(2, min(n_pcs, scaled.shape[1], scaled.shape[0])))
    pca_embedding = pca.fit_transform(scaled)
    if pca_embedding.shape[1] > 2:
        pca_embedding = pca_embedding[:, :2]
    elif pca_embedding.shape[1] == 1:
        pca_embedding = np.column_stack([pca_embedding[:, 0], np.zeros(len(pca_embedding))])

    _plot_embedding(
        pca_embedding,
        truth_labels,
        f"{name} PCA colored by truth",
        output_dir / f"{name}_pca_truth.png",
        "PC",
    )
    _plot_embedding(
        pca_embedding,
        clusters.astype(str),
        f"{name} PCA colored by cluster",
        output_dir / f"{name}_pca_cluster.png",
        "PC",
    )

    if enable_umap and umap is not None:
        reducer = umap.UMAP(random_state=random_seed)
        umap_embedding = reducer.fit_transform(scaled)
        _plot_embedding(
            umap_embedding,
            truth_labels,
            f"{name} UMAP colored by truth",
            output_dir / f"{name}_umap_truth.png",
            "UMAP",
        )
        _plot_embedding(
            umap_embedding,
            clusters.astype(str),
            f"{name} UMAP colored by cluster",
            output_dir / f"{name}_umap_cluster.png",
            "UMAP",
        )
    elif enable_umap:
        _log("UMAP is unavailable; skipping UMAP plots.")

    metrics = {
        "representation": name,
        "n_cells": int(features.shape[0]),
        "n_features": int(features.shape[1]),
        "n_clusters": n_clusters,
        "ari": float(adjusted_rand_score(truth_labels, clusters)),
        "nmi": float(normalized_mutual_info_score(truth_labels, clusters)),
        "silhouette": float(silhouette_score(scaled, clusters)),
    }
    return metrics, clusters


def run_analysis(
    expr_path: Path,
    metadata_path: Path,
    weighted_run_dir: Path,
    local_knn_run_dir: Path,
    output_dir: Path,
    *,
    random_seed: int = 42,
    enable_umap: bool = True,
    n_pcs: int = 2,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(metadata_path)
    if "cell_id" not in metadata.columns or "cell_class_name" not in metadata.columns:
        raise ValueError("Metadata must contain `cell_id` and `cell_class_name`.")
    metadata["cell_id"] = metadata["cell_id"].astype(str)
    metadata = metadata.drop_duplicates(subset=["cell_id"]).set_index("cell_id")
    truth_labels = metadata["cell_class_name"].astype(str)

    weighted_genes = _load_gene_names(weighted_run_dir)
    local_genes = _load_gene_names(local_knn_run_dir)
    if weighted_genes != local_genes:
        raise ValueError("Weighted and local-KNN runs do not share the same gene list.")

    expr_df = build_expression_df(expr_path, metadata, weighted_genes)
    weighted_ckm = build_pooled_ckm_df(weighted_run_dir, metadata)
    local_ckm = build_pooled_ckm_df(local_knn_run_dir, metadata)

    representations = {
        "expr": expr_df,
        "ckm_weighted": weighted_ckm,
        "ckm_local_knn": local_ckm,
    }

    metrics_rows: list[dict[str, float | int | str]] = []
    assignment_frame = pd.DataFrame(
        {
            "cell_id": metadata.index,
            "cell_class_name": truth_labels.values,
        }
    ).set_index("cell_id")

    for name, frame in representations.items():
        metrics, clusters = analyze_representation(
            name,
            frame,
            truth_labels,
            output_dir,
            random_seed=random_seed,
            n_pcs=n_pcs,
            enable_umap=enable_umap,
        )
        metrics_rows.append(metrics)
        assignment_frame[clusters.name] = clusters

    metrics_df = pd.DataFrame(metrics_rows)
    metrics_df.to_csv(output_dir / "clustering_metrics.csv", index=False)
    assignment_frame.reset_index().to_csv(output_dir / "cell_assignments.csv", index=False)
    _log(f"saved analysis outputs to {output_dir}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compare seqFISH expression and CKM clustering.")
    parser.add_argument("--expr-path", type=Path, required=True)
    parser.add_argument("--metadata-path", type=Path, required=True)
    parser.add_argument("--weighted-run-dir", type=Path, required=True)
    parser.add_argument("--local-knn-run-dir", type=Path, required=True)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("config/seqfish/analysis_ckm_compare"),
    )
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--n-pcs", type=int, default=2)
    parser.add_argument("--umap", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    run_analysis(
        expr_path=args.expr_path,
        metadata_path=args.metadata_path,
        weighted_run_dir=args.weighted_run_dir,
        local_knn_run_dir=args.local_knn_run_dir,
        output_dir=args.output_dir,
        random_seed=args.random_seed,
        enable_umap=args.umap,
        n_pcs=args.n_pcs,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
