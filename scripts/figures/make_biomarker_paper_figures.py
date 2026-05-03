from __future__ import annotations

import gzip
import os
import textwrap
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp") / "matplotlib-cscn"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = REPO_ROOT / "data"
OUT_ROOT = REPO_ROOT / "results" / "paper_figures"


DATASETS = {
    "SCP259": {
        "title": "SCP259",
        "out_dir": OUT_ROOT / "SCP259",
        "biomarker_csv": DATA_ROOT / "SCP259" / "Biomarkers_SCP259_inflamed_vs_healthy_crypt_prolif_epi.csv",
        "biomarker_go": DATA_ROOT / "SCP259" / "enrichment_results" / "Biomarkers" / "SCP259_inflamed_vs_healthy_crypt_prolif_epi_Gene_Biomarkers_GO_results.csv",
        "biomarker_kegg": DATA_ROOT / "SCP259" / "enrichment_results" / "Biomarkers" / "SCP259_inflamed_vs_healthy_crypt_prolif_epi_Gene_Biomarkers_KEGG_results.csv",
        "deseq_go": DATA_ROOT / "SCP259" / "enrichment_results" / "DESeq2_only_genes" / "SCP259_inflamed_vs_healthy_crypt_prolif_epi_Gene_DESeq2_only_genes_GO_results.csv",
        "deseq_kegg": DATA_ROOT / "SCP259" / "enrichment_results" / "DESeq2_only_genes" / "SCP259_inflamed_vs_healthy_crypt_prolif_epi_Gene_DESeq2_only_genes_KEGG_results.csv",
        "kind": "uc",
    },
    "GSE159115": {
        "title": "GSE159115",
        "out_dir": OUT_ROOT / "GSE159115",
        "biomarker_csv": DATA_ROOT / "GSE159115" / "Biomarkers_GSE159115_ccrcc_tumor_vs_ptb_ptc_normal.csv",
        "biomarker_go": DATA_ROOT / "GSE159115" / "enrichment_results" / "Biomarkers" / "GSE159115_ccrcc_tumor_vs_ptb_ptc_normal_Gene_Biomarkers_GO_results.csv",
        "biomarker_kegg": DATA_ROOT / "GSE159115" / "enrichment_results" / "Biomarkers" / "GSE159115_ccrcc_tumor_vs_ptb_ptc_normal_Gene_Biomarkers_KEGG_results.csv",
        "deseq_go": DATA_ROOT / "GSE159115" / "enrichment_results" / "DESeq2_only_genes" / "GSE159115_ccrcc_tumor_vs_ptb_ptc_normal_Gene_DESeq2_only_genes_GO_results.csv",
        "deseq_kegg": DATA_ROOT / "GSE159115" / "enrichment_results" / "DESeq2_only_genes" / "GSE159115_ccrcc_tumor_vs_ptb_ptc_normal_Gene_DESeq2_only_genes_KEGG_results.csv",
        "kind": "ccrcc",
    },
    "GSE138852": {
        "title": "GSE138852",
        "out_dir": OUT_ROOT / "GSE138852",
        "biomarker_csv": DATA_ROOT / "GSE138852" / "Biomarkers.csv",
        "biomarker_go": DATA_ROOT / "GSE138852" / "enrichment_results" / "Biomarkers" / "GSE138852_Gene_Biomarkers_GO_results.csv",
        "biomarker_kegg": None,
        "deseq_go": DATA_ROOT / "GSE138852" / "enrichment_results" / "DESeq2_only_genes" / "GSE138852_Gene_DESeq2_only_genes_GO_results.csv",
        "deseq_kegg": DATA_ROOT / "GSE138852" / "enrichment_results" / "DESeq2_only_genes" / "GSE138852_Gene_DESeq2_only_genes_KEGG_results.csv",
        "kind": "ad",
    },
}


def ensure_exists(path: Path | None) -> None:
    if path is not None and not path.exists():
        raise FileNotFoundError(path)


def wrap(s: str, width: int = 24) -> str:
    return "\n".join(textwrap.wrap(str(s), width=width, break_long_words=False))


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def add_panel_label(ax, label: str) -> None:
    ax.text(-0.08, 1.03, label, transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")


def category_matrix(kind: str, biomarkers: pd.DataFrame) -> pd.DataFrame:
    if kind == "uc":
        category_map = {
            "Secretory": {"TFF1", "S100P", "PLA2G2A", "PHGR1", "PI3", "AGR2", "REG4", "CEACAM5", "FAM3C", "SPINK5"},
            "Immune": {"LYZ", "PDIA3", "HLA-DMA"},
            "ECM": {"BSG", "TIMP1", "CEACAM5", "TSPAN13"},
            "Sulfur": {"TST", "MPST"},
            "Mito": {"COX5A", "ATP5G3"},
            "ER stress": {"PDIA3", "PDIA4", "HSPB1", "AGR2"},
        }
    elif kind == "ccrcc":
        category_map = {
            "Metabolic": {"PFKL", "SLC25A25"},
            "Immune": {"CD68", "CRACR2B"},
            "Epithelial": {"TM4SF18", "DCDC2"},
        }
    elif kind == "ad":
        category_map = {
            "Synapse": {"NRXN1", "NTRK2", "GPM6A"},
            "Neurite": {"CTNNA2", "NTRK2", "GPM6A"},
            "Glial/metabolic": {"SLC1A2", "GPC5"},
            "Stress": {"SPP1"},
        }
    else:
        category_map = {}
    genes = biomarkers["gene"].tolist()
    df = pd.DataFrame(0, index=genes, columns=list(category_map.keys()), dtype=int)
    for cat, members in category_map.items():
        df.loc[df.index.intersection(members), cat] = 1
    return df


def plot_ace_panel(ax, biomarkers: pd.DataFrame, panel_label: str = "A") -> None:
    add_panel_label(ax, panel_label)
    df = biomarkers.sort_values("ACE", ascending=True).copy()
    y = np.arange(len(df))
    colors = np.where(df["ACE"] >= 0, "#D95F02", "#1F77B4")
    ax.hlines(y, 0, df["ACE"], color="#C7C7C7", linewidth=2)
    ax.scatter(df["ACE"], y, s=78, c=colors, edgecolor="white", linewidth=0.8, zorder=3)
    ax.axvline(0, color="#333333", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(df["gene"])
    ax.set_xlabel("ACE")
    ax.set_ylabel("Biomarker gene")
    ax.set_title("Biomarker effect sizes", loc="left", pad=4, fontweight="bold")
    xpad = max(0.08, abs(df["ACE"]).max() * 0.12)
    ax.set_xlim(df["ACE"].min() - xpad, df["ACE"].max() + xpad)
    ax.grid(axis="x", color="#E5E5E5", linewidth=0.8)


def plot_category_panel(ax, kind: str, biomarkers: pd.DataFrame, panel_label: str = "B") -> None:
    add_panel_label(ax, panel_label)
    mat = category_matrix(kind, biomarkers)
    if mat.shape[1] == 0:
        ax.set_axis_off()
        ax.text(0.5, 0.5, "No category annotations", ha="center", va="center")
        return
    cmap = matplotlib.colors.ListedColormap(["#F5F5F5", "#2E8B57"])
    ax.imshow(mat.values, aspect="auto", cmap=cmap, interpolation="nearest")
    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=30, ha="right", rotation_mode="anchor", fontsize=8)
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=8)
    ax.set_title("Functional category annotation", loc="left", pad=4, fontweight="bold")
    ax.tick_params(axis="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlim(-0.5, mat.shape[1] - 0.5)
    ax.set_ylim(mat.shape[0] - 0.5, -0.5)


def prepare_dotplot_df(path: Path | None, top_n: int, kind: str) -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path)
    pcol = "p.adjust" if "p.adjust" in df.columns else "pvalue"
    df = df.copy()
    df["p"] = pd.to_numeric(df[pcol], errors="coerce")
    df["Count"] = pd.to_numeric(df["Count"], errors="coerce")
    df["GeneRatio_num"] = df["GeneRatio"].astype(str).str.split("/").apply(
        lambda x: float(x[0]) / float(x[1]) if len(x) == 2 and float(x[1]) != 0 else np.nan
    )
    df = df.dropna(subset=["p", "Description", "Count", "GeneRatio_num"])
    df = df.sort_values("p").head(top_n).copy()
    df["label"] = df["Description"].map(lambda x: wrap(x, 28))
    df["logp"] = -np.log10(df["p"].clip(lower=1e-300))
    if kind == "KEGG" and len(df) > 5:
        df = df.head(5)
    return df


def plot_dotplot(ax, df: pd.DataFrame, title: str, panel_label: str) -> None:
    add_panel_label(ax, panel_label)
    if df.empty:
        ax.set_axis_off()
        ax.text(0.5, 0.5, "No significant terms", ha="center", va="center", fontsize=11)
        ax.set_title(title, loc="left", pad=4, fontweight="bold")
        return
    y = np.arange(len(df))
    denom = max(1, df["Count"].max() - df["Count"].min())
    sizes = 40 + (df["Count"] - df["Count"].min()) / denom * 220
    sc = ax.scatter(
        df["GeneRatio_num"],
        y,
        s=sizes,
        c=df["logp"],
        cmap="Reds",
        edgecolor="#333333",
        linewidth=0.4,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(df["label"])
    ax.invert_yaxis()
    ax.set_title(title, loc="left", pad=4, fontweight="bold")
    ax.set_xlabel("GeneRatio")
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.8)
    cb = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("-log10(p.adjust)", fontsize=9)


def save_figure(fig, out_dir: Path, stem: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.png", bbox_inches="tight")


def build_panel_figure(dataset_name: str, cfg: dict) -> None:
    biomarkers = pd.read_csv(cfg["biomarker_csv"]).sort_values("ACE", key=lambda s: s.abs(), ascending=False).reset_index(drop=True)
    fig, axes = plt.subplots(1, 2, figsize=(15, 10), gridspec_kw={"width_ratios": [1.15, 1.0]})
    plot_ace_panel(axes[0], biomarkers, "A")
    plot_category_panel(axes[1], cfg["kind"], biomarkers, "B")
    fig.suptitle(f"{cfg['title']} biomarker panel", y=0.98, fontsize=15, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_figure(fig, cfg["out_dir"], f"{dataset_name.lower()}_figure_1_biomarker_panel")
    plt.close(fig)


def build_enrichment_figure(dataset_name: str, cfg: dict, mode: str) -> None:
    if mode == "biomarker":
        go_df = prepare_dotplot_df(cfg["biomarker_go"], top_n=8, kind="GO")
        kegg_df = prepare_dotplot_df(cfg["biomarker_kegg"], top_n=5, kind="KEGG")
        title = f"{cfg['title']} biomarker enrichment"
        stem = f"{dataset_name.lower()}_figure_2_biomarker_enrichment"
    else:
        go_df = prepare_dotplot_df(cfg["deseq_go"], top_n=8, kind="GO")
        kegg_df = prepare_dotplot_df(cfg["deseq_kegg"], top_n=5, kind="KEGG")
        title = f"{cfg['title']} DESeq2-only enrichment"
        stem = f"{dataset_name.lower()}_figure_3_deseq_only_enrichment"

    fig, axes = plt.subplots(1, 2, figsize=(15, 7), gridspec_kw={"width_ratios": [1.3, 0.9]})
    plot_dotplot(axes[0], go_df, "GO enrichment", "A")
    plot_dotplot(axes[1], kegg_df, "KEGG enrichment", "B")
    fig.suptitle(title, y=0.98, fontsize=15, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_figure(fig, cfg["out_dir"], stem)
    plt.close(fig)


def write_blueprint(dataset_name: str, cfg: dict) -> None:
    out = cfg["out_dir"] / f"{dataset_name.lower()}_figure_blueprint.md"
    lines = [
        f"# {cfg['title']} Figure Blueprint",
        "",
        "## Figure 1",
        "- Ranked biomarker ACE plot plus functional category annotation.",
        "- Main message: the final biomarker panel is compact and biologically structured.",
        "",
        "## Figure 2",
        "- Biomarker GO and KEGG enrichment.",
        "- Main message: CSCN biomarker functions capture the dataset-specific mechanistic axis.",
        "",
        "## Figure 3",
        "- DESeq2-only GO and KEGG enrichment.",
        "- Main message: DESeq2-only genes provide the broader differential-expression background.",
    ]
    out.write_text("\n".join(lines) + "\n")


def main():
    setup_style()
    for dataset_name, cfg in DATASETS.items():
        for key in ["biomarker_csv", "biomarker_go", "deseq_go", "deseq_kegg"]:
            ensure_exists(cfg[key])
        build_panel_figure(dataset_name, cfg)
        build_enrichment_figure(dataset_name, cfg, "biomarker")
        build_enrichment_figure(dataset_name, cfg, "deseq")
        write_blueprint(dataset_name, cfg)
        print(f"Saved figures to {cfg['out_dir']}")


if __name__ == "__main__":
    main()
