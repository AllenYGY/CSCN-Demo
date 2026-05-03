from __future__ import annotations

import argparse
import os
import textwrap
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp") / "matplotlib-cscn"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = REPO_ROOT / "data"
OUT_DIR = REPO_ROOT / "results" / "paper_figures" / "cross_dataset"

JOURNAL_CMAP = LinearSegmentedColormap.from_list(
    "journal_blue",
    ["#F7F6F2", "#DCE5EA", "#A8BEC9", "#6E8FA5", "#2F5068"],
)
TEXT_COLOR = "#1F2933"
SUBTLE_LINE = "#C9CED3"
POINT_EDGE = "#3A4A57"

DATASETS = [
    {
        "key":
        "BreastTumer",
        "display":
        "BreastTumer",
        "disease":
        "Breast tumor",
        "biomarker_count":
        62,
        "go_csv":
        DATA_ROOT / "BreastTumer" / "enrichment_results" / "Biomarkers" /
        "BreastTumer_Gene_Biomarkers_GO_results.csv",
        "kegg_csv":
        DATA_ROOT / "BreastTumer" / "enrichment_results" / "Biomarkers" /
        "BreastTumer_Gene_Biomarkers_KEGG_results.csv",
        "deseq_go_csv":
        DATA_ROOT / "BreastTumer" / "enrichment_results" /
        "DESeq2_only_genes" /
        "BreastTumer_Gene_DESeq2_only_genes_GO_results.csv",
        "deseq_kegg_csv":
        DATA_ROOT / "BreastTumer" / "enrichment_results" /
        "DESeq2_only_genes" /
        "BreastTumer_Gene_DESeq2_only_genes_KEGG_results.csv",
    },
    {
        "key":
        "SCP259",
        "display":
        "SCP259",
        "disease":
        "Ulcerative colitis",
        "biomarker_count":
        24,
        "go_csv":
        DATA_ROOT / "SCP259" / "enrichment_results" / "Biomarkers" /
        "SCP259_inflamed_vs_healthy_crypt_prolif_epi_Gene_Biomarkers_GO_results.csv",
        "kegg_csv":
        DATA_ROOT / "SCP259" / "enrichment_results" / "Biomarkers" /
        "SCP259_inflamed_vs_healthy_crypt_prolif_epi_Gene_Biomarkers_KEGG_results.csv",
        "deseq_go_csv":
        DATA_ROOT / "SCP259" / "enrichment_results" / "DESeq2_only_genes" /
        "SCP259_inflamed_vs_healthy_crypt_prolif_epi_Gene_DESeq2_only_genes_GO_results.csv",
        "deseq_kegg_csv":
        DATA_ROOT / "SCP259" / "enrichment_results" / "DESeq2_only_genes" /
        "SCP259_inflamed_vs_healthy_crypt_prolif_epi_Gene_DESeq2_only_genes_KEGG_results.csv",
    },
    {
        "key":
        "GSE159115",
        "display":
        "GSE159115",
        "disease":
        "ccRCC",
        "biomarker_count":
        6,
        "go_csv":
        DATA_ROOT / "GSE159115" / "enrichment_results" / "Biomarkers" /
        "GSE159115_ccrcc_tumor_vs_ptb_ptc_normal_Gene_Biomarkers_GO_results.csv",
        "kegg_csv":
        DATA_ROOT / "GSE159115" / "enrichment_results" / "Biomarkers" /
        "GSE159115_ccrcc_tumor_vs_ptb_ptc_normal_Gene_Biomarkers_KEGG_results.csv",
        "deseq_go_csv":
        DATA_ROOT / "GSE159115" / "enrichment_results" / "DESeq2_only_genes" /
        "GSE159115_ccrcc_tumor_vs_ptb_ptc_normal_Gene_DESeq2_only_genes_GO_results.csv",
        "deseq_kegg_csv":
        DATA_ROOT / "GSE159115" / "enrichment_results" / "DESeq2_only_genes" /
        "GSE159115_ccrcc_tumor_vs_ptb_ptc_normal_Gene_DESeq2_only_genes_KEGG_results.csv",
    },
    {
        "key":
        "GSE138852",
        "display":
        "GSE138852",
        "disease":
        "Alzheimer's disease",
        "biomarker_count":
        7,
        "go_csv":
        DATA_ROOT / "GSE138852" / "enrichment_results" / "Biomarkers" /
        "GSE138852_Gene_Biomarkers_GO_results.csv",
        "kegg_csv":
        None,
        "deseq_go_csv":
        DATA_ROOT / "GSE138852" / "enrichment_results" / "DESeq2_only_genes" /
        "GSE138852_Gene_DESeq2_only_genes_GO_results.csv",
        "deseq_kegg_csv":
        DATA_ROOT / "GSE138852" / "enrichment_results" / "DESeq2_only_genes" /
        "GSE138852_Gene_DESeq2_only_genes_KEGG_results.csv",
    },
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=
        "Generate stitched cross-dataset biomarker enrichment figures.")
    parser.add_argument(
        "--exclude",
        action="append",
        default=[],
        help="Dataset display name to exclude. Can be passed multiple times.",
    )
    parser.add_argument(
        "--prefix",
        default="cross_dataset",
        help="Output filename prefix. Default: cross_dataset",
    )
    parser.add_argument(
        "--source",
        choices=["biomarker", "deseq_only"],
        default="biomarker",
        help="Which enrichment source to visualize. Default: biomarker",
    )
    return parser.parse_args()


def ensure_exists(path: Path | None) -> None:
    if path is not None and not path.exists():
        raise FileNotFoundError(path)


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.titlesize": 20,
        "axes.labelsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 13,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "text.color": TEXT_COLOR,
        "axes.labelcolor": TEXT_COLOR,
        "axes.titlecolor": TEXT_COLOR,
        "xtick.color": TEXT_COLOR,
        "ytick.color": TEXT_COLOR,
    })


def wrap(text: str, width: int = 24) -> str:
    return "\n".join(
        textwrap.wrap(str(text), width=width, break_long_words=False))


def add_dataset_title(ax, title: str) -> None:
    return


def dataset_title(ds: dict) -> str:
    return ds["display"]


def select_datasets(excluded: list[str]) -> list[dict]:
    excluded_set = {x.strip() for x in excluded if x and x.strip()}
    selected = [ds for ds in DATASETS if ds["display"] not in excluded_set]
    if not selected:
        raise ValueError("No datasets remain after exclusion.")
    return selected


def prepare_df(path: Path | None, kind: str, top_n: int = 8) -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path).copy()
    pcol = "p.adjust" if "p.adjust" in df.columns else "pvalue"
    df["p"] = pd.to_numeric(df[pcol], errors="coerce")
    df["Count"] = pd.to_numeric(df["Count"], errors="coerce")
    df = df.dropna(subset=["p", "Description", "Count"]).copy()
    df = df.sort_values("p").head(6 if kind == "GO" else min(top_n, 5)).copy()
    df["label"] = df["Description"].map(
        lambda x: wrap(x, 34 if kind == "GO" else 28))
    df["logp"] = -np.log10(df["p"].clip(lower=1e-300))
    return df


def draw_panel(
        ax,
        df: pd.DataFrame,
        title: str,
        panel_label: str,
        norm: Normalize,
        *,
        rotate_y_labels: bool = False,
        dot_x: float = -0.045,
        x_limits: tuple[float, float] = (-0.20, 0.30),
) -> None:
    add_dataset_title(ax, title)
    if df.empty:
        ax.set_axis_off()
        ax.text(0.5,
                0.5,
                "No significant terms",
                ha="center",
                va="center",
                fontsize=12,
                color="#6B7280")
        return

    y = np.arange(len(df))
    denom = max(1, df["logp"].max() - df["logp"].min())
    sizes = 110 + (df["logp"] - df["logp"].min()) / denom * 150
    ax.scatter(
        np.full(len(df), dot_x),
        y,
        s=sizes,
        c=df["logp"],
        cmap=JOURNAL_CMAP,
        norm=norm,
        edgecolor=POINT_EDGE,
        linewidth=0.5,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(
        df["label"],
        fontsize=20,
        rotation=15 if rotate_y_labels else 0,
        ha="right" if rotate_y_labels else "right",
        rotation_mode="anchor" if rotate_y_labels else None,
    )
    ax.invert_yaxis()
    ax.set_xlim(*x_limits)
    ax.set_xticks([])
    ax.set_xlabel("")
    ax.grid(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)


def make_single_panel_figures(kind: str, out_stem: str, datasets: list[dict],
                              source: str) -> None:
    if source == "biomarker":
        key = "go_csv" if kind == "GO" else "kegg_csv"
    else:
        key = "deseq_go_csv" if kind == "GO" else "deseq_kegg_csv"

    rotate_y_labels = (kind == "GO")
    if source == "biomarker" and kind == "GO":
        dot_x = -0.06
        x_limits = (-0.18, 0.18)
    elif rotate_y_labels:
        dot_x = -0.15
        x_limits = (-0.22, 0.18)
    else:
        dot_x = -0.07
        x_limits = (-0.24, 0.22)

    for ds in datasets:
        df = prepare_df(ds[key], kind)
        fig, ax = plt.subplots(1, 1, figsize=(5.2, 7.2))
        fig.patch.set_facecolor("white")
        if df.empty:
            norm = Normalize(vmin=0, vmax=1)
        else:
            vmin = float(df["logp"].min())
            vmax = float(df["logp"].max())
            if np.isclose(vmin, vmax):
                vmin = 0.0
                vmax = max(1.0, vmax)
            norm = Normalize(vmin=vmin, vmax=vmax)
        draw_panel(
            ax,
            df,
            dataset_title(ds),
            "A",
            norm,
            rotate_y_labels=rotate_y_labels,
            dot_x=dot_x,
            x_limits=x_limits,
        )
        ax.set_facecolor("white")
        fig.subplots_adjust(left=0.18, right=0.3, bottom=0.10, top=0.90)
        fig.savefig(OUT_DIR / f"{out_stem}_{ds['display']}.pdf",
                    bbox_inches="tight")
        fig.savefig(OUT_DIR / f"{out_stem}_{ds['display']}.png",
                    bbox_inches="tight")
        plt.close(fig)


def make_figure(kind: str, out_stem: str, title: str, datasets: list[dict],
                source: str) -> None:
    if source == "biomarker":
        key = "go_csv" if kind == "GO" else "kegg_csv"
    else:
        key = "deseq_go_csv" if kind == "GO" else "deseq_kegg_csv"
    data_frames = [prepare_df(ds[key], kind) for ds in datasets]

    fig, axes = plt.subplots(1,
                             len(datasets),
                             figsize=(4.6 * len(datasets), 7.2),
                             gridspec_kw={"wspace": 0.32})
    if len(datasets) == 1:
        axes = [axes]
    fig.patch.set_facecolor("white")
    rotate_y_labels = (source == "deseq_only" and kind == "GO")
    dot_x = -0.075 if rotate_y_labels else -0.045
    x_limits = (-0.24, 0.30) if rotate_y_labels else (-0.20, 0.30)
    for idx, (ax, ds, df) in enumerate(zip(axes, datasets, data_frames)):
        if df.empty:
            norm = Normalize(vmin=0, vmax=1)
        else:
            vmin = float(df["logp"].min())
            vmax = float(df["logp"].max())
            if np.isclose(vmin, vmax):
                vmin = 0.0
                vmax = max(1.0, vmax)
            norm = Normalize(vmin=vmin, vmax=vmax)
        draw_panel(
            ax,
            df,
            dataset_title(ds),
            chr(ord("A") + idx),
            norm,
            rotate_y_labels=rotate_y_labels,
            dot_x=dot_x,
            x_limits=x_limits,
        )
        ax.set_facecolor("white")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(left=0.075,
                        right=0.985,
                        bottom=0.09,
                        top=0.94,
                        wspace=0.34)

    fig.savefig(OUT_DIR / f"{out_stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{out_stem}.png", bbox_inches="tight")
    plt.close(fig)


def write_term_tables(datasets: list[dict], prefix: str, source: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    go_tables = []
    kegg_tables = []
    lines = [
        "# Cross-dataset stitched term summary",
        "",
        "Each dataset keeps its own top enriched terms.",
        "The figures do not force a shared biological-process axis.",
        f"Source: {source}",
        "",
    ]
    for ds in datasets:
        lines.append(f"## {ds['display']}")
        go_key = "go_csv" if source == "biomarker" else "deseq_go_csv"
        kegg_key = "kegg_csv" if source == "biomarker" else "deseq_kegg_csv"
        go_df = prepare_df(ds[go_key], "GO")
        go_df.insert(0, "dataset", ds["display"])
        go_tables.append(go_df[[
            "dataset", "Description", "p", "p.adjust" if "p.adjust" in
            go_df.columns else "p", "Count", "geneID"
        ] if not go_df.empty else ["dataset"]])
        lines.append("### GO")
        if go_df.empty:
            lines.append("- No significant GO terms")
        else:
            for _, row in go_df.iterrows():
                lines.append(
                    f"- {row['Description']} (padj={row['p']:.2e}, Count={int(row['Count'])})"
                )

        kegg_df = prepare_df(ds[kegg_key], "KEGG")
        if not kegg_df.empty:
            kegg_df.insert(0, "dataset", ds["display"])
            kegg_tables.append(kegg_df[[
                "dataset", "Description", "p",
                "p.adjust" if "p.adjust" in kegg_df.columns else "p", "Count",
                "geneID"
            ]])
        lines.append("### KEGG")
        if kegg_df.empty:
            lines.append("- No significant KEGG terms")
        else:
            for _, row in kegg_df.iterrows():
                lines.append(
                    f"- {row['Description']} (padj={row['p']:.2e}, Count={int(row['Count'])})"
                )
        lines.append("")

    if go_tables:
        pd.concat(go_tables, ignore_index=True).to_csv(
            OUT_DIR / f"{prefix}_biomarker_go_terms.csv", index=False)
    if kegg_tables:
        pd.concat(kegg_tables, ignore_index=True).to_csv(
            OUT_DIR / f"{prefix}_biomarker_kegg_terms.csv", index=False)
    (OUT_DIR /
     f"{prefix}_biomarker_process_mapping.md").write_text("\n".join(lines) +
                                                          "\n")


def main() -> None:
    args = parse_args()
    datasets = select_datasets(args.exclude)
    setup_style()
    for ds in datasets:
        if args.source == "biomarker":
            ensure_exists(ds["go_csv"])
            ensure_exists(ds["kegg_csv"])
        else:
            ensure_exists(ds["deseq_go_csv"])
            ensure_exists(ds["deseq_kegg_csv"])
    source_label = "biomarker" if args.source == "biomarker" else "DESeq2-only"
    make_figure("GO", f"{args.prefix}_go_summary",
                f"Cross-dataset {source_label} GO term panels", datasets,
                args.source)
    make_figure("KEGG", f"{args.prefix}_kegg_summary",
                f"Cross-dataset {source_label} KEGG pathway panels", datasets,
                args.source)
    make_single_panel_figures("GO", f"{args.prefix}_go_single", datasets,
                              args.source)
    write_term_tables(datasets, args.prefix, args.source)
    print(f"Saved cross-dataset figures to {OUT_DIR}")


if __name__ == "__main__":
    main()
