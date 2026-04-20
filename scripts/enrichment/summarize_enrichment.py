import argparse
import csv
import math
from pathlib import Path
from collections import defaultdict, Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"


def safe_float(val):
    try:
        return float(val)
    except Exception:
        return None


def read_enrichment(file_path):
    if not file_path.exists():
        return []
    rows = []
    with file_path.open(newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            # clusterProfiler uses p.adjust
            p = r.get("p.adjust") or r.get("pvalue") or r.get("pvalue")
            p_val = safe_float(p)
            if p_val is None:
                continue
            rows.append(
                {
                    "Description": r.get("Description", ""),
                    "p": p_val,
                    "ONTOLOGY": r.get("ONTOLOGY", ""),
                }
            )
    return rows


def top_terms(rows, n=5):
    rows = [r for r in rows if r.get("Description")]
    rows.sort(key=lambda r: r["p"])
    return rows[:n]


def get_groups_by_specific_edges(summary_path, top_n=5):
    groups = []
    if not summary_path.exists():
        return groups
    rows = []
    with summary_path.open(newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            try:
                cnt = int(float(r.get("specific_edges_filtered", 0)))
            except Exception:
                cnt = 0
            group = f"Day{r['time']}_{r['cell_type']}"
            rows.append((cnt, group))
    rows.sort(reverse=True)
    for cnt, group in rows[:top_n]:
        groups.append((group, cnt))
    return groups


def group_to_display(group):
    # Day54_cortical_interneuron -> Day54 cortical interneuron
    return group.replace("_", " ")


def group_to_safe(group):
    return (
        group.replace(" ", "_")
        .replace("/", "_")
        .replace("\\", "_")
    )


def build_term_matrix(groups, role, kind, method_dir, top_terms_n=10):
    per_group_terms = {}
    term_counts = Counter()

    for group, _cnt in groups:
        group_safe = group_to_safe(group)
        file_name = f"{method_dir.name}_{group_safe}_{role}_{kind}_results.csv"
        file_path = (
            method_dir
            / "enrichment_results"
            / f"{group_safe}__{role}"
            / file_name
        )
        rows = read_enrichment(file_path)
        best = {}
        for r in rows:
            term = r["Description"]
            p = r["p"]
            if term not in best or p < best[term]:
                best[term] = p
        per_group_terms[group] = best
        for term in best.keys():
            term_counts[term] += 1

    # pick top terms by frequency then best p
    if not term_counts:
        return [], [], None
    sorted_terms = sorted(
        term_counts.items(), key=lambda x: (-x[1], x[0])
    )
    selected_terms = [t for t, _ in sorted_terms[:top_terms_n]]

    matrix = []
    for group, _cnt in groups:
        row = []
        for term in selected_terms:
            p = per_group_terms.get(group, {}).get(term)
            if p is None or p <= 0:
                row.append(0.0)
            else:
                row.append(-math.log10(p))
        matrix.append(row)

    return selected_terms, [g for g, _ in groups], np.array(matrix)


def draw_heatmap(matrix, row_labels, col_labels, title, out_path):
    if matrix is None or matrix.size == 0:
        return
    fig, ax = plt.subplots(figsize=(max(8, len(col_labels) * 0.6), max(6, len(row_labels) * 0.5)))
    im = ax.imshow(matrix, cmap="Reds")
    ax.set_xticks(range(len(col_labels)))
    ax.set_yticks(range(len(row_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels([group_to_display(r) for r in row_labels], fontsize=8)
    ax.set_title(title, fontsize=12, fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def draw_bar_counts(counts, title, out_path):
    # counts: list of (group, sender_count, receiver_count)
    if not counts:
        return
    labels = [group_to_display(g) for g, _, _ in counts]
    sender = [s for _, s, _ in counts]
    receiver = [r for _, _, r in counts]

    x = np.arange(len(labels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 0.6), 6))
    ax.bar(x - width / 2, sender, width, label="sender")
    ax.bar(x + width / 2, receiver, width, label="receiver")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Significant terms")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-set", default="E-GEOD-93593")
    parser.add_argument("--methods", default="kTotal,kWithin")
    parser.add_argument("--top-groups", type=int, default=5)
    parser.add_argument("--top-terms", type=int, default=10)
    args = parser.parse_args()

    data_set = args.data_set
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]

    for method in methods:
        method_dir = DATA_DIR / data_set / "visualizations" / method
        if not method_dir.exists():
            print(f"Missing method dir: {method_dir}")
            continue

        summary_path = method_dir / "specific_edges_filtered" / "specific_edges_summary.csv"
        groups = get_groups_by_specific_edges(summary_path, top_n=args.top_groups)
        if not groups:
            print(f"No groups found for {method}")
            continue

        out_dir = method_dir / "enrichment_summary"
        out_dir.mkdir(parents=True, exist_ok=True)

        # explanation table
        table_path = out_dir / "explanation_table.csv"
        with table_path.open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "method",
                "group",
                "role",
                "go_top_terms",
                "kegg_top_terms",
            ])
            for group, _cnt in groups:
                for role in ["sender", "receiver"]:
                    group_safe = group_to_safe(group)
                    go_file = (
                        method_dir
                        / "enrichment_results"
                        / f"{group_safe}__{role}"
                        / f"{method}_{group_safe}_{role}_GO_results.csv"
                    )
                    kegg_file = (
                        method_dir
                        / "enrichment_results"
                        / f"{group_safe}__{role}"
                        / f"{method}_{group_safe}_{role}_KEGG_results.csv"
                    )
                    go_terms = top_terms(read_enrichment(go_file), n=5)
                    kegg_terms = top_terms(read_enrichment(kegg_file), n=5)

                    go_str = "; ".join([f"{t['Description']} (p={t['p']:.2e})" for t in go_terms])
                    kegg_str = "; ".join([f"{t['Description']} (p={t['p']:.2e})" for t in kegg_terms])

                    writer.writerow([method, group, role, go_str, kegg_str])

        # heatmaps
        for kind in ["GO", "KEGG"]:
            for role in ["sender", "receiver"]:
                terms, group_labels, matrix = build_term_matrix(
                    groups, role, kind, method_dir, top_terms_n=args.top_terms
                )
                if matrix is None or matrix.size == 0:
                    continue
                heat_path = out_dir / f"{kind}_{role}_heatmap.png"
                title = f"{kind} ({role}) - top groups"
                draw_heatmap(matrix, group_labels, terms, title, heat_path)

        # sender/receiver bar charts
        for kind in ["GO", "KEGG"]:
            counts = []
            for group, _cnt in groups:
                group_safe = group_to_safe(group)
                sender_file = (
                    method_dir
                    / "enrichment_results"
                    / f"{group_safe}__sender"
                    / f"{method}_{group_safe}_sender_{kind}_results.csv"
                )
                receiver_file = (
                    method_dir
                    / "enrichment_results"
                    / f"{group_safe}__receiver"
                    / f"{method}_{group_safe}_receiver_{kind}_results.csv"
                )
                sender_terms = read_enrichment(sender_file)
                receiver_terms = read_enrichment(receiver_file)
                counts.append((group, len(sender_terms), len(receiver_terms)))

            bar_path = out_dir / f"{kind}_sender_receiver_counts.png"
            title = f"{kind} sender vs receiver (top groups)"
            draw_bar_counts(counts, title, bar_path)

        print(f"Saved enrichment summary for {method} -> {out_dir}")


if __name__ == "__main__":
    main()
