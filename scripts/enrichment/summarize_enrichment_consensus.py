import argparse
import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"


def safe_float(val):
    try:
        return float(val)
    except Exception:
        return None


def read_top_terms(path, n=5):
    if not path.exists():
        return []
    rows = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            p = r.get("p.adjust") or r.get("pvalue") or r.get("pvalue")
            p_val = safe_float(p)
            if p_val is None:
                continue
            desc = r.get("Description", "")
            if not desc:
                continue
            rows.append(
                {
                    "Description": desc,
                    "p": p_val,
                    "GeneRatio": r.get("GeneRatio", ""),
                    "Count": r.get("Count", ""),
                }
            )
    rows.sort(key=lambda r: r["p"])
    return rows[:n]


def parse_group(group):
    # group like Day54_cortical_interneuron
    if not group.startswith("Day"):
        return None, None
    try:
        day_part, cell_type = group.split("_", 1)
    except ValueError:
        return None, None
    time_val = day_part.replace("Day", "")
    return time_val, cell_type


def write_summary(out_path, rows):
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "method",
                "time",
                "cell_type",
                "top_terms",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-set", default="E-GEOD-93593")
    parser.add_argument("--methods", default="kTotal,kWithin")
    parser.add_argument("--top-terms", type=int, default=5)
    args = parser.parse_args()

    data_set = args.data_set
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]

    for method in methods:
        base = DATA_DIR / data_set / "visualizations" / method / "enrichment_results_consensus"
        if not base.exists():
            print(f"Missing: {base}")
            continue

        go_rows = []
        kegg_rows = []
        for group_dir in sorted(p for p in base.iterdir() if p.is_dir() and p.name.startswith("Day")):
            group = group_dir.name
            time_val, cell_type = parse_group(group)
            if time_val is None:
                continue

            go_path = group_dir / f"{method}_{group}_consensus_GO_results.csv"
            kegg_path = group_dir / f"{method}_{group}_consensus_KEGG_results.csv"

            go_terms = read_top_terms(go_path, n=args.top_terms)
            kegg_terms = read_top_terms(kegg_path, n=args.top_terms)

            go_str = "; ".join([f"{t['Description']} (p={t['p']:.2e})" for t in go_terms])
            kegg_str = "; ".join([f"{t['Description']} (p={t['p']:.2e})" for t in kegg_terms])

            go_rows.append(
                {
                    "method": method,
                    "time": time_val,
                    "cell_type": cell_type.replace("_", " "),
                    "top_terms": go_str,
                }
            )
            kegg_rows.append(
                {
                    "method": method,
                    "time": time_val,
                    "cell_type": cell_type.replace("_", " "),
                    "top_terms": kegg_str,
                }
            )

        if go_rows:
            write_summary(base / "consensus_GO_summary.csv", go_rows)
        if kegg_rows:
            write_summary(base / "consensus_KEGG_summary.csv", kegg_rows)

        print(f"Saved summaries for {method} -> {base}")


if __name__ == "__main__":
    main()
