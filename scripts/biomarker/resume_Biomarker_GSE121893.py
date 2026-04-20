from __future__ import annotations

import argparse
import os
import pickle
import re
import sys
from functools import lru_cache
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


BASE_DATA_SET = "GSE121893"
DEFAULT_CASE_GROUP = "dHF"
DEFAULT_CONTROL_GROUP = "N"
DEFAULT_REGION = "all"
DEFAULT_CELL_COMPARTMENT = "all"
DEFAULT_MAX_WORKERS = min(4, os.cpu_count() or 1)
DEFAULT_PROGRESS_INTERVAL = 10
DEFAULT_SIGNIFICANCE_LEVEL = 0.05
DEFAULT_USE_BITMAP = True
RESULT_PATTERN = re.compile(r"result_(\d+)\.pkl$")


@lru_cache(maxsize=1)
def load_runtime_modules():
    from biomarker.causal import run_causal_analysis
    from biomarker.cscn import CSCN
    from biomarker.datasets import (
        build_expression_df,
        load_gene_names,
        load_saved_group_graphs,
    )
    from biomarker.graph_utils import (
        get_global_graph,
        identify_biomarkers_from_group_graphs,
        map_node_id_to_gene,
    )

    return {
        "run_causal_analysis": run_causal_analysis,
        "CSCN": CSCN,
        "build_expression_df": build_expression_df,
        "load_gene_names": load_gene_names,
        "load_saved_group_graphs": load_saved_group_graphs,
        "get_global_graph": get_global_graph,
        "identify_biomarkers_from_group_graphs": identify_biomarkers_from_group_graphs,
        "map_node_id_to_gene": map_node_id_to_gene,
    }


def log(message):
    print(f"[{BASE_DATA_SET}] {message}")


def log_stage(title):
    print()
    print(f"=== {title} ===")


def validate_required_file(path: Path, description: str):
    if not path.exists():
        raise FileNotFoundError(f"Missing {description}: {path}")
    log(f"{description}: {path}")


def slugify(value: str):
    return "".join(char.lower() if char.isalnum() else "_" for char in value)


def default_run_name(
    case_group: str = DEFAULT_CASE_GROUP,
    control_group: str = DEFAULT_CONTROL_GROUP,
    region: str = DEFAULT_REGION,
    cell_compartment: str = DEFAULT_CELL_COMPARTMENT,
):
    return (
        f"{BASE_DATA_SET}_{slugify(case_group)}_vs_{slugify(control_group)}_"
        f"{region.lower()}_{cell_compartment.lower()}"
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Resume the current GSE121893 dHF biomarker run by only computing "
            "missing DAGs for one group, then finalize the biomarker table."
        )
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=REPO_ROOT / "data" / BASE_DATA_SET,
        help="Dataset directory that already contains prepared matrices and partial DAG outputs.",
    )
    parser.add_argument(
        "--run-name",
        default=default_run_name(),
        help="Existing run name prefix. Default matches the current dHF vs N all/all run.",
    )
    parser.add_argument(
        "--group",
        default=DEFAULT_CASE_GROUP,
        help="Group to resume. Default: dHF.",
    )
    parser.add_argument(
        "--case-group",
        default=DEFAULT_CASE_GROUP,
        help="Case cohort label used during final biomarker aggregation.",
    )
    parser.add_argument(
        "--control-group",
        default=DEFAULT_CONTROL_GROUP,
        help="Control cohort label used during final biomarker aggregation.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=DEFAULT_MAX_WORKERS,
        help="Worker count for missing DAG tasks. Default is conservative for long resumes.",
    )
    parser.add_argument(
        "--progress-interval",
        type=int,
        default=DEFAULT_PROGRESS_INTERVAL,
        help="Emit a progress line every N completed DAGs during the resume step.",
    )
    parser.add_argument(
        "--significance-level",
        type=float,
        default=DEFAULT_SIGNIFICANCE_LEVEL,
        help="Conditional independence significance level passed to CSCN.",
    )
    parser.add_argument(
        "--skip-finalize",
        action="store_true",
        help="Only resume the requested group and stop before biomarker aggregation.",
    )
    return parser.parse_args()


def discover_valid_result_ids(dag_dir: Path, expected_total: int):
    valid_ids = set()
    invalid_paths = []

    if not dag_dir.exists():
        return valid_ids, invalid_paths

    for path in dag_dir.glob("result_*.pkl"):
        match = RESULT_PATTERN.match(path.name)
        if match is None:
            continue
        task_id = int(match.group(1))
        if not 0 <= task_id < expected_total:
            continue

        try:
            with open(path, "rb") as handle:
                pickle.load(handle)
        except Exception:
            invalid_paths.append(path)
            continue

        valid_ids.add(task_id)

    return valid_ids, invalid_paths


def ensure_cscn_object(
    data_dir: Path,
    run_name: str,
    group: str,
    dag_dir: Path,
    significance_level: float,
):
    CSCN = load_runtime_modules()["CSCN"]
    cscn_path = data_dir / f"{run_name}_{group}_cscn"
    if cscn_path.exists():
        return cscn_path

    cscn = CSCN(
        output_dir=str(dag_dir),
        sigmoid_score=0.1,
        significance_level=significance_level,
        max_cond_vars=20,
        use_bitmap=DEFAULT_USE_BITMAP,
        debug=False,
        show_progress=False,
        progress_interval=100,
    )
    CSCN.save_to_file(cscn, cscn_path)
    log(f"created CSCN placeholder object: {cscn_path}")
    return cscn_path


def resume_group_cscn(
    data_dir: Path,
    run_name: str,
    group: str,
    matrix: np.ndarray,
    max_workers: int,
    progress_interval: int,
    significance_level: float,
):
    CSCN = load_runtime_modules()["CSCN"]
    dag_dir = data_dir / "DAG" / run_name / group
    dag_dir.mkdir(parents=True, exist_ok=True)
    cscn_path = data_dir / f"{run_name}_{group}_cscn"

    total = int(matrix.shape[0])
    valid_ids, invalid_paths = discover_valid_result_ids(dag_dir, expected_total=total)
    missing_ids = [task_id for task_id in range(total) if task_id not in valid_ids]

    log_stage(f"Resume CSCN {group}")
    log(f"group: {group}")
    log(f"matrix shape: {matrix.shape}")
    log(f"DAG dir: {dag_dir}")
    log(f"CSCN object path: {cscn_path}")
    log(f"max_workers: {max_workers}")
    log(f"significance level: {significance_level}")
    log(f"use bitmap: {DEFAULT_USE_BITMAP}")
    log(f"valid existing DAGs: {len(valid_ids)}/{total}")
    if invalid_paths:
        log(f"invalid DAG files to overwrite: {len(invalid_paths)}")
        for path in invalid_paths[:5]:
            log(f"invalid DAG example: {path.name}")
    log(f"missing DAGs to compute: {len(missing_ids)}")

    cscn = CSCN(
        output_dir=str(dag_dir),
        sigmoid_score=0.1,
        significance_level=significance_level,
        max_cond_vars=20,
        use_bitmap=DEFAULT_USE_BITMAP,
        debug=False,
        show_progress=False,
        progress_interval=progress_interval,
    )
    cscn.run_core(matrix, usingNMF=False)

    if missing_ids:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(cscn.run_pc_and_save, task_id): task_id
                for task_id in missing_ids
            }
            completed_new = 0
            completed_total = len(valid_ids)
            for future in as_completed(futures):
                task_id = futures[future]
                filename = future.result()
                completed_new += 1
                completed_total += 1
                if progress_interval and (
                    completed_new == 1
                    or completed_total % progress_interval == 0
                    or completed_total == total
                ):
                    print(
                        f"[CSCN-Resume] {group}: saved DAG {completed_total}/{total} -> "
                        f"{os.path.basename(filename)}"
                    )
    else:
        log(f"all {total} DAGs are already present for {group}; nothing to compute")

    CSCN.save_to_file(cscn, cscn_path)
    log(f"saved CSCN object: {cscn_path}")

    final_ids, final_invalid = discover_valid_result_ids(dag_dir, expected_total=total)
    if final_invalid:
        raise RuntimeError(
            f"{group} still has unreadable DAG files after resume: "
            f"{[path.name for path in final_invalid[:5]]}"
        )
    if len(final_ids) != total:
        raise RuntimeError(
            f"{group} DAG resume incomplete: expected {total}, found {len(final_ids)}"
        )
    log(f"{group} DAGs complete: {len(final_ids)}/{total}")
    return dag_dir, cscn_path


def resolve_used_genes_path(output_dir: Path, run_name: str):
    candidates = sorted(output_dir.glob(f"{run_name}_top*_genes_used.csv"))
    if not candidates:
        raise FileNotFoundError(
            f"Could not find a used-genes CSV matching {run_name}_top*_genes_used.csv in {output_dir}"
        )
    if len(candidates) > 1:
        log(
            f"multiple used-genes files found; using the first sorted candidate: {candidates[0].name}"
        )
    return candidates[0]


def finalize_biomarkers(
    data_dir: Path,
    output_dir: Path,
    run_name: str,
    control_group: str,
    case_group: str,
):
    runtime = load_runtime_modules()
    build_expression_df = runtime["build_expression_df"]
    load_gene_names = runtime["load_gene_names"]
    load_saved_group_graphs = runtime["load_saved_group_graphs"]
    CSCN = runtime["CSCN"]
    map_node_id_to_gene = runtime["map_node_id_to_gene"]
    get_global_graph = runtime["get_global_graph"]
    identify_biomarkers_from_group_graphs = runtime["identify_biomarkers_from_group_graphs"]
    run_causal_analysis = runtime["run_causal_analysis"]

    log_stage("Finalize Biomarkers")
    used_genes_path = resolve_used_genes_path(output_dir, run_name)
    validate_required_file(used_genes_path, "used-gene list")

    group_to_label = {
        control_group: 0,
        case_group: 1,
    }
    matrices = {}
    for group in group_to_label:
        matrix_path = output_dir / f"{run_name}_{group}.npy"
        validate_required_file(matrix_path, f"{group} normalized matrix")
        matrices[group] = np.load(matrix_path)
        log(f"{group} matrix shape: {matrices[group].shape}")

        dag_dir = data_dir / "DAG" / run_name / group
        expected_total = int(matrices[group].shape[0])
        valid_ids, invalid_paths = discover_valid_result_ids(dag_dir, expected_total)
        if invalid_paths:
            raise RuntimeError(
                f"{group} has unreadable DAG files: {[path.name for path in invalid_paths[:5]]}"
            )
        if len(valid_ids) != expected_total:
            raise RuntimeError(
                f"{group} DAG count mismatch: expected {expected_total}, found {len(valid_ids)}"
            )
        ensure_cscn_object(data_dir, run_name, group, dag_dir)

    gene_names = load_gene_names(used_genes_path)
    n_genes = int(matrices[control_group].shape[1])
    if len(gene_names) != n_genes:
        raise ValueError(
            f"Gene count mismatch: used-genes file has {len(gene_names)} rows but matrices have {n_genes} columns"
        )

    expression_df = build_expression_df(matrices, gene_names, group_to_label)
    log(f"combined expression dataframe shape: {expression_df.shape}")

    group_graphs = load_saved_group_graphs(
        data_dir=data_dir,
        dataset_name=run_name,
        groups=group_to_label.keys(),
        gene_names=gene_names,
        cscn_cls=CSCN,
        map_node_id_to_gene_fn=map_node_id_to_gene,
        get_global_graph_fn=get_global_graph,
    )
    for group, graph in group_graphs.items():
        log(
            f"group graph {group}: nodes={graph.number_of_nodes()}, "
            f"edges={graph.number_of_edges()}"
        )

    biomarkers_df = identify_biomarkers_from_group_graphs(
        group_graphs=group_graphs,
        expression_df=expression_df,
        gene_names=gene_names,
        outcome="DISEASE",
        sink_node_name="DISEASE",
        confounder_method="classic",
        include_n_confounders=True,
        sort_by_abs_ace=True,
        skip_missing_genes=True,
        run_causal_analysis_fn=run_causal_analysis,
    )
    biomarker_path = data_dir / f"Biomarkers_{run_name}.csv"
    biomarkers_df.to_csv(biomarker_path, index=False)
    log(f"saved biomarkers: {biomarker_path}")
    log(f"biomarker count: {len(biomarkers_df)}")
    if not biomarkers_df.empty:
        log("top 5 biomarkers:")
        print(biomarkers_df.head(5).to_string(index=False))


def main():
    args = parse_args()
    if args.case_group == args.control_group:
        raise ValueError("--case-group and --control-group must be different")
    if args.group not in {args.case_group, args.control_group}:
        raise ValueError("--group must match either --case-group or --control-group")
    if args.max_workers <= 0:
        raise ValueError("--max-workers must be a positive integer")
    if args.progress_interval <= 0:
        raise ValueError("--progress-interval must be a positive integer")

    data_dir = args.data_dir.resolve()
    output_dir = data_dir / "output_deseq"
    matrix_path = output_dir / f"{args.run_name}_{args.group}.npy"

    log_stage("Configuration")
    log(f"dataset dir: {data_dir}")
    log(f"output dir: {output_dir}")
    log(f"run name: {args.run_name}")
    log(f"resume group: {args.group}")
    log(f"case group: {args.case_group}")
    log(f"control group: {args.control_group}")
    log(f"max workers: {args.max_workers}")
    log(f"progress interval: {args.progress_interval}")
    log(f"significance level: {args.significance_level}")
    log(f"skip finalize: {args.skip_finalize}")

    validate_required_file(matrix_path, f"{args.group} normalized matrix")
    matrix = np.load(matrix_path)
    log(f"loaded {args.group} matrix shape: {matrix.shape}")

    resume_group_cscn(
        data_dir=data_dir,
        run_name=args.run_name,
        group=args.group,
        matrix=matrix,
        max_workers=args.max_workers,
        progress_interval=args.progress_interval,
        significance_level=args.significance_level,
    )

    if args.skip_finalize:
        log("skip-finalize enabled; stopping after the resume step")
        return

    finalize_biomarkers(
        data_dir=data_dir,
        output_dir=output_dir,
        run_name=args.run_name,
        control_group=args.control_group,
        case_group=args.case_group,
    )


if __name__ == "__main__":
    main()
