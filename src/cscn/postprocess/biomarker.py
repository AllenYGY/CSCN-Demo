from __future__ import annotations

import json
from pathlib import Path

from biomarker.datasets import build_expression_df
from biomarker.graph_utils import identify_biomarkers_from_group_graphs

from ..aggregate import build_group_graph, load_group_dags, load_node_map, map_group_dags_to_genes


def run_biomarker(config, layout, prepared_run) -> Path:
    case_group = config.biomarker.case_group
    control_group = config.biomarker.control_group
    if not case_group or not control_group:
        raise ValueError(
            "Biomarker mode requires `biomarker.case_group` and `biomarker.control_group`."
        )
    if case_group == control_group:
        raise ValueError("Biomarker case_group and control_group must be different.")
    if case_group not in prepared_run.groups or control_group not in prepared_run.groups:
        raise ValueError(
            f"Biomarker groups must exist in the prepared run. Available groups: {sorted(prepared_run.groups)}"
        )

    group_to_label = {
        control_group: 0,
        case_group: 1,
    }
    matrices = {
        group_key: prepared_run.groups[group_key].matrix
        for group_key in (control_group, case_group)
    }
    expression_df = build_expression_df(
        matrices=matrices,
        gene_names=prepared_run.gene_names,
        group_to_label=group_to_label,
    )

    node_map = load_node_map(layout.run_dir)
    group_graphs = {}
    for group_key in (control_group, case_group):
        dags = load_group_dags(layout.dag_group_dir(group_key))
        mapped_dags = map_group_dags_to_genes(dags, node_map)
        group_graphs[group_key] = build_group_graph(mapped_dags)

    biomarkers_df = identify_biomarkers_from_group_graphs(
        group_graphs=group_graphs,
        expression_df=expression_df,
        gene_names=prepared_run.gene_names,
        outcome=config.biomarker.outcome_column,
        sink_node_name=config.biomarker.outcome_column,
        confounder_method=config.biomarker.confounder_method,
        include_n_confounders=True,
        sort_by_abs_ace=True,
        skip_missing_genes=True,
    )

    layout.biomarker_dir.mkdir(parents=True, exist_ok=True)
    biomarker_path = layout.biomarker_dir / "biomarkers.csv"
    summary_path = layout.biomarker_dir / "causal_summary.json"
    biomarkers_df.to_csv(biomarker_path, index=False)
    summary_path.write_text(
        json.dumps(
            {
                "case_group": case_group,
                "control_group": control_group,
                "outcome_column": config.biomarker.outcome_column,
                "biomarker_count": int(len(biomarkers_df)),
                "gene_count": int(len(prepared_run.gene_names)),
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    return biomarker_path
