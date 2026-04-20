from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


def safe_component(value: str) -> str:
    normalized = str(value).strip()
    if not normalized:
        raise ValueError("Run name cannot be empty.")
    return normalized.replace("/", "_").replace("\\", "_")


@dataclass(frozen=True)
class RunLayout:
    base_output_dir: Path
    run_name: str
    run_dir: Path
    inputs_dir: Path
    matrices_dir: Path
    dags_dir: Path
    objects_dir: Path
    consensus_dir: Path
    biomarker_dir: Path
    logs_dir: Path
    config_snapshot_path: Path
    group_manifest_path: Path
    genes_path: Path
    cell_metadata_path: Path
    summary_path: Path

    @classmethod
    def from_config(cls, config) -> "RunLayout":
        run_name = safe_component(config.run_name)
        run_dir = config.run.output_dir / run_name
        return cls(
            base_output_dir=config.run.output_dir,
            run_name=run_name,
            run_dir=run_dir,
            inputs_dir=run_dir / "inputs",
            matrices_dir=run_dir / "matrices",
            dags_dir=run_dir / "dags",
            objects_dir=run_dir / "objects",
            consensus_dir=run_dir / "consensus",
            biomarker_dir=run_dir / "biomarker",
            logs_dir=run_dir / "logs",
            config_snapshot_path=run_dir / "config.snapshot.yaml",
            group_manifest_path=run_dir / "inputs" / "groups.csv",
            genes_path=run_dir / "inputs" / "genes.csv",
            cell_metadata_path=run_dir / "inputs" / "cell_metadata.csv",
            summary_path=run_dir / "inputs" / "run_summary.json",
        )

    def ensure_dirs(self) -> None:
        for path in (
            self.run_dir,
            self.inputs_dir,
            self.matrices_dir,
            self.dags_dir,
            self.objects_dir,
            self.consensus_dir,
            self.biomarker_dir,
            self.logs_dir,
        ):
            path.mkdir(parents=True, exist_ok=True)

    def matrix_path(self, group_key: str) -> Path:
        return self.matrices_dir / f"{safe_component(group_key)}.npy"

    def group_cells_path(self, group_key: str) -> Path:
        return self.matrices_dir / f"{safe_component(group_key)}_cells.csv"

    def dag_group_dir(self, group_key: str) -> Path:
        return self.dags_dir / safe_component(group_key)

    def cscn_object_path(self, group_key: str) -> Path:
        return self.objects_dir / f"{safe_component(group_key)}_cscn.pkl"

    def consensus_csv_path(self, group_key: str) -> Path:
        return self.consensus_dir / f"{safe_component(group_key)}.csv"
