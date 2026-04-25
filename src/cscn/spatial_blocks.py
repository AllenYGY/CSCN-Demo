from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.cluster import DBSCAN


@dataclass(frozen=True)
class SpatialBlock:
    block_id: str
    cluster_id: int
    center_x: float
    center_y: float
    min_x: float
    max_x: float
    min_y: float
    max_y: float
    cell_ids: tuple[str, ...]


@dataclass(frozen=True)
class CellBlockAssignment:
    cell_id: str
    primary_block_id: str | None
    covering_block_ids: tuple[str, ...]
    halo_block_ids: tuple[str, ...]


def _resolve_density_clusters(
    coords: np.ndarray,
    *,
    eps: float,
    min_samples: int,
) -> list[np.ndarray]:
    labels = DBSCAN(eps=eps, min_samples=min_samples).fit_predict(coords)
    clusters = [np.flatnonzero(labels == label) for label in sorted(set(labels)) if label >= 0]
    if clusters:
        return clusters
    return [np.arange(coords.shape[0], dtype=np.int64)]


def _window_members(
    coords: np.ndarray,
    *,
    center_x: float,
    center_y: float,
    width: float,
    height: float,
) -> np.ndarray:
    half_w = width / 2.0
    half_h = height / 2.0
    mask = (
        (coords[:, 0] >= center_x - half_w)
        & (coords[:, 0] <= center_x + half_w)
        & (coords[:, 1] >= center_y - half_h)
        & (coords[:, 1] <= center_y + half_h)
    )
    return np.flatnonzero(mask)


def generate_adaptive_blocks(
    *,
    cell_ids: list[str],
    coords: np.ndarray,
    eps: float,
    min_samples: int,
    min_cells_per_block: int,
    overlap_min: float,
) -> list[SpatialBlock]:
    if coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("Adaptive blocks require 2D coordinates.")
    if coords.shape[0] == 0:
        return []

    clusters = _resolve_density_clusters(coords, eps=eps, min_samples=min_samples)
    blocks: list[SpatialBlock] = []
    seen_memberships: set[tuple[int, ...]] = set()

    for cluster_idx, cluster_members in enumerate(clusters):
        cluster_coords = coords[cluster_members]
        min_x = float(np.min(cluster_coords[:, 0]))
        max_x = float(np.max(cluster_coords[:, 0]))
        min_y = float(np.min(cluster_coords[:, 1]))
        max_y = float(np.max(cluster_coords[:, 1]))
        width = max(max_x - min_x, float(eps))
        height = max(max_y - min_y, float(eps))
        center_x = (min_x + max_x) / 2.0
        center_y = (min_y + max_y) / 2.0
        shift_x = width * max(1.0 - float(overlap_min), 0.0)
        shift_y = height * max(1.0 - float(overlap_min), 0.0)
        centers = [
            (center_x, center_y),
            (center_x - shift_x, center_y),
            (center_x + shift_x, center_y),
            (center_x, center_y - shift_y),
            (center_x, center_y + shift_y),
        ]

        for local_idx, (window_center_x, window_center_y) in enumerate(centers):
            members = _window_members(
                coords,
                center_x=window_center_x,
                center_y=window_center_y,
                width=width,
                height=height,
            )
            membership_key = tuple(int(item) for item in members.tolist())
            if len(members) < min_cells_per_block or membership_key in seen_memberships:
                continue
            seen_memberships.add(membership_key)
            blocks.append(
                SpatialBlock(
                    block_id=f"cluster_{cluster_idx}_block_{local_idx}",
                    cluster_id=cluster_idx,
                    center_x=float(window_center_x),
                    center_y=float(window_center_y),
                    min_x=float(window_center_x - width / 2.0),
                    max_x=float(window_center_x + width / 2.0),
                    min_y=float(window_center_y - height / 2.0),
                    max_y=float(window_center_y + height / 2.0),
                    cell_ids=tuple(cell_ids[index] for index in members.tolist()),
                )
            )

    if blocks:
        return blocks

    if coords.shape[0] < min_cells_per_block:
        return []

    # Fall back to one block covering all cells if every candidate was filtered out.
    global_min_x = float(np.min(coords[:, 0]))
    global_max_x = float(np.max(coords[:, 0]))
    global_min_y = float(np.min(coords[:, 1]))
    global_max_y = float(np.max(coords[:, 1]))
    return [
        SpatialBlock(
            block_id="cluster_0_block_0",
            cluster_id=0,
            center_x=(global_min_x + global_max_x) / 2.0,
            center_y=(global_min_y + global_max_y) / 2.0,
            min_x=global_min_x,
            max_x=global_max_x,
            min_y=global_min_y,
            max_y=global_max_y,
            cell_ids=tuple(cell_ids),
        )
    ]


def build_cell_block_assignments(
    *,
    cell_ids: list[str],
    coords: np.ndarray,
    blocks: list[SpatialBlock],
    halo_neighbor_blocks: int,
) -> dict[str, CellBlockAssignment]:
    if not blocks:
        return {
            cell_id: CellBlockAssignment(
                cell_id=cell_id,
                primary_block_id=None,
                covering_block_ids=(),
                halo_block_ids=(),
            )
            for cell_id in cell_ids
        }

    centers = np.asarray([[block.center_x, block.center_y] for block in blocks], dtype=float)
    assignments: dict[str, CellBlockAssignment] = {}
    for index, cell_id in enumerate(cell_ids):
        point = coords[index]
        covering: list[str] = []
        covering_indexes: list[int] = []
        for block_index, block in enumerate(blocks):
            if (
                block.min_x <= point[0] <= block.max_x
                and block.min_y <= point[1] <= block.max_y
            ):
                covering.append(block.block_id)
                covering_indexes.append(block_index)

        distances = np.linalg.norm(centers - point[:2], axis=1)
        ranked_indexes = np.argsort(distances).tolist()
        if covering_indexes:
            primary_index = min(covering_indexes, key=lambda item: distances[item])
        else:
            primary_index = int(ranked_indexes[0])
        halo_indexes = [item for item in ranked_indexes if item != primary_index][
            : max(int(halo_neighbor_blocks), 0)
        ]
        assignments[cell_id] = CellBlockAssignment(
            cell_id=cell_id,
            primary_block_id=blocks[primary_index].block_id,
            covering_block_ids=tuple(
                block_id
                for _, block_id in sorted(
                    zip(
                        [distances[item] for item in covering_indexes],
                        [blocks[item].block_id for item in covering_indexes],
                    )
                )
            ),
            halo_block_ids=tuple(blocks[item].block_id for item in halo_indexes),
        )
    return assignments
