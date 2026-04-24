import os
import pickle
import logging
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.decomposition import NMF
from sklearn.neighbors import NearestNeighbors

from .kdt import KDT


warnings.filterwarnings(
    "ignore",
    message=r"`pgmpy\.estimators\.StructureScore` is deprecated.*",
    category=FutureWarning,
)
warnings.filterwarnings(
    "ignore",
    message=r"PC is deprecated\..*",
    category=FutureWarning,
)
from pgmpy.estimators import PC


logging.getLogger("pgmpy").setLevel(logging.WARNING)


class CSCN:
    def __init__(
        self,
        output_dir="results",
        sigmoid_score=0.1,
        pc_var="stable",
        significance_level=0.01,
        max_cond_vars=20,
        use_bitmap=True,
        debug=False,
        show_progress=False,
        progress_interval=100,
        spatial_enabled=False,
        spatial_strategy="weighted_counts",
        spatial_mode="knn",
        spatial_k=8,
        spatial_radius=None,
        spatial_kernel="gaussian",
        spatial_bandwidth=None,
        spatial_lambda_expr=0.2,
        spatial_min_effective_neighbors=15,
    ):
        self.output_dir = output_dir
        self.sigmoid_score = sigmoid_score
        self.significance_level = significance_level
        self.max_cond_vars = max_cond_vars
        self.pc_var = pc_var
        self.use_bitmap = use_bitmap
        self.debug = debug
        self.show_progress = show_progress
        self.progress_interval = progress_interval
        self.using_nmf = False
        self.ran_cache = {}
        self.bits_cache = {}
        self.data = None
        self.raw_data = None
        self.df = None
        self.kdtree = None
        self.loadings = None
        self.ckm = None
        self.spatial_enabled = spatial_enabled
        self.spatial_strategy = spatial_strategy
        self.spatial_mode = spatial_mode
        self.spatial_k = spatial_k
        self.spatial_radius = spatial_radius
        self.spatial_kernel = spatial_kernel
        self.spatial_bandwidth = spatial_bandwidth
        self.spatial_lambda_expr = spatial_lambda_expr
        self.spatial_min_effective_neighbors = spatial_min_effective_neighbors
        self.spatial_coords = None
        self.spatial_neighbor_cache = {}
        self.spatial_weight_cache = {}
        os.makedirs(output_dir, exist_ok=True)

    def clear_cache(self):
        self.ran_cache.clear()
        self.bits_cache.clear()
        self.spatial_neighbor_cache.clear()
        self.spatial_weight_cache.clear()

    def get_ran_with_indices(self, gene_id, key_cell_idx, sigmoid_score=None):
        if sigmoid_score is None:
            sigmoid_score = self.sigmoid_score

        cache_key = (gene_id, key_cell_idx)
        if cache_key in self.ran_cache:
            return self.ran_cache[cache_key]

        gene_expressions = self.df.iloc[:, gene_id].values
        target_expression = gene_expressions[key_cell_idx]
        sorted_indices = np.argsort(gene_expressions)
        sorted_expressions = gene_expressions[sorted_indices]
        sorted_idx = np.where(sorted_expressions == target_expression)[0][0]

        total_cells = len(gene_expressions)
        window_size = int(sigmoid_score * total_cells)

        right_overflow = 0
        left_overflow = 0
        lower_bound = max(sorted_idx - window_size - right_overflow, 0)
        upper_bound = min(sorted_idx + left_overflow + window_size, total_cells - 1)
        expression_range = (
            sorted_expressions[lower_bound],
            sorted_expressions[upper_bound],
        )
        indices_in_range = sorted_indices[lower_bound : upper_bound + 1]

        result = {
            "expression_range": expression_range,
            "original_indices": indices_in_range.tolist(),
        }
        self.ran_cache[cache_key] = result
        return result

    def get_kdt_counts(self, genes, key_cell_idx, sigmoid_score):
        if len(genes) == 0:
            return 0
        template = self.df.iloc[key_cell_idx].tolist()
        ranges = []
        for idx, _ in enumerate(template):
            if idx in genes:
                ranges.append(
                    self.get_ran_with_indices(
                        idx, key_cell_idx, sigmoid_score
                    )["expression_range"]
                )
            else:
                ranges.append((1, 0))
        return self.kdtree.query_cnt(ranges)

    def get_bits(self, gene_id, key_cell_idx, sigmoid_score):
        n = self.data.shape[0] + 1

        def to_bitset(iterable):
            bitset = 0
            for x in iterable:
                if 0 <= x < n:
                    bitset |= 1 << x
                else:
                    raise ValueError(f"Element {x} is out of range [0, {n - 1}]")
            return bitset

        cache_key = (gene_id, key_cell_idx)
        if cache_key in self.bits_cache:
            return self.bits_cache[cache_key]
        original_indices = self.get_ran_with_indices(
            gene_id, key_cell_idx, sigmoid_score
        )["original_indices"]
        bitset = to_bitset(original_indices)
        self.bits_cache[cache_key] = bitset
        return bitset

    def get_bits_counts(self, genes, key_cell_idx, sigmoid_score=None):
        def intersection_size(bitsets):
            if not bitsets:
                return 0
            res = bitsets[0]
            for bitset in bitsets[1:]:
                res &= bitset
            return res.bit_count()

        if sigmoid_score is None:
            sigmoid_score = self.sigmoid_score
        bitsets = []
        for idx in genes:
            bitsets.append(self.get_bits(idx, key_cell_idx, sigmoid_score))
        return intersection_size(bitsets)

    def _bitset_to_indices(self, bitset):
        indices = []
        idx = 0
        while bitset:
            if bitset & 1:
                indices.append(idx)
            bitset >>= 1
            idx += 1
        return indices

    def build_spatial_neighbors(self):
        self.spatial_neighbor_cache.clear()
        if not self.spatial_enabled or self.spatial_coords is None:
            return
        coords = np.asarray(self.spatial_coords, dtype=np.float64)
        if coords.ndim != 2 or coords.shape[0] == 0 or coords.shape[1] < 2:
            return

        n_cells = coords.shape[0]
        nn = NearestNeighbors()
        nn.fit(coords[:, :2])
        if self.spatial_mode == "knn":
            n_neighbors = min(max(int(self.spatial_k), 1), n_cells)
            distances, indices = nn.kneighbors(coords[:, :2], n_neighbors=n_neighbors)
            for idx in range(n_cells):
                self.spatial_neighbor_cache[idx] = (
                    np.asarray(indices[idx], dtype=np.int64),
                    np.asarray(distances[idx], dtype=np.float64),
                )
            return

        if self.spatial_mode == "radius":
            radius = float(self.spatial_radius) if self.spatial_radius is not None else 0.0
            indices_list, distances_list = nn.radius_neighbors(
                coords[:, :2], radius=radius, sort_results=True
            )
            for idx in range(n_cells):
                self.spatial_neighbor_cache[idx] = (
                    np.asarray(indices_list[idx], dtype=np.int64),
                    np.asarray(distances_list[idx], dtype=np.float64),
                )

    def get_local_subset_indices(self, key_cell_idx):
        if not self.spatial_enabled or self.spatial_coords is None:
            return np.arange(self.data.shape[0], dtype=np.int64)
        cached = self.spatial_neighbor_cache.get(key_cell_idx)
        if cached is None:
            return np.arange(self.data.shape[0], dtype=np.int64)
        indices = np.asarray(cached[0], dtype=np.int64)
        if len(indices) == 0:
            return np.asarray([key_cell_idx], dtype=np.int64)
        return indices

    def build_local_df_for_key_cell(self, key_cell_idx):
        subset_indices = self.get_local_subset_indices(key_cell_idx)
        local_matrix = np.asarray(self.raw_data[subset_indices], dtype=np.float64)
        local_key_matches = np.where(subset_indices == key_cell_idx)[0]
        if len(local_key_matches) == 0:
            raise ValueError(f"Key cell {key_cell_idx} is missing from its local subset.")
        local_key_idx = int(local_key_matches[0])
        return pd.DataFrame(local_matrix, index=range(local_matrix.shape[0])), local_key_idx

    def _build_local_worker(self, subset_indices):
        local_worker = CSCN(
            output_dir=self.output_dir,
            sigmoid_score=self.sigmoid_score,
            pc_var=self.pc_var,
            significance_level=self.significance_level,
            max_cond_vars=self.max_cond_vars,
            use_bitmap=self.use_bitmap,
            debug=self.debug,
            show_progress=self.show_progress,
            progress_interval=self.progress_interval,
            spatial_enabled=False,
            spatial_strategy="weighted_counts",
            spatial_mode=self.spatial_mode,
            spatial_k=self.spatial_k,
            spatial_radius=self.spatial_radius,
            spatial_kernel=self.spatial_kernel,
            spatial_bandwidth=self.spatial_bandwidth,
            spatial_lambda_expr=self.spatial_lambda_expr,
            spatial_min_effective_neighbors=self.spatial_min_effective_neighbors,
        )
        local_spatial_coords = None
        if self.spatial_coords is not None:
            local_spatial_coords = np.asarray(self.spatial_coords[subset_indices], dtype=np.float64)
        local_worker.run_core(
            np.asarray(self.raw_data[subset_indices], dtype=np.float64),
            usingNMF=self.using_nmf,
            spatial_coords=local_spatial_coords,
        )
        return local_worker

    def get_spatial_weights(self, key_cell_idx):
        if key_cell_idx in self.spatial_weight_cache:
            return self.spatial_weight_cache[key_cell_idx]

        n_cells = self.data.shape[0]
        if not self.spatial_enabled or self.spatial_coords is None:
            weight_lookup = np.ones(n_cells, dtype=np.float64)
            result = (np.arange(n_cells, dtype=np.int64), weight_lookup)
            self.spatial_weight_cache[key_cell_idx] = result
            return result

        cached = self.spatial_neighbor_cache.get(key_cell_idx)
        if cached is None:
            weight_lookup = np.ones(n_cells, dtype=np.float64)
            result = (np.arange(n_cells, dtype=np.int64), weight_lookup)
            self.spatial_weight_cache[key_cell_idx] = result
            return result

        indices, distances = cached
        if len(indices) <= 1:
            weight_lookup = np.ones(n_cells, dtype=np.float64)
            result = (np.arange(n_cells, dtype=np.int64), weight_lookup)
            self.spatial_weight_cache[key_cell_idx] = result
            return result

        if self.spatial_kernel == "binary":
            neighbor_weights = np.ones(len(indices), dtype=np.float64)
        else:
            bandwidth = self.spatial_bandwidth
            if bandwidth is None:
                positive = distances[distances > 0]
                bandwidth = float(np.median(positive)) if len(positive) else 1.0
            bandwidth = max(float(bandwidth), 1e-12)
            neighbor_weights = np.exp(-(distances**2) / (2 * (bandwidth**2)))
            max_weight = np.max(neighbor_weights) if len(neighbor_weights) else 1.0
            if max_weight > 0:
                neighbor_weights = neighbor_weights / max_weight

        weight_lookup = np.full(n_cells, float(self.spatial_lambda_expr), dtype=np.float64)
        weight_lookup[indices] = (
            float(self.spatial_lambda_expr)
            + (1.0 - float(self.spatial_lambda_expr)) * neighbor_weights
        )
        result = (indices, weight_lookup)
        self.spatial_weight_cache[key_cell_idx] = result
        return result

    def get_effective_sample_size(self, weights):
        weights = np.asarray(weights, dtype=np.float64)
        denom = np.sum(weights**2)
        if denom <= 0:
            return 0.0
        total = np.sum(weights)
        return float((total**2) / denom)

    def get_weighted_conditional_counts(self, genes, key_cell_idx, sigmoid_score=None):
        if sigmoid_score is None:
            sigmoid_score = self.sigmoid_score
        _, weight_lookup = self.get_spatial_weights(key_cell_idx)
        if len(genes) == 0:
            return float(np.sum(weight_lookup))

        intersected = None
        for idx in genes:
            gene_bits = self.get_bits(idx, key_cell_idx, sigmoid_score)
            intersected = gene_bits if intersected is None else (intersected & gene_bits)
            if intersected == 0:
                return 0.0

        matched_indices = self._bitset_to_indices(intersected)
        if not matched_indices:
            return 0.0
        return float(np.sum(weight_lookup[matched_indices]))

    def _should_use_spatial_counts(self, key_cell_idx):
        if (
            not self.spatial_enabled
            or self.spatial_coords is None
            or self.spatial_strategy != "weighted_counts"
        ):
            return False, None, None
        neighbor_indices, weight_lookup = self.get_spatial_weights(key_cell_idx)
        n_eff = self.get_effective_sample_size(weight_lookup)
        if len(neighbor_indices) <= 1 or n_eff < self.spatial_min_effective_neighbors:
            return False, weight_lookup, n_eff
        return True, weight_lookup, n_eff

    def get_conditional_counts(self, genes, key_cell_idx, sigmoid_score=None):
        if sigmoid_score is None:
            sigmoid_score = self.sigmoid_score

        if self.debug:
            kdt_count = self.get_kdt_counts(genes, key_cell_idx, sigmoid_score)
            bit_count = self.get_bits_counts(genes, key_cell_idx, sigmoid_score)
            if kdt_count != bit_count:
                print(
                    f"kdt and bitset not equal: kdt{kdt_count} != bits{bit_count}, {genes}",
                    flush=True,
                )
            if self.use_bitmap:
                return bit_count
            return kdt_count

        if self.use_bitmap:
            return self.get_bits_counts(genes, key_cell_idx, sigmoid_score)
        return self.get_kdt_counts(genes, key_cell_idx, sigmoid_score)

    def conditional_independence_test(
        self,
        X,
        Y,
        Z,
        data,
        independencies,
        key_cell_idx=0,
        significance_level=None,
        sigmoid_score=None,
    ):
        if significance_level is None:
            significance_level = self.significance_level
        if sigmoid_score is None:
            sigmoid_score = self.sigmoid_score

        try:
            use_spatial_counts, weight_lookup, n_eff = self._should_use_spatial_counts(
                key_cell_idx
            )
            condition_set = set(Z)
            count_z = (
                self.get_weighted_conditional_counts(
                    condition_set, key_cell_idx, sigmoid_score
                )
                if use_spatial_counts
                else self.get_conditional_counts(condition_set, key_cell_idx, sigmoid_score)
            )
            if count_z <= 1:
                return False

            x_condition_set = condition_set.copy()
            x_condition_set.add(X)
            count_x_z = (
                self.get_weighted_conditional_counts(
                    x_condition_set, key_cell_idx, sigmoid_score
                )
                if use_spatial_counts
                else self.get_conditional_counts(x_condition_set, key_cell_idx, sigmoid_score)
            )

            y_condition_set = condition_set.copy()
            y_condition_set.add(Y)
            count_y_z = (
                self.get_weighted_conditional_counts(
                    y_condition_set, key_cell_idx, sigmoid_score
                )
                if use_spatial_counts
                else self.get_conditional_counts(y_condition_set, key_cell_idx, sigmoid_score)
            )

            xy_condition_set = condition_set.copy()
            xy_condition_set.add(X)
            xy_condition_set.add(Y)
            count_xy_z = (
                self.get_weighted_conditional_counts(
                    xy_condition_set, key_cell_idx, sigmoid_score
                )
                if use_spatial_counts
                else self.get_conditional_counts(xy_condition_set, key_cell_idx, sigmoid_score)
            )

            rho = ((count_z * count_xy_z) - (count_x_z * count_y_z)) / (count_z**2)

            epsilon = 1e-10
            variance_numerator = (
                count_x_z
                * count_y_z
                * (count_z - count_x_z)
                * (count_z - count_y_z)
            )
            variance_numerator = np.clip(variance_numerator, epsilon, None)
            if use_spatial_counts:
                if n_eff is None or n_eff <= 1 + epsilon:
                    return False
                variance_denominator = (count_z**4) * max(n_eff - 1, epsilon)
            else:
                variance_denominator = (count_z**4) * (count_z - 1)
            std_deviation = np.sqrt(
                variance_numerator / (variance_denominator + epsilon)
            )

            z_score = np.divide(
                rho,
                std_deviation,
                out=np.zeros_like(rho),
                where=(std_deviation > epsilon),
            )
            p_value = 2 * (1 - norm.cdf(abs(z_score)))

            if std_deviation <= epsilon and np.abs(rho) > epsilon:
                p_value = 0

            return p_value > significance_level

        except Exception as e:
            print(f"ICT ERROR: {str(e)}", flush=True)
            return True

    def run_pc(self, key_cell_idx):
        if self.spatial_enabled and self.spatial_strategy == "local_knn_subset":
            subset_indices = self.get_local_subset_indices(key_cell_idx)
            local_worker = self._build_local_worker(subset_indices)
            local_key_idx = int(np.where(subset_indices == key_cell_idx)[0][0])
            return local_worker.run_pc(local_key_idx)
        pc = PC(self.df)
        return pc.estimate(
            variant=self.pc_var,
            ci_test=self.conditional_independence_test,
            significance_level=self.significance_level,
            max_cond_vars=self.max_cond_vars,
            return_type="dag",
            key_cell_idx=key_cell_idx,
            sigmoid_score=self.sigmoid_score,
            show_progress=self.show_progress,
        )

    def run_pc_and_save(self, task_id):
        result = self.run_pc(task_id)
        filename = os.path.join(self.output_dir, f"result_{task_id}.pkl")
        with open(filename, "wb") as handle:
            pickle.dump(result, handle)
        return filename

    def run_pc_concurrently(
        self,
        max_workers=None,
        progress_interval=None,
        progress_label=None,
    ):
        total = len(self.df)
        if progress_interval is None:
            progress_interval = self.progress_interval

        results = [None] * total
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(self.run_pc_and_save, i): i for i in range(total)
            }
            completed = 0
            for future in as_completed(futures):
                task_id = futures[future]
                filename = future.result()
                results[task_id] = filename
                completed += 1

                if progress_interval and (
                    completed == 1
                    or completed % progress_interval == 0
                    or completed == total
                ):
                    label = f"{progress_label}: " if progress_label else ""
                    print(
                        f"[CSCN] {label}saved DAG {completed}/{total} -> "
                        f"{os.path.basename(filename)}",
                        flush=True,
                    )
        return results

    def load_all_dags(self):
        dags = []
        for filename in os.listdir(self.output_dir):
            if filename.endswith(".pkl"):
                with open(os.path.join(self.output_dir, filename), "rb") as handle:
                    dags.append(
                        (
                            int(filename[7:].split(".pkl")[0]),
                            pickle.load(handle),
                        )
                    )
        return dags

    @staticmethod
    def save_to_file(instance, filename):
        with open(filename, "wb") as handle:
            pickle.dump(instance, handle)

    @staticmethod
    def load_from_file(filename):
        with open(filename, "rb") as handle:
            return pickle.load(handle)

    def _resolve_ckm_projection(self):
        if self.data is None:
            raise RuntimeError("CKM requires CSCN data. Run run_core() first.")

        feature_count = int(self.data.shape[1])
        if self.loadings is None:
            return np.eye(feature_count, dtype=float)

        loadings = np.asarray(self.loadings, dtype=float)
        if loadings.ndim != 2:
            raise ValueError("CSCN loadings must be a 2-dimensional matrix.")

        if int(loadings.shape[0]) == feature_count:
            return loadings
        if int(loadings.shape[1]) == feature_count:
            return loadings.T
        raise ValueError(
            "CKM projection could not align CSCN loadings with the current feature space."
        )

    @staticmethod
    def _prepare_graph_for_ckm(dag, node_count: int):
        graph = nx.DiGraph()
        graph.add_nodes_from(range(node_count))
        graph.add_edges_from(dag.edges())
        invalid_nodes = sorted(
            node for node in graph.nodes() if int(node) < 0 or int(node) >= node_count
        )
        if invalid_nodes:
            raise ValueError(f"CKM DAG contains invalid node ids: {invalid_nodes[:10]}")
        return graph

    def compute_ckm(
        self,
        dags=None,
        alpha: float = 0.05,
        beta_transform: str = "log1p",
        save_path: str | None = None,
        strict: bool = True,
    ):
        if self.data is None:
            raise RuntimeError("CKM requires CSCN data. Run run_core() first.")

        projection = self._resolve_ckm_projection()
        latent_dim = int(self.data.shape[1])
        gene_dim = int(projection.shape[1])
        n_cells = int(self.data.shape[0])
        ckm = np.zeros((n_cells, gene_dim), dtype=float)

        if dags is None:
            dags = self.load_all_dags()
        dags = list(dags)
        if not dags:
            raise RuntimeError("No DAGs were provided for CKM computation.")

        dag_map = {int(cell_idx): dag for cell_idx, dag in dags}
        missing = sorted(set(range(n_cells)) - set(dag_map))
        if strict and missing:
            raise ValueError(
                f"CKM is missing DAGs for {len(missing)} cells, e.g. {missing[:10]}"
            )

        for cell_idx, dag in sorted(dag_map.items()):
            if cell_idx < 0 or cell_idx >= n_cells:
                raise ValueError(f"CKM received out-of-range cell index: {cell_idx}")

            base_beta = np.asarray(self.data[cell_idx], dtype=float)
            if beta_transform == "log1p":
                beta = np.log1p(base_beta)
            elif beta_transform == "identity":
                beta = base_beta
            else:
                raise ValueError(f"Unsupported CKM beta_transform: {beta_transform}")

            graph = self._prepare_graph_for_ckm(dag, latent_dim)
            katz = nx.katz_centrality(
                graph,
                beta={node_idx: float(beta[node_idx]) for node_idx in range(latent_dim)},
                alpha=float(alpha),
            )
            katz_vec = np.asarray(
                [float(katz[node_idx]) for node_idx in range(latent_dim)],
                dtype=float,
            )
            ckm[cell_idx, :] = np.log1p(katz_vec @ projection)

        self.ckm = ckm
        if save_path:
            np.save(save_path, ckm)
        return ckm

    def run_core(self, data, usingNMF=False, spatial_coords=None):
        self.using_nmf = usingNMF
        self.raw_data = np.asarray(data, dtype=np.float64)
        _, col = data.shape

        if usingNMF:
            n_components = min(col, 100)
            nmf = NMF(n_components=n_components, random_state=42, max_iter=50000)
            factors = nmf.fit_transform(self.raw_data)
            loadings = nmf.components_.T
            self.data = factors
            self.loadings = loadings
        else:
            self.data = self.raw_data
            self.loadings = np.eye(col)
        self.df = pd.DataFrame(self.data, index=range(self.data.shape[0]))
        self.kdtree = KDT(list(self.data))
        self.spatial_coords = (
            None
            if spatial_coords is None
            else np.asarray(spatial_coords, dtype=np.float64)
        )
        self.clear_cache()
        self.build_spatial_neighbors()
