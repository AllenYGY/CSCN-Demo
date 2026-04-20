import os
import pickle
import logging
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.decomposition import NMF

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
        self.ran_cache = {}
        self.bits_cache = {}
        self.data = None
        self.df = None
        self.kdtree = None
        self.loadings = None
        os.makedirs(output_dir, exist_ok=True)

    def clear_cache(self):
        self.ran_cache.clear()
        self.bits_cache.clear()

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

    def get_conditional_counts(self, genes, key_cell_idx, sigmoid_score=None):
        if sigmoid_score is None:
            sigmoid_score = self.sigmoid_score

        if self.debug:
            kdt_count = self.get_kdt_counts(genes, key_cell_idx, sigmoid_score)
            bit_count = self.get_bits_counts(genes, key_cell_idx, sigmoid_score)
            if kdt_count != bit_count:
                print(
                    f"kdt and bitset not equal: kdt{kdt_count} != bits{bit_count}, {genes}"
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
            condition_set = set(Z)
            count_z = self.get_conditional_counts(
                condition_set, key_cell_idx, sigmoid_score
            )
            if count_z <= 1:
                return False

            x_condition_set = condition_set.copy()
            x_condition_set.add(X)
            count_x_z = self.get_conditional_counts(
                x_condition_set, key_cell_idx, sigmoid_score
            )

            y_condition_set = condition_set.copy()
            y_condition_set.add(Y)
            count_y_z = self.get_conditional_counts(
                y_condition_set, key_cell_idx, sigmoid_score
            )

            xy_condition_set = condition_set.copy()
            xy_condition_set.add(X)
            xy_condition_set.add(Y)
            count_xy_z = self.get_conditional_counts(
                xy_condition_set, key_cell_idx, sigmoid_score
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
            print(f"ICT ERROR: {str(e)}")
            return True

    def run_pc(self, key_cell_idx):
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
                        f"{os.path.basename(filename)}"
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

    def run_core(self, data, usingNMF=False):
        _, col = data.shape

        if usingNMF:
            n_components = min(col, 100)
            nmf = NMF(n_components=n_components, random_state=42, max_iter=50000)
            factors = nmf.fit_transform(data)
            loadings = nmf.components_.T
            self.data = factors
            self.loadings = loadings
        else:
            self.data = data
            self.loadings = np.eye(col)
        self.df = pd.DataFrame(self.data, index=range(self.data.shape[0]))
        self.kdtree = KDT(list(self.data))
        self.clear_cache()
