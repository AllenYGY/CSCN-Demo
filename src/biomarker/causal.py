import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class SimpleCausalInference:
    def __init__(
        self,
        dag,
        data,
        treatment_node,
        outcome_node,
        confounders=None,
        random_state=42,
    ):
        self.random_state = random_state
        np.random.seed(random_state)

        self.dag = dag
        self.data = data
        self.treatment_node = treatment_node
        self.outcome_node = outcome_node

        if confounders is None:
            self.confounders = self.auto_identify_confounders()
        else:
            self.confounders = confounders

        self._validate_inputs()

        print(f"DAG: {len(dag.nodes())} nodes, {len(dag.edges())} edges")
        print(f"Data: {data.shape}")
        print(f"Treatment: {treatment_node}")
        print(f"Outcome: {outcome_node}")
        print(f"Confounders: {self.confounders}")

    def auto_identify_confounders(self):
        if self.dag is None:
            return []
        treatment_parents = set(self.dag.predecessors(self.treatment_node))
        outcome_parents = set(self.dag.predecessors(self.outcome_node))
        outcome_parents.discard(self.treatment_node)
        return list(treatment_parents.intersection(outcome_parents))

    def _validate_inputs(self):
        missing_vars = []
        if self.treatment_node not in self.data.columns:
            missing_vars.append(self.treatment_node)
        if self.outcome_node not in self.data.columns:
            missing_vars.append(self.outcome_node)
        for confounder in self.confounders:
            if confounder not in self.data.columns:
                missing_vars.append(confounder)

        if missing_vars:
            raise ValueError(f"Variables not found in data: {missing_vars}")

    def compute_ate(self):
        print("\n=== Computing ATE ===")

        adjustment_set = self.confounders
        print(f"Adjustment set: {adjustment_set}")

        x0 = self.data[self.data[self.outcome_node] == 0][self.treatment_node].mean()
        x1 = self.data[self.data[self.outcome_node] == 1][self.treatment_node].mean()
        print(f"Intervention levels: x0={x0:.3f}, x1={x1:.3f}")

        ate_adjustment = self._adjustment_formula(adjustment_set, x0, x1)
        ate_gformula = self._g_formula(adjustment_set, x0, x1)

        print(f"\n=== Results ===")
        print(f"Adjustment Formula ATE: {ate_adjustment:.4f}")
        print(f"G-formula ATE: {ate_gformula:.4f}")
        print(f"Difference: {abs(ate_adjustment - ate_gformula):.4f}")

        return {"adjustment_formula": ate_adjustment, "g_formula": ate_gformula}

    def _adjustment_formula(self, adjustment_set, x0, x1):
        print("\n--- Adjustment Formula ---")

        if not adjustment_set:
            print("No confounders, computing naive estimate")
            x0 = self.data[self.data[self.outcome_node] == 0][self.treatment_node].mean()
            x1 = self.data[self.data[self.outcome_node] == 1][self.treatment_node].mean()

            y0_mean = self.data[x0][self.outcome_node].mean()
            y1_mean = self.data[x1][self.outcome_node].mean()
            return y1_mean - y0_mean

        n_bins = 3
        data_binned = self.data.copy()

        for confounder in adjustment_set:
            try:
                bins = pd.qcut(
                    self.data[confounder],
                    q=n_bins,
                    labels=False,
                    duplicates="drop",
                )
                data_binned[f"{confounder}_bin"] = bins
            except Exception:
                median_val = self.data[confounder].median()
                data_binned[f"{confounder}_bin"] = (
                    self.data[confounder] > median_val
                ).astype(int)

        bin_cols = [f"{confounder}_bin" for confounder in adjustment_set]
        causal_effects = []

        for x_val in [x0, x1]:
            total_expectation = 0
            total_probability = 0

            unique_strata = data_binned[bin_cols].drop_duplicates()

            for _, stratum in unique_strata.iterrows():
                mask = True
                for col in bin_cols:
                    mask &= data_binned[col] == stratum[col]

                stratum_data = data_binned[mask]

                if len(stratum_data) < 3:
                    continue

                stratum_prob = len(stratum_data) / len(self.data)
                conditional_expectation = self._estimate_conditional_expectation(
                    stratum_data, x_val
                )

                total_expectation += conditional_expectation * stratum_prob
                total_probability += stratum_prob

            if total_probability > 0:
                causal_effects.append(total_expectation / total_probability)
            else:
                causal_effects.append(0)

        ate = causal_effects[1] - causal_effects[0]
        print(f"E[Y|do(X=x0={x0:.3f})] = {causal_effects[0]:.4f}")
        print(f"E[Y|do(X=x1={x1:.3f})] = {causal_effects[1]:.4f}")

        return ate

    def _estimate_conditional_expectation(self, stratum_data, x_val):
        X = stratum_data[self.treatment_node].values
        Y = stratum_data[self.outcome_node].values

        if len(X) == 0:
            return 0

        if len(X) >= 5:
            distances = np.abs(X - x_val)
            k = min(5, len(X))
            nearest_indices = np.argsort(distances)[:k]
            return Y[nearest_indices].mean()
        return Y.mean()

    def _g_formula(self, adjustment_set, x0, x1):
        print("\n--- G-formula ---")

        if not adjustment_set:
            print("No confounders, computing naive estimate")
            x0 = self.data[self.data[self.outcome_node] == 0][self.treatment_node].mean()
            x1 = self.data[self.data[self.outcome_node] == 1][self.treatment_node].mean()
            y0_mean = self.data[x0][self.outcome_node].mean()
            y1_mean = self.data[x1][self.outcome_node].mean()
            return y1_mean - y0_mean

        data_stratified = self.data.copy()
        n_bins = 3

        for confounder in adjustment_set:
            try:
                data_stratified[f"{confounder}_stratum"] = pd.qcut(
                    self.data[confounder],
                    q=n_bins,
                    labels=False,
                    duplicates="drop",
                )
            except Exception:
                median_val = self.data[confounder].median()
                data_stratified[f"{confounder}_stratum"] = (
                    self.data[confounder] > median_val
                ).astype(int)

        stratum_cols = [f"{confounder}_stratum" for confounder in adjustment_set]

        def compute_expectation(x_val):
            total_expectation = 0
            total_weight = 0

            unique_strata = data_stratified[stratum_cols].drop_duplicates()

            for _, stratum in unique_strata.iterrows():
                mask = True
                for col in stratum_cols:
                    mask &= data_stratified[col] == stratum[col]

                stratum_data = data_stratified[mask]

                if len(stratum_data) < 3:
                    continue

                stratum_prob = len(stratum_data) / len(self.data)
                conditional_expectation = (
                    self._estimate_conditional_expectation_gformula(stratum_data, x_val)
                )

                total_expectation += conditional_expectation * stratum_prob
                total_weight += stratum_prob

            return total_expectation / total_weight if total_weight > 0 else 0

        prob_low = compute_expectation(x0)
        prob_high = compute_expectation(x1)

        ate = prob_high - prob_low
        print(f"E[Y|do(X=x0={x0:.3f})] = {prob_low:.4f}")
        print(f"E[Y|do(X=x1={x1:.3f})] = {prob_high:.4f}")

        return ate

    def _estimate_conditional_expectation_gformula(self, stratum_data, x_val):
        X = stratum_data[self.treatment_node].values
        Y = stratum_data[self.outcome_node].values

        if len(X) == 0:
            return 0

        if len(X) < 5:
            return Y.mean()

        try:
            from sklearn.linear_model import LinearRegression

            model = LinearRegression()
            model.fit(X.reshape(-1, 1), Y)
            prediction = model.predict([[x_val]])[0]

            if set(Y).issubset({0, 1}):
                prediction = np.clip(prediction, 0, 1)

            return prediction

        except Exception:
            distances = np.abs(X - x_val)
            k = min(3, len(X))
            nearest_indices = np.argsort(distances)[:k]
            return Y[nearest_indices].mean()

    def _visualize_results(self, x0, x1, ate_adj, ate_gf):
        plt.close("all")
        fig, axes = plt.subplots(1, 2, figsize=(16, 6), facecolor="whitesmoke")

        axes[0].hist(
            self.data[self.treatment_node],
            bins=30,
            alpha=0.6,
            color="mediumseagreen",
            edgecolor="white",
        )
        axes[0].axvline(
            x0,
            color="dodgerblue",
            linestyle="--",
            linewidth=2,
            label=f"x0 = {x0:.3f}",
        )
        axes[0].axvline(
            x1,
            color="darkorange",
            linestyle="--",
            linewidth=2,
            label=f"x1 = {x1:.3f}",
        )
        axes[0].set_xlabel(f"{self.treatment_node}", fontsize=14, color="dimgray")
        axes[0].set_ylabel("Frequency", fontsize=14, color="dimgray")
        axes[0].set_title(
            f"Treatment Gene {self.treatment_node} Distribution",
            fontsize=16,
            fontweight="bold",
            color="dimgray",
        )
        axes[0].legend(fontsize=12, frameon=False)
        axes[0].grid(True, alpha=0.4, linestyle="--")
        axes[0].tick_params(axis="both", which="major", labelsize=12, color="dimgray")
        axes[0].spines["top"].set_visible(False)
        axes[0].spines["right"].set_visible(False)

        methods = ["Adjustment\nFormula", "G-formula"]
        ates = [ate_adj, ate_gf]
        colors = ["skyblue", "lightcoral"]

        axes[1].bar(
            methods,
            ates,
            color=colors,
            alpha=0.8,
            edgecolor="black",
            linewidth=1.2,
        )
        axes[1].set_ylabel("Average Treatment Effect", fontsize=14, color="dimgray")
        axes[1].set_title(
            "ATE Comparison",
            fontsize=16,
            fontweight="bold",
            color="dimgray",
        )
        axes[1].grid(True, alpha=0.4, linestyle="--", axis="y")
        axes[1].tick_params(axis="both", which="major", labelsize=12, color="dimgray")
        axes[1].spines["top"].set_visible(False)
        axes[1].spines["right"].set_visible(False)
        axes[1].axhline(y=0, color="black", linestyle="-", alpha=0.3)

        plt.tight_layout()
        plt.show()


def run_causal_analysis(dag, data, treatment, outcome, confounders=None):
    try:
        framework = SimpleCausalInference(
            dag=dag,
            data=data,
            treatment_node=treatment,
            outcome_node=outcome,
            confounders=confounders,
        )
        results = framework.compute_ate()
        return {"success": True, "results": results, "framework": framework}
    except Exception as e:
        print(f"Error in causal analysis: {e}")
        return {"success": False, "error": str(e), "results": None}
