"""
modules/differential_expression.py
===================================
Differential Expression Analysis (DEA): Normal vs Tumor.

What this module does:
    For every gene in the matrix, it compares expression values between
    Normal samples and Tumor samples using a statistical test, then
    corrects for the fact that we are testing thousands of genes at once
    (multiple testing correction).

Statistics computed per gene:
    log2FoldChange  — Direction and magnitude of expression change.
                      Positive = higher in tumor (upregulated / oncogene).
                      Negative = lower in tumor  (silenced / tumor suppressor).
    p-value         — Probability of observing this difference by chance.
    padj (FDR)      — p-value corrected for multiple testing (Benjamini-Hochberg).
    significant     — True if padj < threshold AND |log2FC| > threshold.
    regulation      — 'up', 'down', or 'not_significant'.

Output files:
    de_results_full.csv        — All genes with statistics.
    de_results_significant.csv — Only statistically significant DEGs.
    de_summary.json            — Machine-readable summary.
"""

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


class DifferentialExpression:
    """
    Differential expression analysis between Normal and Tumor conditions.

    Args:
        output_dir (str | Path): Where to save result files.
        config     (dict):       Pipeline configuration dictionary.
    """

    def __init__(self, output_dir="results/de_results/", config: dict = None):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.config = config or {}
        de_cfg = self.config.get("differential_expression", {})

        # Significance thresholds (configurable)
        self.padj_threshold = de_cfg.get("padj_threshold", 0.05)
        self.lfc_threshold  = de_cfg.get("lfc_threshold",   1.0)
        self.test_method    = de_cfg.get("test_method",   "welch")

    # ── Public API ──────────────────────────────────────────────────────────

    def run(self, norm_df: pd.DataFrame, metadata_df: pd.DataFrame) -> pd.DataFrame:
        """
        Run differential expression analysis.

        Args:
            norm_df     (pd.DataFrame): Normalized & log-transformed matrix (genes x samples).
            metadata_df (pd.DataFrame): Sample metadata with 'condition' column.

        Returns:
            pd.DataFrame: Full results table sorted by adjusted p-value.
        """
        # Split samples by condition
        normal_samples = metadata_df[metadata_df["condition"] == "normal"].index.tolist()
        tumor_samples  = metadata_df[metadata_df["condition"] == "tumor"].index.tolist()

        # Keep only samples present in the expression matrix
        normal_samples = [s for s in normal_samples if s in norm_df.columns]
        tumor_samples  = [s for s in tumor_samples  if s in norm_df.columns]

        logger.info(f"  Normal samples : {len(normal_samples)} — {normal_samples}")
        logger.info(f"  Tumor samples  : {len(tumor_samples)}  — {tumor_samples}")
        logger.info(f"  Statistical test: {self.test_method}")

        normal_data = norm_df[normal_samples]
        tumor_data  = norm_df[tumor_samples]

        # Compute statistics gene-by-gene
        results = self._compute_statistics(normal_data, tumor_data)

        # Correct for multiple testing (Benjamini-Hochberg FDR)
        results = self._apply_fdr_correction(results)

        # Classify genes as up, down, or not significant
        results = self._classify_genes(results)

        # Write output files
        self._save_results(results)

        # Print summary
        n_sig  = results["significant"].sum()
        n_up   = (results["regulation"] == "up").sum()
        n_down = (results["regulation"] == "down").sum()
        logger.info(
            f"  Significant DEGs : {n_sig:,} "
            f"(up: {n_up}, down: {n_down}) "
            f"[padj < {self.padj_threshold}, |LFC| > {self.lfc_threshold}]"
        )
        return results

    # ── Statistics ──────────────────────────────────────────────────────────

    def _compute_statistics(
        self,
        normal_df: pd.DataFrame,
        tumor_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Compute per-gene statistics.

        For each gene:
          - log2 Fold Change  = mean(tumor) - mean(normal)  [in log2 space]
          - p-value           = Welch t-test or Mann-Whitney U test
        """
        records = []
        for gene in normal_df.index:
            n_vals = normal_df.loc[gene].values.astype(float)
            t_vals = tumor_df.loc[gene].values.astype(float)

            mean_n = float(np.mean(n_vals))
            mean_t = float(np.mean(t_vals))
            std_n  = float(np.std(n_vals, ddof=1))
            std_t  = float(np.std(t_vals, ddof=1))

            # In log2 space, subtraction = log2 fold change
            lfc = mean_t - mean_n

            if self.test_method == "welch":
                # Welch t-test: does not assume equal variance between groups
                stat, pval = stats.ttest_ind(t_vals, n_vals, equal_var=False)
            else:
                # Mann-Whitney U: non-parametric alternative (rank-based)
                stat, pval = stats.mannwhitneyu(
                    t_vals, n_vals, alternative="two-sided"
                )

            records.append({
                "gene":        gene,
                "log2FC":      lfc,
                "mean_normal": mean_n,
                "mean_tumor":  mean_t,
                "std_normal":  std_n,
                "std_tumor":   std_t,
                "stat":        float(stat) if not np.isnan(stat) else 0.0,
                "pvalue":      float(pval) if not np.isnan(pval) else 1.0,
            })

        results = pd.DataFrame(records).set_index("gene")
        results["pvalue"] = results["pvalue"].fillna(1.0).clip(0, 1)
        return results

    def _apply_fdr_correction(self, results: pd.DataFrame) -> pd.DataFrame:
        """
        Apply Benjamini-Hochberg False Discovery Rate (FDR) correction.

        Why do we need this?
            If we test 5,000 genes at p < 0.05, we expect ~250 false positives
            by pure chance. FDR correction controls the *proportion* of false
            discoveries among all genes we call significant.

        The corrected p-value (padj) represents: if we call all genes with
        padj < 0.05 significant, then at most 5% of them are false positives.
        """
        pvals          = results["pvalue"].values
        results["padj"] = self._benjamini_hochberg(pvals)
        return results

    @staticmethod
    def _benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
        """
        Benjamini-Hochberg procedure for FDR control.

        Steps:
            1. Sort p-values from smallest to largest and track their ranks.
            2. Adjusted p-value[i] = p-value[i] x (n / rank[i])
            3. Enforce monotonicity: walk backward so padj[i] <= padj[i+1].
            4. Clip to [0, 1].
        """
        n     = len(pvals)
        order = np.argsort(pvals)
        ranks = np.empty(n, dtype=int)
        ranks[order] = np.arange(1, n + 1)

        padj         = pvals * n / ranks

        # Enforce monotonicity (step-down adjustment)
        ordered_padj = padj[order]
        for i in range(n - 2, -1, -1):
            ordered_padj[i] = min(ordered_padj[i], ordered_padj[i + 1])

        result         = np.empty(n)
        result[order]  = ordered_padj
        return result.clip(0, 1)

    def _classify_genes(self, results: pd.DataFrame) -> pd.DataFrame:
        """
        Label genes as upregulated, downregulated, or not significant.

        A gene is significant if BOTH conditions are met:
            1. padj < padj_threshold    (usually 0.05)
            2. |log2FC| > lfc_threshold (usually 1.0, meaning >= 2x change)
        """
        sig_mask = (
            (results["padj"]  < self.padj_threshold) &
            (results["log2FC"].abs() > self.lfc_threshold)
        )
        results["significant"] = sig_mask
        results["regulation"]  = "not_significant"
        results.loc[sig_mask & (results["log2FC"] > 0), "regulation"] = "up"
        results.loc[sig_mask & (results["log2FC"] < 0), "regulation"] = "down"

        # Sort by adjusted p-value (most significant first)
        return results.sort_values("padj")

    # ── Output ──────────────────────────────────────────────────────────────

    def _save_results(self, results: pd.DataFrame):
        """Save full results, significant genes only, and a JSON summary."""
        # Full table (all genes)
        full_path = self.output_dir / "de_results_full.csv"
        results.to_csv(full_path)

        # Significant genes only
        sig      = results[results["significant"]]
        sig_path = self.output_dir / "de_results_significant.csv"
        sig.to_csv(sig_path)

        # JSON summary
        summary = {
            "total_genes_tested": int(len(results)),
            "significant_degs":   int(len(sig)),
            "upregulated":        int((results["regulation"] == "up").sum()),
            "downregulated":      int((results["regulation"] == "down").sum()),
            "padj_threshold":     self.padj_threshold,
            "lfc_threshold":      self.lfc_threshold,
            "test_method":        self.test_method,
            "top_10_upregulated":   (
                sig[sig["regulation"] == "up"].head(10).index.tolist()
            ),
            "top_10_downregulated": (
                sig[sig["regulation"] == "down"].head(10).index.tolist()
            ),
        }
        json_path = self.output_dir / "de_summary.json"
        with open(json_path, "w") as f:
            json.dump(summary, f, indent=2)

        logger.info(f"  Full results  -> {full_path}  ({len(results):,} genes)")
        logger.info(f"  Significant   -> {sig_path}  ({len(sig):,} DEGs)")
        logger.info(f"  Summary JSON  -> {json_path}")
