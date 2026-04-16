"""
modules/qc.py
=============
Quality Control (QC) for RNAseq count matrices.

Steps performed:
    1. Low-expression gene filtering
       Remove genes that are not detected (or barely detected) across samples.
       These genes carry no useful signal and add noise to the analysis.

    2. Zero-variance gene filtering
       Remove genes with identical expression across all samples.
       These cannot be differentially expressed by definition.

    3. Outlier sample detection
       Flag samples whose average expression deviates greatly from others.
       Outliers can indicate failed sequencing or mislabeled samples.

    4. QC visualization panel (4 subplots saved as PNG)

    5. QC metrics export (JSON)
"""

import json
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)


class QualityControl:
    """
    Quality control and filtering for RNAseq count matrices.

    Args:
        output_dir  (Path): Directory where QC plots and metrics are saved.
        config      (dict): Pipeline configuration dictionary.
    """

    def __init__(self, output_dir="results/qc/", config: dict = None):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.config   = config or {}
        qc_cfg        = self.config.get("qc", {})

        # Thresholds (all configurable in pipeline_config.json)
        self.min_count   = qc_cfg.get("min_count",   10)  # minimum reads to keep a gene
        self.min_samples = qc_cfg.get("min_samples",  3)   # in how many samples
        self.outlier_sd  = qc_cfg.get("outlier_sd",   3)   # standard deviations threshold

        self.metrics = {}  # Collected during run(), exported as JSON

    # ── Public API ──────────────────────────────────────────────────────────

    def run(self, counts_df: pd.DataFrame, metadata_df: pd.DataFrame):
        """
        Run the full QC pipeline.

        Args:
            counts_df   (pd.DataFrame): Raw count matrix (genes x samples).
            metadata_df (pd.DataFrame): Sample metadata with 'condition' column.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: Filtered count matrix and metadata.
        """
        logger.info(f"  Input : {counts_df.shape[0]:,} genes x {counts_df.shape[1]} samples")

        counts_df = self._filter_low_expression(counts_df)
        counts_df = self._filter_zero_variance(counts_df)
        counts_df, metadata_df = self._detect_outliers(counts_df, metadata_df)

        self._compute_metrics(counts_df)
        self._plot_qc_panel(counts_df, metadata_df)
        self._save_metrics()

        logger.info(f"  Output: {counts_df.shape[0]:,} genes x {counts_df.shape[1]} samples")
        return counts_df, metadata_df

    # ── Filtering steps ─────────────────────────────────────────────────────

    def _filter_low_expression(self, counts_df: pd.DataFrame) -> pd.DataFrame:
        """
        Keep only genes with >= min_count reads in at least min_samples samples.

        Example: min_count=10, min_samples=3 means:
            Keep a gene if at least 3 samples have 10+ reads for that gene.
        """
        mask      = (counts_df >= self.min_count).sum(axis=1) >= self.min_samples
        filtered  = counts_df[mask]
        n_removed = len(counts_df) - len(filtered)
        logger.info(f"  Low-expression filter : removed {n_removed:,} genes "
                    f"(< {self.min_count} counts in < {self.min_samples} samples)")
        self.metrics["genes_removed_low_expr"] = int(n_removed)
        return filtered

    def _filter_zero_variance(self, counts_df: pd.DataFrame) -> pd.DataFrame:
        """Remove genes with zero variance (same value across all samples)."""
        var       = counts_df.var(axis=1)
        filtered  = counts_df[var > 0]
        n_removed = len(counts_df) - len(filtered)
        logger.info(f"  Zero-variance filter  : removed {n_removed:,} genes")
        self.metrics["genes_removed_zero_var"] = int(n_removed)
        return filtered

    def _detect_outliers(
        self,
        counts_df: pd.DataFrame,
        metadata_df: pd.DataFrame
    ):
        """
        Detect outlier samples based on log2 mean expression deviation.

        Samples whose mean expression deviates more than outlier_sd standard
        deviations from the group mean are flagged as potential outliers.
        Outliers are NOT automatically removed — they are logged for manual review.
        """
        log_counts    = np.log2(counts_df + 1)
        sample_means  = log_counts.mean(axis=0)
        overall_mean  = sample_means.mean()
        overall_std   = sample_means.std()

        z_scores = (sample_means - overall_mean) / (overall_std + 1e-9)
        outliers = z_scores[z_scores.abs() > self.outlier_sd].index.tolist()

        self.metrics["potential_outlier_samples"] = outliers
        if outliers:
            logger.warning(f"  Potential outlier samples detected: {outliers}")
            logger.warning("  These samples were kept — review manually before excluding.")
        else:
            logger.info("  No outlier samples detected.")

        return counts_df, metadata_df

    # ── Metrics ─────────────────────────────────────────────────────────────

    def _compute_metrics(self, counts_df: pd.DataFrame):
        """Compute descriptive metrics for the filtered matrix."""
        total_counts = counts_df.sum(axis=0)
        self.metrics.update({
            "final_genes":         int(counts_df.shape[0]),
            "final_samples":       int(counts_df.shape[1]),
            "mean_library_size":   float(total_counts.mean()),
            "median_library_size": float(total_counts.median()),
            "min_library_size":    float(total_counts.min()),
            "max_library_size":    float(total_counts.max()),
        })

    def _save_metrics(self):
        out = self.output_dir / "qc_metrics.json"
        with open(out, "w") as f:
            json.dump(self.metrics, f, indent=2, default=str)
        logger.info(f"  QC metrics  -> {out}")

    # ── Visualization ────────────────────────────────────────────────────────

    def _plot_qc_panel(self, counts_df: pd.DataFrame, metadata_df: pd.DataFrame):
        """
        Generate a 2x2 QC summary panel:
            [0,0] Library sizes bar chart
            [0,1] Log2 expression distribution boxplots
            [1,0] Expression density overlays
            [1,1] Sample-to-sample correlation heatmap
        """
        fig = plt.figure(figsize=(16, 12))
        fig.suptitle(
            "Quality Control Report  —  RNAseq Pipeline",
            fontsize=16, fontweight="bold"
        )
        gs         = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.35)
        colors     = self._condition_colors(metadata_df, counts_df.columns)
        log_counts = np.log2(counts_df + 1)

        # ── Panel 1: Library sizes ──
        ax1 = fig.add_subplot(gs[0, 0])
        lib_sizes = counts_df.sum(axis=0) / 1e6
        ax1.bar(range(len(lib_sizes)), lib_sizes, color=colors, edgecolor="white")
        ax1.axhline(lib_sizes.mean(), color="crimson", linestyle="--", linewidth=1.5,
                    label=f"Mean: {lib_sizes.mean():.1f}M")
        ax1.set_title("Library Sizes", fontweight="bold")
        ax1.set_xlabel("Sample")
        ax1.set_ylabel("Total Counts (Millions)")
        ax1.set_xticks(range(len(lib_sizes)))
        ax1.set_xticklabels(lib_sizes.index, rotation=45, ha="right", fontsize=7)
        ax1.legend(fontsize=9)
        ax1.set_facecolor("#F8F9FA")

        # ── Panel 2: Log2 expression boxplots ──
        ax2 = fig.add_subplot(gs[0, 1])
        bplot = ax2.boxplot(
            log_counts.values,
            patch_artist=True,
            notch=False,
            medianprops=dict(color="black", linewidth=2),
            flierprops=dict(marker="o", markersize=2, alpha=0.3),
        )
        for patch, color in zip(bplot["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
        ax2.set_title("Log\u2082 Expression Distribution", fontweight="bold")
        ax2.set_xlabel("Sample")
        ax2.set_ylabel("log\u2082(count + 1)")
        ax2.set_xticks(range(1, len(log_counts.columns) + 1))
        ax2.set_xticklabels(log_counts.columns, rotation=45, ha="right", fontsize=7)
        ax2.set_facecolor("#F8F9FA")

        # ── Panel 3: Expression density ──
        ax3 = fig.add_subplot(gs[1, 0])
        for i, col in enumerate(log_counts.columns):
            ax3.hist(
                log_counts[col].values,
                bins=60, density=True, alpha=0.35, color=colors[i]
            )
        ax3.set_title("Log\u2082 Expression Density", fontweight="bold")
        ax3.set_xlabel("log\u2082(count + 1)")
        ax3.set_ylabel("Density")
        ax3.set_facecolor("#F8F9FA")
        from matplotlib.patches import Patch
        ax3.legend(handles=[
            Patch(facecolor="#2E86AB", alpha=0.7, label="Normal"),
            Patch(facecolor="#E84855", alpha=0.7, label="Tumor"),
        ])

        # ── Panel 4: Sample correlation heatmap ──
        ax4 = fig.add_subplot(gs[1, 1])
        corr = log_counts.corr()
        sns.heatmap(
            corr, ax=ax4, cmap="RdYlGn", vmin=0.8, vmax=1.0,
            annot=(len(corr) <= 12), fmt=".2f",
            xticklabels=corr.columns, yticklabels=corr.index,
            cbar_kws={"shrink": 0.8, "label": "Pearson r"},
            linewidths=0.2, linecolor="#EEEEEE",
        )
        ax4.set_title("Sample-to-Sample Correlation", fontweight="bold")
        ax4.set_xticklabels(ax4.get_xticklabels(), rotation=45, ha="right", fontsize=7)
        ax4.set_yticklabels(ax4.get_yticklabels(), rotation=0, fontsize=7)

        out = self.output_dir / "qc_summary.png"
        fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        logger.info(f"  QC plot     -> {out}")

    def _condition_colors(self, metadata_df: pd.DataFrame, sample_names) -> list:
        """Map sample names to colors based on their condition label."""
        palette = {"normal": "#2E86AB", "tumor": "#E84855"}
        colors  = []
        for s in sample_names:
            if s in metadata_df.index and "condition" in metadata_df.columns:
                cond = str(metadata_df.loc[s, "condition"]).lower()
                colors.append(palette.get(cond, "#999999"))
            else:
                colors.append("#999999")
        return colors
