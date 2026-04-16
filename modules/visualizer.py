"""
modules/visualizer.py
=====================
Visualization module for RNAseq differential expression analysis.

Plots generated:
    1. PCA Plot
       Principal Component Analysis of samples.
       Reduces thousands of gene dimensions to 2-3 axes.
       Good separation of Normal vs Tumor clusters = clean experiment.

    2. Volcano Plot
       X-axis: log2 Fold Change (effect size).
       Y-axis: -log10(adjusted p-value) (statistical confidence).
       Top-right = strongly upregulated AND significant.
       Top-left  = strongly downregulated AND significant.

    3. Heatmap
       Expression of the top differentially expressed genes across all samples.
       Color = z-score (how far above/below average is this gene in this sample).
       Hierarchical clustering groups genes and samples with similar patterns.

    4. MA Plot
       X-axis (A) = average expression of the gene.
       Y-axis (M) = log2 fold change.
       Ideal: cloud of points centered at M=0.
       Systematic deviation = normalization problem.

    5. Gene Expression Boxplots
       Distribution of expression for each top gene, split by condition.
       Individual data points shown as jittered dots.
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

logger = logging.getLogger(__name__)

# ── Global color palette ────────────────────────────────────────────────────────
PALETTE = {
    "normal":    "#2E86AB",  # Blue  — normal samples
    "tumor":     "#E84855",  # Red   — tumor samples
    "up":        "#E84855",  # Red   — upregulated in tumor
    "down":      "#2E86AB",  # Blue  — downregulated in tumor
    "ns":        "#BBBBBB",  # Gray  — not significant
    "highlight": "#F4A261",  # Orange — annotation highlights
    "bg":        "#F8F9FA",  # Light background for axes
}

# Global matplotlib style
plt.rcParams.update({
    "font.family":        "DejaVu Sans",
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "figure.dpi":         100,
})


class Visualizer:
    """
    Generate all publication-quality visualizations for the RNAseq pipeline.

    Args:
        output_dir (str | Path): Directory where PNG files are saved.
    """

    def __init__(self, output_dir="results/plots/"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    # ── 1. PCA Plot ─────────────────────────────────────────────────────────

    def plot_pca(self, norm_df: pd.DataFrame, metadata_df: pd.DataFrame):
        """
        PCA scatter plot + scree plot.

        The PCA is computed on all genes simultaneously. If Normal and Tumor
        samples cluster separately on PC1 or PC2, the gene expression profiles
        are genuinely different — a good sign for the analysis.
        """
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle(
            "Principal Component Analysis (PCA)  —  Sample Overview",
            fontsize=14, fontweight="bold"
        )

        # Prepare: transpose so rows=samples, scale features
        X        = norm_df.T.values
        scaler   = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        n_comp   = min(5, X.shape[0], X.shape[1])
        pca      = PCA(n_components=n_comp)
        comps    = pca.fit_transform(X_scaled)
        var_exp  = pca.explained_variance_ratio_ * 100

        samples    = norm_df.columns.tolist()
        colors     = []
        conditions = []
        for s in samples:
            if s in metadata_df.index and "condition" in metadata_df.columns:
                c = str(metadata_df.loc[s, "condition"]).lower()
            else:
                c = "unknown"
            conditions.append(c)
            colors.append(PALETTE.get(c, "#888888"))

        # Left panel: PC1 vs PC2
        ax1 = axes[0]
        for i, (sample, cond, color) in enumerate(zip(samples, conditions, colors)):
            ax1.scatter(
                comps[i, 0], comps[i, 1],
                c=color, s=120, edgecolors="white", linewidth=1.5, zorder=3
            )
            ax1.annotate(
                sample,
                (comps[i, 0], comps[i, 1]),
                textcoords="offset points", xytext=(6, 4),
                fontsize=7, alpha=0.85
            )
        ax1.set_xlabel(f"PC1  ({var_exp[0]:.1f}% variance)", fontsize=11)
        ax1.set_ylabel(f"PC2  ({var_exp[1]:.1f}% variance)", fontsize=11)
        ax1.set_title("Sample Clustering — PC1 vs PC2", fontsize=11)
        ax1.set_facecolor(PALETTE["bg"])
        ax1.grid(True, alpha=0.3, linewidth=0.6)
        ax1.legend(
            handles=[
                mpatches.Patch(color=PALETTE["normal"], label="Normal"),
                mpatches.Patch(color=PALETTE["tumor"],  label="Tumor"),
            ],
            framealpha=0.9, fontsize=10
        )

        # Right panel: Scree plot
        ax2 = axes[1]
        x_labels = [f"PC{i+1}" for i in range(n_comp)]
        bars     = ax2.bar(
            x_labels, var_exp[:n_comp],
            color="#5C7AEA", edgecolor="white", linewidth=1.5
        )
        ax2.plot(
            x_labels, np.cumsum(var_exp[:n_comp]),
            "o-", color=PALETTE["highlight"], linewidth=2.2,
            label="Cumulative %"
        )
        for bar, val in zip(bars, var_exp[:n_comp]):
            ax2.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.5,
                f"{val:.1f}%",
                ha="center", va="bottom", fontsize=9
            )
        ax2.set_xlabel("Principal Component", fontsize=11)
        ax2.set_ylabel("Variance Explained (%)", fontsize=11)
        ax2.set_title("Scree Plot", fontsize=11)
        ax2.set_facecolor(PALETTE["bg"])
        ax2.grid(True, alpha=0.3, axis="y", linewidth=0.6)
        ax2.legend(fontsize=10)

        fig.tight_layout()
        out = self.output_dir / "pca_plot.png"
        fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        logger.info(f"  PCA plot       -> {out}")

    # ── 2. Volcano Plot ─────────────────────────────────────────────────────

    def plot_volcano(self, results_df: pd.DataFrame, config: dict = None):
        """
        Volcano plot: log2FoldChange vs -log10(adjusted p-value).

        Points far right & high  = upregulated oncogenes (tumor > normal).
        Points far left & high   = silenced tumor suppressors (tumor < normal).
        The top 15 most significant genes are labeled automatically.
        """
        cfg      = (config or {}).get("differential_expression", {})
        padj_thr = cfg.get("padj_threshold", 0.05)
        lfc_thr  = cfg.get("lfc_threshold",  1.0)

        df              = results_df.copy()
        df["-log10padj"] = -np.log10(df["padj"].clip(1e-300))

        color_map = {
            "up":              PALETTE["up"],
            "down":            PALETTE["down"],
            "not_significant": PALETTE["ns"],
        }
        colors = df["regulation"].map(color_map).fillna(PALETTE["ns"])
        sizes  = np.where(df["significant"], 55, 12)

        fig, ax = plt.subplots(figsize=(11, 8))
        ax.scatter(
            df["log2FC"], df["-log10padj"],
            c=colors, s=sizes, alpha=0.72, edgecolors="none", zorder=2
        )

        # Significance threshold lines
        ax.axhline(
            -np.log10(padj_thr), color="#666666",
            linestyle="--", linewidth=1, alpha=0.75,
            label=f"padj = {padj_thr}"
        )
        ax.axvline(
             lfc_thr, color="#666666", linestyle="--", linewidth=1, alpha=0.75)
        ax.axvline(
            -lfc_thr, color="#666666", linestyle="--", linewidth=1, alpha=0.75)

        # Annotate top 15 significant genes
        top = df[df["significant"]].nsmallest(15, "padj")
        for gene, row in top.iterrows():
            ax.annotate(
                gene,
                (row["log2FC"], row["-log10padj"]),
                textcoords="offset points", xytext=(5, 4),
                fontsize=7.5, fontweight="bold", color="#222222",
                arrowprops=dict(arrowstyle="-", color="#AAAAAA", lw=0.7),
            )

        n_up   = (df["regulation"] == "up").sum()
        n_down = (df["regulation"] == "down").sum()

        ax.set_xlabel("log\u2082 Fold Change  (Tumor / Normal)", fontsize=12)
        ax.set_ylabel("\u2212log\u2081\u2080 (adjusted p-value)", fontsize=12)
        ax.set_title(
            "Volcano Plot  —  Differential Gene Expression\n"
            "Normal vs Tumor",
            fontsize=14, fontweight="bold"
        )
        ax.set_facecolor(PALETTE["bg"])
        ax.grid(True, alpha=0.2, linewidth=0.6)
        ax.legend(
            handles=[
                mpatches.Patch(color=PALETTE["up"],   label=f"Upregulated in Tumor  (n={n_up})"),
                mpatches.Patch(color=PALETTE["down"], label=f"Downregulated in Tumor (n={n_down})"),
                mpatches.Patch(color=PALETTE["ns"],   label="Not Significant"),
            ],
            loc="upper left", framealpha=0.9, fontsize=10
        )

        fig.tight_layout()
        out = self.output_dir / "volcano_plot.png"
        fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        logger.info(f"  Volcano plot   -> {out}")

    # ── 3. Heatmap ──────────────────────────────────────────────────────────

    def plot_heatmap(
        self,
        norm_df: pd.DataFrame,
        results_df: pd.DataFrame,
        metadata_df: pd.DataFrame,
        config: dict = None,
        top_n: int = 50
    ):
        """
        Clustered heatmap of the top differentially expressed genes.

        Expression values are z-score normalized per gene so that colors
        reflect relative (not absolute) expression levels:
            Blue  = below average expression for this gene
            White = average
            Red   = above average expression for this gene
        """
        sig_genes = results_df[results_df["significant"]].index.tolist()
        if len(sig_genes) == 0:
            logger.warning("  No significant genes — showing top 50 by p-value instead.")
            sig_genes = results_df.nsmallest(50, "padj").index.tolist()

        genes     = sig_genes[:top_n]
        plot_data = norm_df.loc[genes]

        # Z-score normalization across samples (row-wise)
        z_data = plot_data.subtract(plot_data.mean(axis=1), axis=0)
        z_data = z_data.divide(plot_data.std(axis=1).clip(1e-9), axis=0)

        # Column color bar by condition
        cond_colors = {}
        for col in z_data.columns:
            cond = "unknown"
            if col in metadata_df.index and "condition" in metadata_df.columns:
                cond = str(metadata_df.loc[col, "condition"]).lower()
            cond_colors[col] = PALETTE.get(cond, "#888888")
        col_color_series = pd.Series(cond_colors)

        # Blue-White-Red colormap
        cmap = LinearSegmentedColormap.from_list(
            "rnaseq_bwr", ["#2166AC", "#F7F7F7", "#D6604D"], N=256
        )

        fig_h = max(10, len(genes) * 0.22)
        g = sns.clustermap(
            z_data,
            cmap=cmap,
            center=0, vmin=-3, vmax=3,
            col_colors=col_color_series,
            figsize=(14, fig_h),
            xticklabels=True, yticklabels=True,
            dendrogram_ratio=(0.10, 0.12),
            cbar_pos=(0.02, 0.80, 0.03, 0.15),
            linewidths=0,
        )
        g.ax_heatmap.set_xticklabels(
            g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=8
        )
        g.ax_heatmap.set_yticklabels(
            g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=7
        )
        g.fig.suptitle(
            f"Heatmap  —  Top {len(genes)} Differentially Expressed Genes\n"
            "Z-score normalized expression (row-wise)",
            fontsize=13, fontweight="bold", y=1.01
        )
        g.ax_col_colors.legend(
            handles=[
                mpatches.Patch(color=PALETTE["normal"], label="Normal"),
                mpatches.Patch(color=PALETTE["tumor"],  label="Tumor"),
            ],
            loc="center right", bbox_to_anchor=(1.22, 0.5),
            framealpha=0.9, fontsize=9
        )

        out = self.output_dir / "heatmap.png"
        g.fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(g.fig)
        logger.info(f"  Heatmap        -> {out}")

    # ── 4. MA Plot ──────────────────────────────────────────────────────────

    def plot_ma(self, results_df: pd.DataFrame):
        """
        MA plot: Average expression (A) vs log2 Fold Change (M).

        'MA' stands for Minus/Average. It was introduced in microarray analysis
        and remains a standard diagnostic for RNA-seq normalization quality.

        A well-normalized dataset shows a symmetric cloud centered at M = 0
        for genes with moderate average expression. Systematic up- or
        down-shifts indicate normalization problems.
        """
        df      = results_df.copy()
        df["A"] = (df["mean_normal"] + df["mean_tumor"]) / 2.0  # average expression

        color_map = {
            "up":              PALETTE["up"],
            "down":            PALETTE["down"],
            "not_significant": PALETTE["ns"],
        }
        colors = df["regulation"].map(color_map).fillna(PALETTE["ns"])
        sizes  = np.where(df["significant"], 45, 9)

        fig, ax = plt.subplots(figsize=(10, 7))
        ax.scatter(df["A"], df["log2FC"], c=colors, s=sizes, alpha=0.65, edgecolors="none")

        # Reference lines
        ax.axhline(0,    color="black",  linewidth=1.8, linestyle="-",  zorder=3)
        ax.axhline( 1.0, color="#666666", linewidth=0.9, linestyle="--", alpha=0.7)
        ax.axhline(-1.0, color="#666666", linewidth=0.9, linestyle="--", alpha=0.7)

        n_up   = (df["regulation"] == "up").sum()
        n_down = (df["regulation"] == "down").sum()

        ax.set_xlabel("A  —  Average log\u2082 Expression", fontsize=12)
        ax.set_ylabel("M  —  log\u2082 Fold Change (Tumor / Normal)", fontsize=12)
        ax.set_title(
            "MA Plot  —  Normalization Diagnostic\n"
            "Cloud centered at M=0 indicates good normalization",
            fontsize=13, fontweight="bold"
        )
        ax.set_facecolor(PALETTE["bg"])
        ax.grid(True, alpha=0.2, linewidth=0.6)
        ax.legend(
            handles=[
                mpatches.Patch(color=PALETTE["up"],   label=f"Upregulated ({n_up})"),
                mpatches.Patch(color=PALETTE["down"], label=f"Downregulated ({n_down})"),
                mpatches.Patch(color=PALETTE["ns"],   label="Not significant"),
            ],
            framealpha=0.9, fontsize=10
        )

        fig.tight_layout()
        out = self.output_dir / "ma_plot.png"
        fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        logger.info(f"  MA plot        -> {out}")

    # ── 5. Gene Expression Boxplots ─────────────────────────────────────────

    def plot_gene_expression_boxplots(
        self,
        norm_df: pd.DataFrame,
        metadata_df: pd.DataFrame,
        results_df: pd.DataFrame,
        config: dict = None,
        n_genes: int = 12
    ):
        """
        Individual expression boxplots for the most significant DEGs.

        Shows the full distribution of expression values across Normal and Tumor
        samples for each gene. Each dot is one biological sample.
        Box = interquartile range, line = median.
        """
        sig = results_df[results_df["significant"]]

        # Select top upregulated and top downregulated genes
        top_up   = sig[sig["regulation"] == "up"].head(n_genes // 2).index.tolist()
        top_down = sig[sig["regulation"] == "down"].head(n_genes // 2).index.tolist()
        genes    = (top_up + top_down)[:n_genes]

        if len(genes) == 0:
            genes = results_df.nsmallest(n_genes, "padj").index.tolist()

        n_cols = 4
        n_rows = int(np.ceil(len(genes) / n_cols))

        fig, axes = plt.subplots(
            n_rows, n_cols,
            figsize=(n_cols * 3.5, n_rows * 3.2)
        )
        fig.suptitle(
            "Top Differentially Expressed Genes  —  Normal vs Tumor",
            fontsize=14, fontweight="bold"
        )
        axes_flat = axes.flatten() if hasattr(axes, "flatten") else [axes]

        for i, gene in enumerate(genes):
            ax = axes_flat[i]

            # Build per-sample dataframe for this gene
            rows = []
            for col in norm_df.columns:
                if col in metadata_df.index:
                    cond = str(metadata_df.loc[col, "condition"]).lower()
                    expr = norm_df.loc[gene, col] if gene in norm_df.index else np.nan
                    rows.append({"sample": col, "condition": cond, "expression": expr})
            gene_df = pd.DataFrame(rows).dropna(subset=["expression"])

            normal_vals = gene_df[gene_df["condition"] == "normal"]["expression"].values
            tumor_vals  = gene_df[gene_df["condition"] == "tumor"]["expression"].values

            bp = ax.boxplot(
                [normal_vals, tumor_vals],
                patch_artist=True, widths=0.5, notch=False,
                medianprops=dict(color="black", linewidth=2.2),
                flierprops=dict(marker="o", markersize=2, alpha=0.3),
            )
            bp["boxes"][0].set(facecolor=PALETTE["normal"], alpha=0.75)
            if len(bp["boxes"]) > 1:
                bp["boxes"][1].set(facecolor=PALETTE["tumor"], alpha=0.75)

            # Overlay jittered individual data points
            for j, (vals, cond) in enumerate(
                [(normal_vals, "normal"), (tumor_vals, "tumor")]
            ):
                jitter = np.random.uniform(-0.12, 0.12, len(vals))
                ax.scatter(
                    j + 1 + jitter, vals,
                    color=PALETTE[cond], s=22, alpha=0.85, zorder=4
                )

            # Title with statistics
            if gene in results_df.index:
                padj = results_df.loc[gene, "padj"]
                lfc  = results_df.loc[gene, "log2FC"]
                direction = "↑ UP" if lfc > 0 else "↓ DOWN"
                ax.set_title(
                    f"{gene}   {direction}\n"
                    f"LFC={lfc:+.2f}  padj={padj:.1e}",
                    fontsize=8.5, fontweight="bold"
                )
            else:
                ax.set_title(gene, fontsize=9, fontweight="bold")

            ax.set_xticks([1, 2])
            ax.set_xticklabels(["Normal", "Tumor"], fontsize=9)
            ax.set_ylabel("log\u2082 Expression", fontsize=8)
            ax.set_facecolor(PALETTE["bg"])
            ax.grid(True, alpha=0.2, axis="y", linewidth=0.6)

        # Hide unused subplots
        for j in range(len(genes), len(axes_flat)):
            axes_flat[j].set_visible(False)

        fig.tight_layout()
        out = self.output_dir / "gene_expression_boxplots.png"
        fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        logger.info(f"  Boxplots       -> {out}")
