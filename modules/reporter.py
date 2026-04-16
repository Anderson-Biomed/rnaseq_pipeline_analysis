"""
modules/reporter.py
===================
Generate the final analysis report in Markdown and JSON formats.

Outputs:
    results/analysis_report.md   — Human-readable Markdown report.
    results/analysis_summary.json — Machine-readable JSON summary.
"""

import json
import logging
from datetime import datetime
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


class Reporter:
    """
    Compile all pipeline outputs into a final analysis report.

    Args:
        output_dir (str | Path): Root results directory.
    """

    def __init__(self, output_dir="results/"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate(
        self,
        counts_df: pd.DataFrame,
        norm_df: pd.DataFrame,
        metadata_df: pd.DataFrame,
        de_results: pd.DataFrame,
        config: dict,
        geo_id: str = None
    ):
        """Build and save the Markdown report and JSON summary."""
        sig        = de_results[de_results["significant"]]
        n_up       = int((de_results["regulation"] == "up").sum())
        n_down     = int((de_results["regulation"] == "down").sum())
        cond_cts   = (
            metadata_df["condition"].value_counts().to_dict()
            if "condition" in metadata_df.columns else {}
        )
        de_cfg     = config.get("differential_expression", {})
        norm_method = config.get("normalization", {}).get("method", "cpm").upper()

        # ── JSON Summary ───────────────────────────────────────────────────
        summary = {
            "pipeline":  "RNAseq Analysis Pipeline — Normal vs Tumor",
            "run_date":  datetime.now().isoformat(),
            "dataset":   geo_id or "local",
            "samples": {
                "total":  int(counts_df.shape[1]),
                "normal": int(cond_cts.get("normal", 0)),
                "tumor":  int(cond_cts.get("tumor",  0)),
            },
            "genes": {
                "after_qc":     int(counts_df.shape[0]),
                "tested":       int(len(de_results)),
                "significant":  int(len(sig)),
                "upregulated":  n_up,
                "downregulated": n_down,
            },
            "normalization": norm_method,
            "thresholds": {
                "padj":   de_cfg.get("padj_threshold", 0.05),
                "log2FC": de_cfg.get("lfc_threshold",  1.0),
            },
            "top_upregulated":   sig[sig["regulation"] == "up"].head(10).index.tolist(),
            "top_downregulated": sig[sig["regulation"] == "down"].head(10).index.tolist(),
        }
        json_path = self.output_dir / "analysis_summary.json"
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)

        # ── Markdown Report ────────────────────────────────────────────────
        top_up_table   = sig[sig["regulation"] == "up"].head(10)
        top_down_table = sig[sig["regulation"] == "down"].head(10)

        up_rows = "\n".join(
            f"| {rank} | **{gene}** | +{row['log2FC']:.2f} | {row['padj']:.2e} |"
            for rank, (gene, row) in enumerate(top_up_table.iterrows(), 1)
        )
        down_rows = "\n".join(
            f"| {rank} | **{gene}** | {row['log2FC']:.2f} | {row['padj']:.2e} |"
            for rank, (gene, row) in enumerate(top_down_table.iterrows(), 1)
        )

        md = f"""# RNAseq Analysis Report
## Normal Tissue vs Cancer Tumor — Differential Gene Expression

**Run date:** {datetime.now().strftime("%Y-%m-%d  %H:%M")}
**Dataset:** `{geo_id or "local data"}`
**Pipeline:** RNAseq-Pipeline v1.0

---

## Dataset Overview

| Parameter | Value |
|-----------|-------|
| Total samples | {counts_df.shape[1]} |
| Normal samples | {cond_cts.get("normal", "N/A")} |
| Tumor samples | {cond_cts.get("tumor", "N/A")} |
| Genes after QC | {counts_df.shape[0]:,} |
| Normalization method | {norm_method} |

---

## Differential Expression Results

| Metric | Value |
|--------|-------|
| Genes tested | {len(de_results):,} |
| **Significant DEGs** | **{len(sig):,}** |
| Upregulated in Tumor | {n_up} |
| Downregulated in Tumor | {n_down} |
| FDR threshold (padj) | {de_cfg.get("padj_threshold", 0.05)} |
| Fold-change threshold (\\|log2FC\\|) | {de_cfg.get("lfc_threshold", 1.0)} |

---

## Top Upregulated Genes (Tumor vs Normal)
> These genes show **increased** expression in tumor — candidate oncogenes.

| # | Gene | log2FC | padj |
|---|------|--------|------|
{up_rows}

---

## Top Downregulated Genes (Silenced in Tumor)
> These genes show **decreased** expression in tumor — candidate tumor suppressors.

| # | Gene | log2FC | padj |
|---|------|--------|------|
{down_rows}

---

## Generated Visualizations

| Plot | File | Description |
|------|------|-------------|
| QC Summary | `qc/qc_summary.png` | Library sizes, distributions, sample correlations |
| PCA | `plots/pca_plot.png` | Sample clustering in principal component space |
| Volcano | `plots/volcano_plot.png` | Statistical significance vs fold change |
| Heatmap | `plots/heatmap.png` | Top DEG expression across all samples |
| MA Plot | `plots/ma_plot.png` | Normalization diagnostic |
| Boxplots | `plots/gene_expression_boxplots.png` | Per-gene Normal vs Tumor expression |

---

## Output File Structure

```
results/
├── analysis_summary.json           <- Machine-readable summary
├── analysis_report.md              <- This report
├── qc/
│   ├── qc_summary.png
│   └── qc_metrics.json
├── de_results/
│   ├── de_results_full.csv         <- All genes with statistics
│   ├── de_results_significant.csv  <- Significant DEGs only
│   └── de_summary.json
└── plots/
    ├── pca_plot.png
    ├── volcano_plot.png
    ├── heatmap.png
    ├── ma_plot.png
    └── gene_expression_boxplots.png
```

---
*Generated by RNAseq Analysis Pipeline — Open Source Bioinformatics*
"""
        md_path = self.output_dir / "analysis_report.md"
        with open(md_path, "w", encoding="utf-8") as f:
            f.write(md)

        logger.info(f"  Markdown report -> {md_path}")
        logger.info(f"  JSON summary    -> {json_path}")
