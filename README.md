# 🧬 RNAseq Analysis Pipeline
### Normal Tissue vs Cancer Tumor — Gene Expression & Silencing

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Beginner Friendly](https://img.shields.io/badge/level-beginner-green.svg)]()
[![Open Source](https://img.shields.io/badge/open-source-brightgreen.svg)]()

A beginner-friendly, open-source bioinformatics pipeline for analyzing RNA sequencing data.
Identifies genes that are **overexpressed** (potential oncogenes) or **silenced** (potential
tumor suppressors) by comparing normal tissue samples against cancerous tumor samples.

---

## 🗂 Project Structure

```
rnaseq_pipeline/
│
├── pipeline.py                     ← Entry point. Run this to start the analysis.
│
├── modules/                        ← Core analysis modules (one task each)
│   ├── __init__.py
│   ├── db_connector.py             ← Downloads data from NCBI GEO database
│   ├── qc.py                       ← Quality control & sample filtering
│   ├── normalization.py            ← CPM / TPM / TMM / DESeq2
│   ├── differential_expression.py  ← Statistical testing (t-test + FDR)
│   ├── visualizer.py               ← All plots (PCA, Volcano, Heatmap, MA, Boxplots)
│   └── reporter.py                 ← Writes final Markdown + JSON report
│
├── utils/                          ← Shared helper utilities
│   ├── __init__.py
│   ├── logger.py                   ← Console + file logging
│   └── config_loader.py            ← Loads pipeline_config.json
│
├── config/
│   └── pipeline_config.json        ← ⚙️  Edit this to change analysis parameters
│
├── data/
│   ├── example_counts.csv          ← Sample count matrix (use as template)
│   ├── example_metadata.csv        ← Sample metadata (use as template)
│   └── cache/                      ← Auto-created when downloading from GEO
│
├── results/                        ← Auto-created when you run the pipeline
│   ├── analysis_report.md
│   ├── analysis_summary.json
│   ├── qc/
│   ├── de_results/
│   └── plots/
│
└── requirements.txt
```

---

## ⚡ Quick Start (5 minutes)

### Step 1 — Clone the repository

```bash
git clone https://github.com/your-username/rnaseq-pipeline.git
cd rnaseq-pipeline
```

### Step 2 — Create a virtual environment

```bash
# Create the environment
python -m venv venv

# Activate it
source venv/bin/activate        # macOS / Linux
venv\Scripts\activate           # Windows
```

> **What is a virtual environment?**
> It's an isolated Python installation just for this project.
> This prevents conflicts with other Python packages on your machine.

### Step 3 — Install dependencies

```bash
pip install -r requirements.txt
```

### Step 4 — Run the pipeline

```bash
# Option A: Download a real cancer dataset from GEO (requires internet)
python pipeline.py --geo-id GSE47756 --output results/folder_name

# Option B: Use the included example data (no internet needed)
python pipeline.py --counts data/example_counts.csv --metadata data/example_metadata.csv --output results/folder_name
```

That's it! Check the `results/` folder for all outputs. Adding a different name for each folder in "folder_name" is recommendable to avoid overwrite your results.

---

## 💾 How to Load Your Own Data

### Your count matrix (`counts.csv`)

This file contains the raw sequencing read counts for each gene in each sample.

```
gene_id,Normal_01,Normal_02,Normal_03,Tumor_01,Tumor_02,Tumor_03
TP53,1200,1150,980,320,280,410
MYC,450,520,480,2100,1980,2300
BRCA1,890,920,850,210,195,240
...
```

Rules:
- First column: gene identifiers (HGNC symbols, Ensembl IDs, or any unique name)
- Column headers: sample names (must exactly match the `sample_id` column in metadata)
- Values: **integer raw read counts** (NOT RPKM, not TPM — the pipeline normalizes for you)

### Your metadata file (`metadata.csv`)

This file describes each sample.

```
sample_id,condition,tissue,batch
Normal_01,normal,breast,batch1
Normal_02,normal,breast,batch1
Tumor_01,tumor,breast,batch1
Tumor_02,tumor,breast,batch2
```

Rules:
- `sample_id`: must exactly match the column headers in your counts file
- `condition`: must be exactly **`normal`** or **`tumor`** (lowercase)
- Other columns (`tissue`, `batch`, etc.) are optional but recommended

### Run with your files

```bash
python pipeline.py \
    --counts path/to/your_counts.csv \
    --metadata path/to/your_metadata.csv \
    --output results/
```

---

## 🌐 How to Download Data from GEO

GEO (Gene Expression Omnibus) is NCBI's public database with thousands of
published RNA-seq experiments. This pipeline connects to it automatically.

### Find a dataset

1. Go to https://www.ncbi.nlm.nih.gov/geo/
2. Search for, e.g., `breast cancer RNA-seq normal tumor`
3. Note the accession ID (format: `GSExxxxx`)

### Download and analyze it

```bash
python pipeline.py --geo-id GSE47756 --output results/
```

The pipeline will:
1. Query the NCBI API for dataset info
2. Download the expression matrix file (~50-200 MB)
3. Cache it locally in `data/cache/` so you don't re-download
4. Auto-detect which samples are Normal and which are Tumor

**Note:** If you are offline or the download fails, the pipeline automatically
falls back to a synthetic demo dataset so you can still explore all features.

### Use cached data (offline / faster)

```bash
python pipeline.py --geo-id GSE47756 --use-cache --output results/
```

### Recommended datasets to start with

| GEO ID | Cancer Type | Samples | Notes |
|--------|-------------|---------|-------|
| `GSE47756` | Breast | 12 | Great for beginners |
| `GSE19804` | Lung | 60 | Lung adenocarcinoma |
| `GSE14520` | Liver | 247 | Hepatocellular carcinoma |
| `GSE22093` | Colorectal | 290 | Colorectal cancer |

---

## ⚙️ Configuration Reference

Edit `config/pipeline_config.json` to change analysis parameters:

```json
{
  "qc": {
    "min_count":   10,   <- Genes with fewer reads than this are filtered out
    "min_samples":  3,   <- Gene must have min_count in at least this many samples
    "outlier_sd":   3    <- Flag samples > N standard deviations from group mean
  },
  "normalization": {
    "method": "cpm"      <- Options: "cpm" | "tpm" | "tmm" | "deseq2"
  },
  "differential_expression": {
    "test_method":    "welch",  <- "welch" or "mannwhitney"
    "padj_threshold":  0.05,   <- Max adjusted p-value for significance
    "lfc_threshold":   1.0     <- Min |log2 fold change| for significance
  }
}
```

### Normalization methods explained

| Method | Best for | Corrects for |
|--------|----------|--------------|
| **CPM** | Simple comparisons, getting started | Library size |
| **TPM** | Genes of very different lengths | Library size + gene length |
| **TMM** | Many samples, composition bias | Library size + composition |
| **DESeq2** | Standard bulk RNA-seq (gold standard) | Library size + composition |

---

## 📊 Understanding the Outputs

### Plots

| File | What it shows | What to look for |
|------|---------------|------------------|
| `qc/qc_summary.png` | 4-panel QC dashboard | Even library sizes, correlated samples |
| `plots/pca_plot.png` | Sample clustering | Normal and Tumor clusters should separate |
| `plots/volcano_plot.png` | All genes: effect vs significance | Points in top-left / top-right corners |
| `plots/heatmap.png` | Top DEG expression across samples | Clear blue/red blocks per condition |
| `plots/ma_plot.png` | Normalization check | Cloud centered at M = 0 |
| `plots/gene_expression_boxplots.png` | Per-gene Normal vs Tumor | Separated boxes = real effect |

### CSV result files

**`de_results/de_results_full.csv`** — one row per gene:

| Column | Meaning |
|--------|---------|
| `gene` | Gene symbol |
| `log2FC` | log2 Fold Change (positive = higher in tumor) |
| `mean_normal` | Average expression in normal samples |
| `mean_tumor` | Average expression in tumor samples |
| `pvalue` | Raw p-value from statistical test |
| `padj` | FDR-adjusted p-value (use this for decisions) |
| `significant` | True if padj < threshold AND |log2FC| > threshold |
| `regulation` | `up` / `down` / `not_significant` |

**`de_results/de_results_significant.csv`** — same columns, filtered to significant genes only.

---

## 🖥️ Full CLI Reference

```bash
python pipeline.py [OPTIONS]

Options:
  --geo-id      TEXT    GEO accession ID to download (e.g. GSE47756)
  --counts      PATH    Path to your local count matrix CSV
  --metadata    PATH    Path to your local metadata CSV
  --config      PATH    Path to config JSON [default: config/pipeline_config.json]
  --output      PATH    Output directory [default: results/]
  --use-cache           Load from local cache instead of re-downloading
  --skip-qc             Skip the quality control step
  --log-level   LEVEL   DEBUG | INFO | WARNING | ERROR [default: INFO]
  --help                Show this help message
```

---

## 🧬 Bioinformatics Concepts Explained

### What is RNA-seq?
RNA sequencing measures all active gene transcripts in a cell at a given moment.
The result is a count matrix: how many RNA molecules were detected for each gene
in each sample. Higher counts = more gene activity.

### What is differential expression?
Comparing counts between two conditions (Normal vs Tumor) to find genes that
behave differently. A gene is **upregulated** if it's more active in tumor
(possible oncogene) or **downregulated** if less active (possible tumor suppressor).

### What is the log2 fold change?
Measures how many times more (or less) a gene is expressed in tumor vs normal:
- **LFC = +2** → 4× higher expression in tumor (2² = 4)
- **LFC = -1** → 2× lower expression in tumor (silenced)
- **LFC = 0** → no change

### What is FDR / adjusted p-value?
When testing 5,000 genes at p < 0.05, about 250 will appear significant by
chance (false positives). The **False Discovery Rate (FDR)** correction controls
the proportion of false positives among all called significant genes.
Always use `padj`, not `pvalue`, to make decisions.

---

## 🔧 Extending the Pipeline

### Add a new normalization method

```python
# In modules/normalization.py:
def _my_method(self, counts_df: pd.DataFrame) -> pd.DataFrame:
    # Your implementation
    return normalized_df
```

Then add `"my_method"` to `SUPPORTED_METHODS` and add an `elif` branch in `normalize()`.

### Add a new visualization

```python
# In modules/visualizer.py:
def plot_my_figure(self, norm_df, metadata_df):
    fig, ax = plt.subplots(figsize=(10, 7))
    # ... plotting code ...
    out = self.output_dir / "my_figure.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  My figure -> {out}")
```

Then call it in `pipeline.py` inside the visualizations step.

### Connect to another database (e.g., TCGA)

```python
# Create modules/tcga_connector.py following the same interface:
class TCGAConnector:
    def fetch(self, project_id: str, use_cache: bool = True):
        # Must return: (counts_df, metadata_df)
        ...
```

---

## 📚 Learning Resources

- 📖 [StatQuest: RNA-seq explained](https://www.youtube.com/watch?v=tlf6wYJrwKY) — best intro video
- 📖 [StatQuest: DESeq2](https://www.youtube.com/watch?v=UFB993xufUU) — differential expression
- 📖 [NCBI GEO Database Guide](https://www.ncbi.nlm.nih.gov/geo/info/overview.html)
- 📖 [Bioconductor RNA-seq Workflow](https://bioconductor.org/packages/release/workflows/html/rnaseqGene.html)
- 📖 [Harvard HBC RNA-seq Training](https://hbctraining.github.io/DGE_workshop/)

---

## 🤝 Contributing

Contributions are welcome! Ideas:

- [ ] TCGA (The Cancer Genome Atlas) database connector
- [ ] Gene Ontology (GO) / KEGG pathway enrichment analysis
- [ ] Streamlit web interface for interactive exploration
- [ ] Batch effect correction (ComBat)
- [ ] Single-cell RNA-seq support

**Steps:**
1. Fork the repository
2. Create a branch: `git checkout -b feature/my-feature`
3. Commit: `git commit -m "Add: description"`
4. Push & open a Pull Request

---

## 📝 License

MIT License — free to use, modify, and distribute.

---
*Built for the open-source bioinformatics community* 🧬
