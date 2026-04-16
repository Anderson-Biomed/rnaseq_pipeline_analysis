"""
modules/db_connector.py
=======================
Database connector for NCBI GEO (Gene Expression Omnibus).

This module handles:
    - Querying the NCBI eUtils API to find datasets by accession ID
    - Downloading expression matrix files (SOFT format) from GEO FTP
    - Parsing raw GEO matrix files into clean pandas DataFrames
    - Auto-detecting Normal vs Tumor conditions from sample metadata
    - Caching downloaded data to disk to avoid re-downloading
    - Generating realistic synthetic demo data when offline

GEO Database: https://www.ncbi.nlm.nih.gov/geo/
"""

import json
import logging
import time
from pathlib import Path

import pandas as pd
import requests

logger = logging.getLogger(__name__)


class GEOConnector:
    """
    Connector for the NCBI GEO (Gene Expression Omnibus) database.

    Download flow:
        1. Query the NCBI eUtils API with the dataset accession ID (GSExxx)
        2. Retrieve dataset metadata (title, number of samples, etc.)
        3. Download the series matrix file from the GEO FTP server
        4. Parse the SOFT-format matrix into a count matrix DataFrame
        5. Parse sample metadata and auto-detect condition labels

    Attributes:
        cache_dir (Path): Local directory where downloaded data is cached.
        base_url  (str):  Base URL for the NCBI eUtils REST API.
    """

    BASE_URL   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    GEO_URL    = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    MIRROR_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/"

    def __init__(self, cache_dir: str = "data/cache/", config: dict = None):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.config = config or {}
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "RNAseq-Pipeline/1.0"})

    # ── Public API ──────────────────────────────────────────────────────────

    def fetch(self, geo_id: str, use_cache: bool = True):
        """
        Download or load from cache a GEO dataset.

        Args:
            geo_id    (str):  GEO accession ID, e.g. "GSE47756".
            use_cache (bool): If True, load from local cache when available.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]:
                counts_df   — Count matrix (genes x samples).
                metadata_df — Sample metadata with 'condition' column.
        """
        logger.info(f"  Fetching GEO dataset: {geo_id}")
        cache_counts   = self.cache_dir / f"{geo_id}_counts.csv"
        cache_metadata = self.cache_dir / f"{geo_id}_metadata.csv"

        if use_cache and cache_counts.exists() and cache_metadata.exists():
            logger.info("  Loading from local cache...")
            counts_df   = pd.read_csv(cache_counts,   index_col=0)
            metadata_df = pd.read_csv(cache_metadata, index_col=0)
            self._log_dataset_info(counts_df, metadata_df)
            return counts_df, metadata_df

        # Try to download from GEO; fall back to synthetic demo data if offline
        try:
            counts_df, metadata_df = self._download_geo(geo_id)
        except Exception as exc:
            logger.warning(f"  Could not download {geo_id}: {exc}")
            logger.warning("  Generating synthetic demo data for demonstration purposes.")
            counts_df, metadata_df = self._generate_demo_data(geo_id)

        # Cache to disk for future runs
        counts_df.to_csv(cache_counts)
        metadata_df.to_csv(cache_metadata)
        logger.info(f"  Data cached to: {self.cache_dir}")
        self._log_dataset_info(counts_df, metadata_df)
        return counts_df, metadata_df

    # ── Download helpers ────────────────────────────────────────────────────

    def _download_geo(self, geo_id: str):
        """Download real expression data from NCBI GEO via eUtils API."""
        logger.info(f"  Querying NCBI GEO for {geo_id}...")

        # Step 1: Search for the accession ID
        search_url = (
            f"{self.BASE_URL}esearch.fcgi"
            f"?db=gds&term={geo_id}[Accession]&retmode=json"
        )
        resp = self.session.get(search_url, timeout=30)
        resp.raise_for_status()
        ids = resp.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            raise ValueError(f"GEO accession {geo_id} not found in NCBI database.")

        # Step 2: Fetch summary metadata
        summary_url = (
            f"{self.BASE_URL}esummary.fcgi"
            f"?db=gds&id={ids[0]}&retmode=json"
        )
        resp = self.session.get(summary_url, timeout=30)
        resp.raise_for_status()
        result = resp.json().get("result", {}).get(ids[0], {})
        logger.info(f"  Title   : {result.get('title', 'N/A')[:80]}")
        logger.info(f"  Samples : {result.get('n_samples', 'N/A')}")

        # Step 3: Build FTP path and download series matrix
        prefix     = geo_id[:5] + "nnn"
        matrix_url = (
            f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}"
            f"/{geo_id}/matrix/{geo_id}_series_matrix.txt.gz"
        )
        logger.info("  Downloading expression matrix...")
        time.sleep(0.5)  # Respect NCBI rate limits

        matrix_resp = self.session.get(matrix_url, timeout=120, stream=True)
        if matrix_resp.status_code != 200:
            raise ConnectionError(f"Matrix file not accessible at: {matrix_url}")

        matrix_file = self.cache_dir / f"{geo_id}_matrix.txt.gz"
        with open(matrix_file, "wb") as f:
            for chunk in matrix_resp.iter_content(chunk_size=8192):
                f.write(chunk)

        return self._parse_geo_matrix(matrix_file, geo_id)

    def _parse_geo_matrix(self, matrix_file: Path, geo_id: str):
        """Parse a GEO SOFT series matrix file into DataFrames."""
        import gzip
        import io

        metadata_rows = []
        data_lines    = []
        in_data       = False

        with gzip.open(matrix_file, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.rstrip()
                if line.startswith("!Sample_"):
                    metadata_rows.append(line)
                elif line.startswith("!series_matrix_table_begin"):
                    in_data = True
                elif line.startswith("!series_matrix_table_end"):
                    in_data = False
                elif in_data:
                    data_lines.append(line)

        # Parse expression matrix
        counts_df = pd.read_csv(
            io.StringIO("\n".join(data_lines)), sep="\t", index_col=0
        )
        counts_df = counts_df.dropna(how="all").fillna(0)
        counts_df = counts_df.apply(pd.to_numeric, errors="coerce").fillna(0)

        # Parse sample metadata
        meta_dict = {}
        for row in metadata_rows:
            parts = row.split("\t")
            key   = parts[0].replace("!", "").strip()
            vals  = [v.strip('"') for v in parts[1:]]
            meta_dict[key] = vals

        metadata_df = pd.DataFrame(meta_dict)
        if "Sample_title" in metadata_df.columns:
            metadata_df.index = metadata_df["Sample_title"]
        metadata_df = self._annotate_conditions(metadata_df)

        # Align samples between matrix and metadata
        common = [c for c in counts_df.columns if c in metadata_df.index]
        if common:
            counts_df   = counts_df[common]
            metadata_df = metadata_df.loc[common]

        return counts_df, metadata_df

    def _annotate_conditions(self, metadata_df: pd.DataFrame) -> pd.DataFrame:
        """
        Auto-detect Normal vs Tumor condition labels from metadata columns.
        Searches for keywords like 'tumor', 'cancer', 'carcinoma' in tissue/source fields.
        """
        metadata_df = metadata_df.copy()
        condition_col = None

        # Look for condition-related columns
        for col in metadata_df.columns:
            if any(k in col.lower() for k in ["tissue", "condition", "type", "source", "status"]):
                condition_col = col
                break

        tumor_keywords = ["tumor", "cancer", "carcinoma", "malignant", "neoplasm"]
        if condition_col:
            metadata_df["condition"] = metadata_df[condition_col].apply(
                lambda x: "tumor" if any(t in str(x).lower() for t in tumor_keywords)
                else "normal"
            )
        else:
            # Default: first half = normal, second half = tumor
            n = len(metadata_df)
            metadata_df["condition"] = ["normal"] * (n // 2) + ["tumor"] * (n - n // 2)

        return metadata_df

    # ── Synthetic demo data ─────────────────────────────────────────────────

    def _generate_demo_data(self, geo_id: str):
        """
        Generate realistic synthetic RNAseq data for demonstration.

        Simulates a breast cancer experiment with:
          - 6 normal samples  /  6 tumor samples
          - 5,000 genes including known oncogenes and tumor suppressors
          - Realistic differential expression patterns for cancer genes
        """
        import numpy as np
        np.random.seed(42)

        n_genes   = 5000
        n_normal  = 6
        n_tumor   = 6
        n_samples = n_normal + n_tumor

        # Mix of real cancer-relevant gene symbols + generic gene IDs
        known_genes = [
            "TP53", "BRCA1", "BRCA2", "MYC", "EGFR", "KRAS", "PIK3CA",
            "PTEN", "VHL", "APC", "RB1", "CDH1", "CDKN2A", "MLH1",
            "SMAD4", "ERBB2", "FGFR1", "VEGFA", "CDK4", "MDM2",
            "BCL2", "BAX", "CASP3", "TNF", "IL6", "TGFB1", "STAT3",
            "NF1", "WT1", "MEN1", "RET", "KIT", "PDGFRA", "ALK",
            "ROS1", "MET", "NRAS", "BRAF", "HRAS", "MAP2K1",
        ]
        filler   = [f"GENE{i:04d}" for i in range(n_genes - len(known_genes))]
        all_genes = known_genes + filler

        # Baseline expression using negative binomial distribution
        base_counts = np.random.negative_binomial(
            10, 0.3, size=(n_genes, n_samples)
        ).astype(float)

        # Introduce biologically realistic differential expression
        upregulated_in_tumor   = ["MYC", "EGFR", "VEGFA", "BCL2", "STAT3", "ERBB2", "CDK4"]
        downregulated_in_tumor = ["TP53", "BRCA1", "PTEN", "RB1", "CDKN2A", "BAX", "CASP3"]

        for gene in upregulated_in_tumor:
            idx = all_genes.index(gene)
            base_counts[idx, n_normal:] *= np.random.uniform(3, 8)

        for gene in downregulated_in_tumor:
            idx = all_genes.index(gene)
            base_counts[idx, n_normal:] *= np.random.uniform(0.1, 0.4)

        # Add biological noise
        base_counts += np.random.normal(0, 5, base_counts.shape).clip(0)
        base_counts  = base_counts.round().clip(0)

        sample_names = (
            [f"Normal_{i+1:02d}" for i in range(n_normal)] +
            [f"Tumor_{i+1:02d}"  for i in range(n_tumor)]
        )

        counts_df = pd.DataFrame(
            base_counts, index=all_genes, columns=sample_names
        )

        metadata_df = pd.DataFrame({
            "sample_id": sample_names,
            "condition": ["normal"] * n_normal + ["tumor"] * n_tumor,
            "tissue":    ["breast"] * n_samples,
            "geo_id":    [geo_id]   * n_samples,
            "batch":     (["batch1"] * 3 + ["batch2"] * 3) * 2,
        }, index=sample_names)

        logger.info(
            f"  Synthetic dataset: {n_genes:,} genes x {n_samples} samples "
            f"({n_normal} normal, {n_tumor} tumor)"
        )
        return counts_df, metadata_df

    # ── Info helpers ────────────────────────────────────────────────────────

    def _log_dataset_info(self, counts_df: pd.DataFrame, metadata_df: pd.DataFrame):
        logger.info(f"  Genes   : {counts_df.shape[0]:,}")
        logger.info(f"  Samples : {counts_df.shape[1]}")
        if "condition" in metadata_df.columns:
            for cond, n in metadata_df["condition"].value_counts().items():
                logger.info(f"    - {cond}: {n} samples")

    def list_recommended_datasets(self) -> list:
        """Return a curated list of GEO datasets recommended for beginners."""
        return [
            {"id": "GSE47756", "tissue": "Breast",      "samples": 12,  "note": "Good for beginners"},
            {"id": "GSE19804", "tissue": "Lung",         "samples": 60,  "note": "Lung adenocarcinoma"},
            {"id": "GSE14520", "tissue": "Liver",        "samples": 247, "note": "Hepatocellular carcinoma"},
            {"id": "GSE22093", "tissue": "Colorectal",   "samples": 290, "note": "Colorectal cancer"},
        ]
