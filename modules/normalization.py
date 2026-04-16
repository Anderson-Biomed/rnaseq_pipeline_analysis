"""
modules/normalization.py
========================
Normalization of RNAseq count matrices.

Why normalize?
    Different samples are sequenced to different depths (library sizes).
    A gene appearing 1,000 times in sample A and 500 times in sample B
    is NOT necessarily more expressed in A — sample A may simply have
    twice as many total reads. Normalization corrects for this.

Methods implemented:
    CPM   — Counts Per Million. Simplest method. Corrects for library size only.
    TPM   — Transcripts Per Million. Corrects for library size AND gene length.
    TMM   — Trimmed Mean of M-values. edgeR-style robust normalization.
    DESeq2— Median-of-ratios. Gold standard for bulk RNA-seq experiments.

Recommendation for beginners: start with CPM or DESeq2.
"""

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class Normalizer:
    """
    Normalize a raw count matrix using one of four methods.

    After normalization the values are log2-transformed:
        log2(normalized_value + 1)
    This stabilizes variance and makes downstream statistics more robust.

    Args:
        method (str): One of 'cpm', 'tpm', 'tmm', 'deseq2'.
    """

    SUPPORTED_METHODS = ["cpm", "tpm", "tmm", "deseq2"]

    def __init__(self, method: str = "cpm"):
        method = method.lower()
        if method not in self.SUPPORTED_METHODS:
            raise ValueError(
                f"Unknown normalization method '{method}'. "
                f"Choose from: {self.SUPPORTED_METHODS}"
            )
        self.method = method

    def normalize(
        self,
        counts_df: pd.DataFrame,
        gene_lengths: pd.Series = None
    ) -> pd.DataFrame:
        """
        Normalize the count matrix.

        Args:
            counts_df    (pd.DataFrame): Raw count matrix (genes x samples).
            gene_lengths (pd.Series):   Gene lengths in base pairs, required for TPM.
                                        Index must match counts_df.index.

        Returns:
            pd.DataFrame: log2(normalized + 1) matrix, same shape as input.
        """
        logger.info(f"  Method : {self.method.upper()}")

        if self.method == "cpm":
            norm_df = self._cpm(counts_df)
        elif self.method == "tpm":
            norm_df = self._tpm(counts_df, gene_lengths)
        elif self.method == "tmm":
            norm_df = self._tmm(counts_df)
        elif self.method == "deseq2":
            norm_df = self._deseq2_median_of_ratios(counts_df)

        # Log2 transform for variance stabilization
        log_norm_df = np.log2(norm_df + 1)

        logger.info(
            f"  Value range after log2 transform: "
            f"[{log_norm_df.min().min():.2f}, {log_norm_df.max().max():.2f}]"
        )
        return log_norm_df

    # ── CPM ─────────────────────────────────────────────────────────────────

    def _cpm(self, counts_df: pd.DataFrame) -> pd.DataFrame:
        """
        Counts Per Million (CPM).

        Formula:
            CPM_i = (count_i / total_library_counts) x 1,000,000

        Corrects for: library size.
        Does NOT correct for: gene length.
        Best for: simple comparisons, short genes, or when gene length is unknown.
        """
        lib_sizes = counts_df.sum(axis=0)  # total counts per sample
        cpm       = counts_df.divide(lib_sizes, axis=1) * 1_000_000
        logger.info(f"  CPM range: [{cpm.min().min():.2f}, {cpm.max().max():.2f}]")
        return cpm

    # ── TPM ─────────────────────────────────────────────────────────────────

    def _tpm(self, counts_df: pd.DataFrame, gene_lengths: pd.Series = None) -> pd.DataFrame:
        """
        Transcripts Per Million (TPM).

        Formula:
            Step 1 — Reads Per Kilobase (RPK):
                RPK_i = count_i / (gene_length_i / 1000)
            Step 2 — Normalize to 1 million:
                TPM_i = (RPK_i / sum_of_all_RPK) x 1,000,000

        Corrects for: library size AND gene length.
        Best for: comparing between samples when gene length varies significantly.
        Note: TPM values sum to exactly 1,000,000 per sample.
        """
        if gene_lengths is None:
            logger.warning(
                "  TPM: No gene lengths provided — "
                "using log-normal simulated lengths (for demo only)."
            )
            np.random.seed(99)
            gene_lengths = pd.Series(
                np.random.lognormal(mean=7, sigma=1, size=len(counts_df)),
                index=counts_df.index
            )

        gene_lengths = gene_lengths.reindex(counts_df.index).fillna(1000)
        rpk          = counts_df.divide(gene_lengths / 1000, axis=0)
        scale        = rpk.sum(axis=0) / 1_000_000
        tpm          = rpk.divide(scale, axis=1)
        logger.info(f"  TPM range: [{tpm.min().min():.2f}, {tpm.max().max():.2f}]")
        return tpm

    # ── TMM ─────────────────────────────────────────────────────────────────

    def _tmm(
        self,
        counts_df: pd.DataFrame,
        trim_m: float = 0.30,
        trim_a: float = 0.05
    ) -> pd.DataFrame:
        """
        Trimmed Mean of M-values (TMM) — edgeR-style normalization.

        Algorithm:
            1. Select a reference sample (closest to 75th percentile library size).
            2. For each other sample compute M values (log fold changes) and
               A values (average log expression) relative to the reference.
            3. Remove extreme M and A values (trimming).
            4. Compute weighted mean of remaining M values as the normalization factor.
            5. Scale all factors so their geometric mean equals 1.

        Corrects for: composition bias (highly expressed genes dominating counts).
        Best for: experiments where a few genes are very highly expressed.
        """
        lib_sizes = counts_df.sum(axis=0)
        # Reference sample: closest to the 75th percentile library size
        ref_col   = (lib_sizes - lib_sizes.quantile(0.75)).abs().idxmin()
        ref       = counts_df[ref_col]

        norm_factors = {}
        for col in counts_df.columns:
            sample = counts_df[col]
            mask   = (ref > 0) & (sample > 0)
            ref_s, sample_s = ref[mask], sample[mask]

            # M = log2 fold change (sample vs reference, library-size adjusted)
            M = (
                np.log2(sample_s / sample_s.sum()) -
                np.log2(ref_s    / ref_s.sum())
            )
            # A = average log expression
            A = 0.5 * (
                np.log2(sample_s / sample_s.sum()) +
                np.log2(ref_s    / ref_s.sum())
            )

            # Trim extreme values
            m_lo, m_hi = M.quantile(trim_m / 2), M.quantile(1 - trim_m / 2)
            a_lo, a_hi = A.quantile(trim_a / 2), A.quantile(1 - trim_a / 2)
            keep = (M >= m_lo) & (M <= m_hi) & (A >= a_lo) & (A <= a_hi)

            if keep.sum() > 0:
                w = 1.0 / (1.0 / sample_s[keep] + 1.0 / ref_s[keep])
                norm_factors[col] = 2 ** float(np.average(M[keep], weights=w))
            else:
                norm_factors[col] = 1.0

        # Scale to geometric mean of 1
        nf          = pd.Series(norm_factors)
        nf          = nf / np.exp(np.log(nf).mean())
        effective   = lib_sizes * nf
        cpm         = counts_df.divide(effective, axis=1) * 1_000_000

        logger.info(f"  TMM factors: min={nf.min():.3f}, max={nf.max():.3f}")
        return cpm

    # ── DESeq2 Median-of-Ratios ──────────────────────────────────────────────

    def _deseq2_median_of_ratios(self, counts_df: pd.DataFrame) -> pd.DataFrame:
        """
        DESeq2-style Median-of-Ratios normalization.

        Algorithm:
            1. For each gene, compute its geometric mean across all samples.
            2. Divide each sample's counts by the gene-wise geometric means
               to get per-sample ratios.
            3. The size factor for a sample is the median of its ratios.
            4. Divide raw counts by the size factor.

        This method is robust to highly differentially expressed genes because
        the geometric mean provides a stable reference that outlier genes
        cannot distort.

        Best for: standard bulk RNA-seq; gold standard alongside TMM.
        """
        # Compute log geometric means; exclude genes with any zero count
        log_counts     = np.log(counts_df.replace(0, np.nan))
        log_geo_means  = log_counts.mean(axis=1)   # mean of logs = log of geometric mean
        valid          = np.isfinite(log_geo_means)

        log_counts_v   = log_counts[valid]
        log_geo_v      = log_geo_means[valid]

        # Per-sample ratio relative to the geometric mean reference
        log_ratios     = log_counts_v.subtract(log_geo_v, axis=0)
        size_factors   = np.exp(log_ratios.median(axis=0))

        normalized     = counts_df.divide(size_factors, axis=1)
        logger.info(
            f"  DESeq2 size factors: "
            f"min={size_factors.min():.3f}, max={size_factors.max():.3f}"
        )
        return normalized
