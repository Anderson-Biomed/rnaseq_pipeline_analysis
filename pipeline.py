"""
pipeline.py
===========
RNAseq Analysis Pipeline — Main Orchestrator
Normal Tissue vs Cancer Tumor: Gene Expression & Silencing

Pipeline steps:
    1. Data acquisition (GEO database or local files)
    2. Quality control (QC)
    3. Normalization
    4. Differential expression analysis
    5. Visualizations
    6. Final report
"""

import argparse
import logging
import sys
from pathlib import Path

from modules.db_connector import GEOConnector
from modules.qc import QualityControl
from modules.normalization import Normalizer
from modules.differential_expression import DifferentialExpression
from modules.visualizer import Visualizer
from modules.reporter import Reporter
from utils.logger import setup_logger
from utils.config_loader import load_config


def parse_args():
    parser = argparse.ArgumentParser(
        description="RNAseq Pipeline: Normal vs Tumor Expression Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run full pipeline downloading from GEO database
  python pipeline.py --geo-id GSE47756 --output results/

  # Run with your own local count matrix + metadata
  python pipeline.py --counts data/my_counts.csv --metadata data/my_metadata.csv --output results/

  # Skip download and use previously cached data
  python pipeline.py --geo-id GSE47756 --use-cache --output results/

  # Use a custom config file
  python pipeline.py --geo-id GSE47756 --config config/my_config.json --output results/
        """
    )
    parser.add_argument(
        "--geo-id",
        type=str,
        help="GEO dataset accession ID (e.g. GSE47756). Downloads from NCBI GEO."
    )
    parser.add_argument(
        "--counts",
        type=str,
        help="Path to local count matrix CSV (genes as rows, samples as columns)."
    )
    parser.add_argument(
        "--metadata",
        type=str,
        help="Path to sample metadata CSV. Must contain a 'condition' column (normal/tumor)."
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config/pipeline_config.json",
        help="Path to JSON configuration file (default: config/pipeline_config.json)."
    )
    parser.add_argument(
        "--output",
        type=str,
        default="results/",
        help="Output directory for all results and plots (default: results/)."
    )
    parser.add_argument(
        "--use-cache",
        action="store_true",
        help="Load previously downloaded data from cache instead of re-downloading."
    )
    parser.add_argument(
        "--skip-qc",
        action="store_true",
        help="Skip quality control step."
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity level (default: INFO)."
    )
    return parser.parse_args()


def run_pipeline(args):
    # ── Setup ──────────────────────────────────────────────────────────────
    logger = setup_logger(args.log_level)
    config = load_config(args.config)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 62)
    logger.info("   RNAseq Analysis Pipeline  —  Normal vs Tumor")
    logger.info("=" * 62)

    # ── Step 1: Data Acquisition ───────────────────────────────────────────
    logger.info("\n[STEP 1/6]  Data Acquisition")

    if args.geo_id:
        connector = GEOConnector(cache_dir="data/cache/", config=config)
        counts_df, metadata_df = connector.fetch(
            geo_id=args.geo_id,
            use_cache=args.use_cache
        )
    elif args.counts and args.metadata:
        import pandas as pd
        counts_df   = pd.read_csv(args.counts,   index_col=0)
        metadata_df = pd.read_csv(args.metadata, index_col=0)
        logger.info(f"  Loaded local data: "
                    f"{counts_df.shape[0]:,} genes x {counts_df.shape[1]} samples")
    else:
        logger.error("  Provide --geo-id OR both --counts + --metadata")
        sys.exit(1)

    # ── Step 2: Quality Control ────────────────────────────────────────────
    if not args.skip_qc:
        logger.info("\n[STEP 2/6]  Quality Control")
        qc = QualityControl(output_dir=output_dir / "qc", config=config)
        counts_df, metadata_df = qc.run(counts_df, metadata_df)
    else:
        logger.info("\n[STEP 2/6]  Quality Control  →  SKIPPED")

    # ── Step 3: Normalization ──────────────────────────────────────────────
    logger.info("\n[STEP 3/6]  Normalization")
    normalizer = Normalizer(method=config["normalization"]["method"])
    norm_df = normalizer.normalize(counts_df)

    # ── Step 4: Differential Expression ───────────────────────────────────
    logger.info("\n[STEP 4/6]  Differential Expression Analysis")
    de = DifferentialExpression(
        output_dir=output_dir / "de_results",
        config=config
    )
    results_df = de.run(norm_df, metadata_df)

    # ── Step 5: Visualizations ─────────────────────────────────────────────
    logger.info("\n[STEP 5/6]  Generating Visualizations")
    viz = Visualizer(output_dir=output_dir / "plots")
    viz.plot_pca(norm_df, metadata_df)
    viz.plot_volcano(results_df, config)
    viz.plot_heatmap(norm_df, results_df, metadata_df, config)
    viz.plot_ma(results_df)
    viz.plot_gene_expression_boxplots(norm_df, metadata_df, results_df, config)

    # ── Step 6: Report ─────────────────────────────────────────────────────
    logger.info("\n[STEP 6/6]  Generating Report")
    reporter = Reporter(output_dir=output_dir)
    reporter.generate(
        counts_df=counts_df,
        norm_df=norm_df,
        metadata_df=metadata_df,
        de_results=results_df,
        config=config,
        geo_id=args.geo_id
    )

    logger.info("\n" + "=" * 62)
    logger.info(f"  Pipeline complete!  Results saved to: {output_dir}/")
    logger.info("=" * 62)


if __name__ == "__main__":
    args = parse_args()
    run_pipeline(args)
