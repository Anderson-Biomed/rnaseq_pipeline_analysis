"""
utils/config_loader.py
======================
Load pipeline configuration from a JSON file.

If the file does not exist, built-in defaults are used.
User-provided values override defaults; missing keys fall back to defaults.
"""

import json
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# Default configuration — all values can be overridden in pipeline_config.json
DEFAULTS = {
    "qc": {
        "min_count":   10,   # Minimum read count to keep a gene
        "min_samples":  3,   # Minimum number of samples that must meet min_count
        "outlier_sd":   3    # Standard deviations threshold for outlier detection
    },
    "normalization": {
        "method": "cpm"      # Options: cpm | tpm | tmm | deseq2
    },
    "differential_expression": {
        "test_method":    "welch",  # Options: welch | mannwhitney
        "padj_threshold":  0.05,   # FDR significance threshold
        "lfc_threshold":   1.0     # Minimum |log2FC| to call a gene significant
    },
    "visualizations": {
        "top_heatmap_genes": 50,
        "top_boxplot_genes": 12,
        "dpi": 150
    }
}


def load_config(path: str) -> dict:
    """
    Load configuration from a JSON file, merged with defaults.

    Args:
        path (str): Path to the JSON configuration file.

    Returns:
        dict: Merged configuration dictionary.
    """
    config_path = Path(path)
    if config_path.exists():
        with open(config_path) as f:
            user_config = json.load(f)
        merged = _deep_merge(DEFAULTS, user_config)
        logger.info(f"  Config loaded from: {config_path}")
        return merged
    else:
        logger.warning(
            f"  Config file not found at '{config_path}' — using built-in defaults."
        )
        return DEFAULTS.copy()


def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge two dicts; override values take precedence."""
    result = base.copy()
    for key, val in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(val, dict):
            result[key] = _deep_merge(result[key], val)
        else:
            result[key] = val
    return result
