"""
utils/logger.py
===============
Logging setup for the RNAseq pipeline.

Writes log messages to both stdout and a log file (results/logs/pipeline.log).
"""

import logging
import sys
from pathlib import Path


def setup_logger(level: str = "INFO") -> logging.Logger:
    """
    Configure the root logger with console + file handlers.

    Args:
        level (str): Logging level — DEBUG, INFO, WARNING, or ERROR.

    Returns:
        logging.Logger: Configured root logger.
    """
    log_dir = Path("results/logs")
    log_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger()
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))
    logger.handlers.clear()

    formatter = logging.Formatter(
        "%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S"
    )

    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # File handler
    fh = logging.FileHandler(log_dir / "pipeline.log", mode="a")
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    return logger
