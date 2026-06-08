"""Logging configuration module.

Author: Y Chen, University of Reading, 2025
"""
import logging

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(levelname)s: Desroziers' diagnostics %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)
