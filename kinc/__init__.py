"""
KINC - Knowledge Independent Network Construction

A Python implementation for constructing gene co-expression networks
from gene expression matrices using correlation and GMM-based clustering.
"""

__version__ = "4.0.0a"
__author__ = "KINC Development Team"

from . import data
from . import analytics
from . import core
from . import parallel

__all__ = ["data", "analytics", "core", "parallel", "__version__"]
