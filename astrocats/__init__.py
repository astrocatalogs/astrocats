"""Astrocats.

Scripts for creating and analyzing catalogs of astronomical data.
"""

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import os
import sys

__version__ = '0.3.42'
__author__ = 'James Guillochon & Luke Kelley'
__license__ = 'MIT'

_PATH_ROOT = os.path.join(os.path.dirname(__file__), "")
_PATH_CATALOG = os.path.join(_PATH_ROOT, "catalog", "")
_PATH_SCHEMA = os.path.join(_PATH_ROOT, "schema", "")

warnings.warn("Adding `pyastroschema` to `sys.path` for easy access... fix this")
_PAS_PATH = "/Users/lzkelley/Research/catalogs/astroschema"
if _PAS_PATH not in sys.path:
    sys.path.append(_PAS_PATH)

from .utils import gitter

__git_version__ = gitter.get_git()
