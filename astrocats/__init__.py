"""Astrocats.

Scripts for creating and analyzing catalogs of astronomical data.
"""

import os

__version__ = '0.3.42'
__author__ = 'James Guillochon & Luke Kelley'
__license__ = 'MIT'

# Set the path for the user's configuration file
_CONFIG_PATH = os.path.join(
    os.path.expanduser('~'), '.config', 'astrocats', 'astrocatsrc')

_PATH_ROOT = os.path.join(os.path.dirname(__file__), "")
_PATH_CATALOG = os.path.join(_PATH_ROOT, "catalog", "")
_PATH_SCHEMA = os.path.join(_PATH_ROOT, "schema", "")
