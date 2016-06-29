"""Manage the import of supernovae data.

Modified from `scripts/import.py`.
"""
import os

from . import supernova

_PATH_SUPERNOVAE = os.path.abspath(os.path.split(os.path.dirname(__file__)))
print("_PATH_SUPERNOVAE = '{}'".format(_PATH_SUPERNOVAE))
