"""Astrocats: Production submodule.
"""
# flake8: noqa  --- ignore imported but unused flake8 warnings
import os

_DIR_THIS = os.path.realpath(os.path.dirname(__file__))
DIR_TEMPLATES = os.path.join(_DIR_THIS, "templates", "")

from . import director
from . import producer
from . import utils
