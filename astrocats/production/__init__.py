"""Astrocats: Production submodule.
"""
# flake8: noqa  --- ignore imported but unused flake8 warnings
import os

_DIR_THIS = os.path.realpath(os.path.dirname(__file__))
DIR_TEMPLATES = os.path.join(_DIR_THIS, "templates", "")

PATH_SDSS_MISSING_HOST_IMAGE = os.path.join(DIR_TEMPLATES, "sdss_missing_host.jpg")
PATH_SDSS_FAILED_HOST_IMAGE = os.path.join(DIR_TEMPLATES, "sdss_failed_host.jpg")

from . import director
from . import producer
from . import utils
from . import html_pro
from . import host_image_pro
