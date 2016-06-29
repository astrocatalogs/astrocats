"""Supernovae specific catalog class.
"""
from collections import OrderedDict
import astrocats.catalog
from astrocats.catalog.utils.imports import read_json_arr, read_json_dict

from . import _PATH_SUPERNOVAE


class Catalog(astrocats.catalog.catalog.Catalog):

    def __init__(self):
        """
        """
        # Initialize super `astrocats.catalog.catalog.Catalog` object
        super().__init__(self)

        self._load_aux()

    def _load_aux(self):
        # Create/Load auxiliary dictionaries
        self.nedd_dict = OrderedDict()
        self.bibauthor_dict = read_json_dict(FILENAME.BIBAUTHORS)
        self.biberror_dict = read_json_dict(FILENAME.BIBERRORS)
        self.extinctions_dict = read_json_dict(FILENAME.EXTINCT)
        # Create/Load auxiliary arrays
        self.nonsneprefixes_dict = read_json_arr(FILENAME.NON_SNE_PREFIXES)
        self.nonsnetypes = read_json_arr(FILENAME.NON_SNE_TYPES)
        return
