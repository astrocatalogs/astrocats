"""Supernovae specific catalog class.
"""
from collections import OrderedDict

import astrocats.catalog
from astrocats.supernovae.supernova import Supernova
from astrocats.catalog.utils.imports import read_json_arr, read_json_dict

from .constants import FILENAME


class Catalog(astrocats.catalog.catalog.Catalog):

    def __init__(self, args):
        """
        """
        self.proto = Supernova
        # Initialize super `astrocats.catalog.catalog.Catalog` object
        super().__init__(args)

        self._load_aux_data()

        return

    def _load_aux_data(self):
        # Create/Load auxiliary dictionaries
        self.nedd_dict = OrderedDict()
        self.bibauthor_dict = read_json_dict(FILENAME.BIBAUTHORS)
        self.biberror_dict = read_json_dict(FILENAME.BIBERRORS)
        self.extinctions_dict = read_json_dict(FILENAME.EXTINCT)
        # Create/Load auxiliary arrays
        self.nonsneprefixes_dict = read_json_arr(FILENAME.NON_SNE_PREFIXES)
        self.nonsnetypes = read_json_arr(FILENAME.NON_SNE_TYPES)
        return

    def _clone_repos(self):
        # Load the local 'supernovae' repository names
        repos_dict = read_json_dict(FILENAME.REPOS)
        all_repos = [repos_dict[x] for x in repos_dict]
        all_repos = [i for x in all_repos for i in x]

        super().clone_repos(all_repos)
        return
