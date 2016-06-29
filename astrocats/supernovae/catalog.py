"""Supernovae specific catalog class.
"""
from collections import OrderedDict
import os

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
        self.FILENAME = FILENAME

        return

    def _load_aux_data(self):
        # Create/Load auxiliary dictionaries
        self.nedd_dict = OrderedDict()
        self.repos_dict = read_json_dict(FILENAME.REPOS)
        self.bibauthor_dict = read_json_dict(FILENAME.BIBAUTHORS)
        self.biberror_dict = read_json_dict(FILENAME.BIBERRORS)
        self.extinctions_dict = read_json_dict(FILENAME.EXTINCT)
        # Create/Load auxiliary arrays
        self.nonsneprefixes_dict = read_json_arr(FILENAME.NON_SNE_PREFIXES)
        self.nonsnetypes = read_json_arr(FILENAME.NON_SNE_TYPES)
        return

    def clone_repos(self):
        # Load the local 'supernovae' repository names
        # all_repos = [self.repos_dict[x] for x in self.repos_dict]
        # all_repos = [i for x in all_repos for i in x]
        all_repos = self._get_input_repo_folders()
        all_repos += self.get_output_repo_folders()
        super()._clone_repos(all_repos)
        return

    def get_repo_file_list(self, normal=True, bones=True):
        repo_folders = self.get_output_repo_folders()
        return super()._get_repo_file_list(
            repo_folders, normal=normal, bones=bones)

    def _get_input_repo_folders(self):
        """
        """
        repo_folders = []
        repo_folders += self.repos_dict['external']
        repo_folders += self.repos_dict['internal']
        repo_folders = [os.path.join(FILENAME.PATH_INPUT, rf)
                        for rf in repo_folders]
        return repo_folders

    def get_output_repo_folders(self):
        """
        """
        repo_folders = []
        repo_folders += self.repos_dict['output']
        repo_folders += self.repos_dict['boneyard']
        repo_folders = [os.path.join(FILENAME.PATH_OUTPUT, rf)
                        for rf in repo_folders]
        return repo_folders

    def get_repo_years(self):
        """
        """
        repo_folders = self.get_output_repo_folders()
        repo_years = [int(repo_folders[x][-4:])
                      for x in range(len(repo_folders) - 1)]
        repo_years[0] -= 1
        return repo_years
