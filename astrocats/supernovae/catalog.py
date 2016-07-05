"""Supernovae specific catalog class.
"""
import os
from collections import OrderedDict
from subprocess import check_output

import astrocats.catalog
from astrocats.catalog.utils.imports import read_json_arr, read_json_dict
from astrocats.supernovae.supernova import Supernova


class Catalog(astrocats.catalog.catalog.Catalog):

    class PATHS(astrocats.catalog.catalog.Catalog.PATHS):

        PATH_BASE = os.path.abspath(os.path.dirname(__file__))

        def __init__(self, catalog):
            super().__init__(catalog)
            # auxiliary datafiles
            self.TYPE_SYNONYMS = os.path.join(
                self.PATH_INPUT, 'type-synonyms.json')
            self.SOURCE_SYNONYMS = os.path.join(
                self.PATH_INPUT, 'source-synonyms.json')
            self.URL_REDIRECTS = os.path.join(
                self.PATH_INPUT, 'url-redirects.json')
            self.NON_SNE_TYPES = os.path.join(
                self.PATH_INPUT, 'non-sne-types.json')
            self.NON_SNE_PREFIXES = os.path.join(
                self.PATH_INPUT, 'non-sne-prefixes.json')
            self.BIBERRORS = os.path.join(self.PATH_INPUT, 'biberrors.json')
            self.ATELS = os.path.join(self.PATH_INPUT, 'atels.json')
            self.CBETS = os.path.join(self.PATH_INPUT, 'cbets.json')
            self.IAUCS = os.path.join(self.PATH_INPUT, 'iaucs.json')
            # cached datafiles
            self.BIBAUTHORS = os.path.join(
                self.PATH_OUTPUT, 'cache', 'bibauthors.json')
            self.EXTINCT = os.path.join(
                self.PATH_OUTPUT, 'cache', 'extinctions.json')

        def get_repo_output_file_list(self, normal=True, bones=True):
            repo_folders = self.get_repo_output_folders()
            return super()._get_repo_file_list(
                repo_folders, normal=normal, bones=bones)

        def get_repo_input_folders(self):
            """
            """
            repo_folders = []
            repo_folders += self.repos_dict['external']
            repo_folders += self.repos_dict['internal']
            repo_folders = [os.path.join(self.PATH_INPUT, rf)
                            for rf in repo_folders]
            return repo_folders

        def get_repo_output_folders(self):
            """
            """
            repo_folders = []
            repo_folders += self.repos_dict['output']
            repo_folders += self.repos_dict['boneyard']
            repo_folders = [os.path.join(self.PATH_OUTPUT, rf)
                            for rf in repo_folders]
            return repo_folders

        def get_repo_years(self):
            """
            """
            repo_folders = self.get_repo_output_folders()
            repo_years = [int(repo_folders[x][-4:])
                          for x in range(len(repo_folders) - 1)]
            repo_years[0] -= 1
            return repo_years

    class SCHEMA:
        HASH = (check_output(['git', 'log', '-n', '1', '--format="%H"',
                              '--',
                              'OSC-JSON-format.md'])
                .decode('ascii').strip().strip('"').strip())
        URL = ('https://github.com/astrocatalogs/sne/blob/' + HASH +
               '/OSC-JSON-format.md')

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
        self.bibauthor_dict = read_json_dict(self.PATHS.BIBAUTHORS)
        self.biberror_dict = read_json_dict(self.PATHS.BIBERRORS)
        self.extinctions_dict = read_json_dict(self.PATHS.EXTINCT)
        self.iaucs_dict = read_json_dict(self.PATHS.IAUCS)
        self.cbets_dict = read_json_dict(self.PATHS.CBETS)
        self.atels_dict = read_json_dict(self.PATHS.ATELS)
        self.source_syns = read_json_dict(self.PATHS.SOURCE_SYNONYMS)
        self.url_redirs = read_json_dict(self.PATHS.URL_REDIRECTS)
        self.type_syns = read_json_dict(self.PATHS.TYPE_SYNONYMS)
        # Create/Load auxiliary arrays
        self.nonsneprefixes_dict = read_json_arr(
            self.PATHS.NON_SNE_PREFIXES)
        self.nonsnetypes = read_json_arr(self.PATHS.NON_SNE_TYPES)
        return

    def clone_repos(self):
        # Load the local 'supernovae' repository names
        all_repos = self.PATHS.get_repo_input_folders()
        all_repos += self.PATHS.get_repo_output_folders()
        super()._clone_repos(all_repos)
        return
