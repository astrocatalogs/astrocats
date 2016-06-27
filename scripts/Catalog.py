import os
import sys
from collections import OrderedDict

from git import Repo

from scripts import PATH

from .utils import pbar, tprint
from .importer.funcs import (get_bibauthor_dict, get_biberror_dict,
                             get_extinctions_dict, get_repos_dict)


class Catalog():
    """
    Object to hold the main catalog dictionary and other catalog globals.
    """

    nedd_dict = OrderedDict()
    repos_dict = get_repos_dict()
    bibauthor_dict = get_bibauthor_dict()
    biberror_dict = get_biberror_dict()
    extinctions_dict = get_extinctions_dict()

    def __init__(self, log):
        all_repos = [self.repos_dict[x] for x in self.repos_dict]
        all_repos = [i for x in all_repos for i in x]
        for repo in pbar(all_repos):
            if not os.path.isdir(PATH.ROOT + "/" + repo):
                try:
                    log.warning('Cloning "' + repo +
                                '" (only needs to be done ' +
                                'once, may take few minutes per repo).')
                    Repo.clone_from("git@github.com:astrocatalogs/" +
                                    repo + ".git",
                                    PATH.ROOT + "/" + repo)
                except:
                    raise
                    sys.exit()
        return
