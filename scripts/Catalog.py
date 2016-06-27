import os
import sys
from collections import OrderedDict

from git import Repo

from scripts import PATH

from .utils import logger, pbar
from .importer.funcs import (get_bibauthor_dict, get_biberror_dict,
                             get_extinctions_dict, get_repos_dict)


class Catalog():
    """
    Object to hold the main catalog dictionary and other catalog globals.
    """

    def __init__(self, args):
        # Store runtime arguments
        self.args = args
        # Load a logger object
        # Determine verbosity ('None' means use default)
        log_stream_level = None
        if args.debug:
            log_stream_level = logger.DEBUG
        elif args.verbose:
            log_stream_level = logger.INFO

        # Destination of log-file ('None' means no file)
        self.log = logger.get_logger(
            stream_level=log_stream_level, tofile=args.log_filename)

        # Make sure repositories are cloned
        self._clone_repos()

        # Load all dictionaries
        self._load_dicts()

        return

    def _clone_repos(self):
        all_repos = [self.repos_dict[x] for x in self.repos_dict]
        all_repos = [i for x in all_repos for i in x]
        for repo in pbar(all_repos):
            if not os.path.isdir(PATH.ROOT + "/" + repo):
                try:
                    self.log.warning('Cloning "' + repo +
                                     '" (only needs to be done ' +
                                     'once, may take few minutes per repo).')
                    Repo.clone_from("git@github.com:astrocatalogs/" +
                                    repo + ".git",
                                    PATH.ROOT + "/" + repo)
                except:
                    self.log.error("CLONING '{}' INTERRUPTED".format(repo))
                    raise
                    sys.exit()
        return

    def _load_dicts(self):
        # Create empty `events` collection
        self.events = OrderedDict()
        # Create/Load auxiliary dictionaries
        self.nedd_dict = OrderedDict()
        self.repos_dict = get_repos_dict()
        self.bibauthor_dict = get_bibauthor_dict()
        self.biberror_dict = get_biberror_dict()
        self.extinctions_dict = get_extinctions_dict()
        return
