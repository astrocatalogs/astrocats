"""
"""

import importlib
import os
import sys
from datetime import datetime

from astrocats import __version__, __git_version__, argshandler
from astrocats.utils import logger

_BASE_PATH_KEY = 'base_path'
_PROFILE = False


def main():
    """Primary entry point for all AstroCats catalogs.

    From this entry point, all internal catalogs can be accessed and their
    public methods executed (for example: import scripts).

    """

    # Initialize: Load command line arguments, log, etc
    # ----------------------------------------------------
    args = argshandler.parse_args()
    if args is None:
        return

    # Create a logging object
    log = load_log(args)

    title_str = "Astrocats, version: {}, SHA: {}".format(__version__, __git_version__)
    log.warning("\n{}\n{}\n".format(title_str, '=' * len(title_str)))

    beg_time = datetime.now()
    log.info("{}\n".format(beg_time.ctime()))

    log.info("command: '{}'".format(args.command))
    log.info("catalog: '{}'".format(args.catalog))

    # Load target catalog module
    # ----------------------------------------

    catalog_path = args.catalog
    catalog_path = os.path.realpath(os.path.abspath(catalog_path))
    log.debug("Importing specified catalog module: '{}'".format(catalog_path))

    if not os.path.exists(catalog_path):
        log.raise_error("`catalog_path` '{}' does not exist!".format(catalog_path))

    if not os.path.isdir(catalog_path):
        log.raise_error("`catalog_path` '{}' is not a directory!".format(catalog_path))

    # Try to import the specified module
    try:
        # catalog_module = importlib.import_module('.' + mod_name, package='astrocats')
        # Get the location of the module to add to `sys.path`
        catalog_path_dir = os.path.dirname(catalog_path)
        catalog_module_name = os.path.basename(catalog_path)

        # Insert the path to this module.  Make it first to avoid builtin package conflicts
        log.debug("adding catalog module path: '{}'".format(catalog_path_dir))
        sys.path.insert(0, catalog_path_dir)

        # Import catalog module
        log.debug("importing catalog module name: '{}'".format(catalog_module_name))
        catalog_module = importlib.import_module(catalog_module_name)

        # Remove added path to avoid unintended imports
        sys.path.pop(0)
        log.info("Loaded: '{}'".format(str(catalog_module)))

    except Exception as err:
        msg = "Import of specified catalog module '{}' failed.".format(catalog_path)
        log.raise_error(msg, err)

    # Run target command in the given catalog module
    # -----------------------------------------------------------

    if _PROFILE:
        msg = "RUNNING IN PROFILE MODE"
        log.warning("")
        log.warning("="*len(msg))
        log.warning(msg)
        log.warning("="*len(msg))
        log.warning("")
        import cProfile
        pr = cProfile.Profile()
        pr.enable()

    run_command(args, log, catalog_module)

    if _PROFILE:
        pr.disable()
        pr.print_stats(sort='time')
        dt_str = beg_time.strftime("%Y%m%d-%H%M%S")
        profile_fname = "{:mod}_{:date}_cprofile-run-stats.dat".format(catalog_module_name, dt_str)
        msg = "Saved profile stats to '{}'".format(profile_fname)
        log.warning("")
        log.warning("="*len(msg))
        log.warning(msg)
        log.warning("="*len(msg))
        log.warning("")
        pr.dump_stats(profile_fname)

    end_time = datetime.now()
    log.warning("\nAll complete at {}".format(end_time))
    log.warning("Duration {}\n".format(end_time - beg_time))
    return


def load_log(args):
    """Load a `logging.Logger` object.

    Arguments
    ---------
    args : `argparse.Namespace` object
        Namespace containing required settings:
        {`args.debug`, `args.verbose`, and `args.log_filename`}.

    Returns
    -------
    log : `logging.Logger` object

    """

    # Determine verbosity ('None' means use default)
    log_stream_level = None
    if args.debug:
        log_stream_level = logger.DEBUG
    elif args.verbose:
        log_stream_level = logger.INFO

    # Create log
    log = logger.get_logger(stream_level=log_stream_level, tofile=args.log_filename)
    log._verbose = args.verbose
    log._debug = args.debug
    return log


def run_command(args, log, catalog_module):
    log.debug(__file__ + ":run_command()")

    try:
        # catalog_module_name = catalog_module.__name__
        # Catalog_Class = importlib.import_module(".Catalog_Class", catalog_module_name)
        print("Catalog_Class = ", catalog_module.Catalog_Class)
        catalog = catalog_module.Catalog_Class(args, log)
        print(catalog)
    except Exception as err:
        raise

    return

    # Load the catalog

    # Data Import
    # -----------
    if args.command == 'import':
        log.log(log_lvl, "Running 'import'.")
        catalog.import_data()

    # Data Export
    # -----------
    elif args.subcommand == 'produce':
        # from . import production
        self.log.log(log_lvl, "Running 'produce'.")
        manager = catalog.Director(catalog, args)
        manager.direct(args)

    # Git Subcommands
    # ---------------
    elif args.subcommand.startswith('git'):
        from . import gitter

        if args.subcommand == 'git-clone':
            self.log.log(log_lvl, "Running 'git clone'.")
            gitter.git_clone_all_repos(catalog)
        elif args.subcommand == 'git-push':
            self.log.log(log_lvl, "Running 'git push'.")
            gitter.git_add_commit_push_all_repos(catalog)
        elif args.subcommand == 'git-pull':
            self.log.log(log_lvl, "Running 'git pull'.")
            gitter.git_pull_all_repos(catalog)
        elif args.subcommand == 'git-reset-local':
            self.log.log(log_lvl, "Running 'git reset' using the local HEAD.")
            gitter.git_reset_all_repos(catalog, hard=True, origin=False, clean=True)
        elif args.subcommand == 'git-reset-origin':
            self.log.log(log_lvl, "Running 'git reset' using 'origin/master'.")
            gitter.git_reset_all_repos(catalog, hard=True, origin=True, clean=True)
        elif args.subcommand == 'git-status':
            self.log.log(log_lvl, "Running 'git status'.")
            gitter.git_status_all_repos(catalog)
        else:
            self.log.error("Unrecognized git subcommand '{}'".format(args.subcommand))

    # Analyze Catalogs
    # ----------------
    elif args.subcommand == 'analyze':
        self.log.log(log_lvl, "Running 'analyze'.")
        from . import analysis
        # Create an `Analysis` instance
        lysis = analysis.analysis.Analysis(catalog, self.log)
        # Pass the command-line arguments to run.
        lysis.analyze(args)



    return


if __name__ == "__main__":
    main()
