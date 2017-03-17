"""
"""
import importlib
import json
import os

from astrocats import _CONFIG_PATH, __version__
from astrocats.catalog.utils import log_raise

_BASE_PATH_KEY = 'base_path'


def main():
    """Primary entry point for all AstroCats catalogs.

    From this entry point, all internal catalogs can be accessed and their
    public methods executed (for example: import scripts).

    """
    from datetime import datetime

    # Initialize Command-Line and User-Config Settings, Log
    # -----------------------------------------------------

    beg_time = datetime.now()
    # Process command-line arguments to determine action
    # If no subcommand (e.g. 'import') is given, returns 'None' --> exit
    args, sub_clargs = load_command_line_args()
    if args is None:
        return

    # Create a logging object
    log = load_log(args)

    # Run configuration/setup interactive script
    if args.command == 'setup':
        setup_user_config(log)
        return

    # Make sure configuration file exists, or that's what we're doing
    # (with the 'setup' subcommand)
    if not os.path.isfile(_CONFIG_PATH):
        raise RuntimeError("'{}' does not exist.  "
                           "Run `astrocats setup` to configure."
                           "".format(_CONFIG_PATH))

    git_vers = get_git()
    title_str = "Astrocats, version: {}, SHA: {}".format(__version__, git_vers)
    log.warning("\n\n{}\n{}\n{}\n".format(title_str, '=' * len(title_str),
                                          beg_time.ctime()))

    # Load the user settings from the home directory
    args = load_user_config(args, log)

    # Choose Catalog and Operation(s) to perform
    # ------------------------------------------
    mod_name = args.command
    log.debug("Importing specified module: '{}'".format(mod_name))
    # Try to import the specified module
    try:
        mod = importlib.import_module('.' + mod_name, package='astrocats')
    except Exception as err:
        log.error("Import of specified module '{}' failed.".format(mod_name))
        log_raise(log, str(err), type(err))

    # Run the `main.main` method of the specified module
    log.debug("Running `main.main()`")
    mod.main.main(args, sub_clargs, log)

    end_time = datetime.now()
    log.warning("\nAll complete at {}, After {}".format(end_time, end_time -
                                                        beg_time))
    return


def setup_user_config(log):
    """Setup a configuration file in the user's home directory.

    Currently this method stores default values to a fixed configuration
    filename.  It should be modified to run an interactive prompt session
    asking for parameters (or at least confirming the default ones).

    Arguments
    ---------
    log : `logging.Logger` object

    """
    log.warning("AstroCats Setup")
    log.warning("Configure filepath: '{}'".format(_CONFIG_PATH))

    # Create path to configuration file as needed
    config_path_dir = os.path.split(_CONFIG_PATH)[0]
    if not os.path.exists(config_path_dir):
        log.debug("Creating config directory '{}'".format(config_path_dir))
        os.makedirs(config_path_dir)

    if not os.path.isdir(config_path_dir):
        log_raise(log, "Configure path error '{}'".format(config_path_dir))

    # Determine default settings

    # Get this containing directory and use that as default data path
    def_base_path = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
    log.warning("Setting '{}' to default path: '{}'".format(_BASE_PATH_KEY,
                                                            def_base_path))
    config = {_BASE_PATH_KEY: def_base_path}

    # Write settings to configuration file
    json.dump(config, open(_CONFIG_PATH, 'w'))
    if not os.path.exists(def_base_path):
        log_raise(log, "Problem creating configuration file.")

    return


def load_user_config(args, log):
    """Load settings from the user's confiuration file, and add them to `args`.

    Settings are loaded from the configuration file in the user's home
    directory.  Those parameters are added (as attributes) to the `args`
    object.

    Arguments
    ---------
    args : `argparse.Namespace`
        Namespace object to which configuration attributes will be added.

    Returns
    -------
    args : `argparse.Namespace`
        Namespace object with added attributes.

    """
    if not os.path.exists(_CONFIG_PATH):
        err_str = (
            "Configuration file does not exists ({}).\n".format(_CONFIG_PATH) +
            "Run `python -m astrocats setup` to configure.")
        log_raise(log, err_str)

    config = json.load(open(_CONFIG_PATH, 'r'))
    setattr(args, _BASE_PATH_KEY, config[_BASE_PATH_KEY])
    log.debug("Loaded configuration: {}: {}".format(_BASE_PATH_KEY, config[
        _BASE_PATH_KEY]))
    return args


def load_command_line_args(clargs=None):
    """Load and parse command-line arguments.

    Arguments
    ---------
    args : str or None
        'Faked' commandline arguments passed to `argparse`.

    Returns
    -------
    args : `argparse.Namespace` object
        Namespace in which settings are stored - default values modified by the
        given command-line arguments.

    """
    import argparse
    git_vers = get_git()

    parser = argparse.ArgumentParser(
        prog='astrocats',
        description='Generate catalogs for astronomical data.')

    parser.add_argument('command', nargs='?', default=None)

    parser.add_argument(
        '--version',
        action='version',
        version='AstroCats v{}, SHA: {}'.format(__version__, git_vers))
    parser.add_argument(
        '--verbose',
        '-v',
        dest='verbose',
        default=False,
        action='store_true',
        help='Print more messages to the screen.')
    parser.add_argument(
        '--debug',
        '-d',
        dest='debug',
        default=False,
        action='store_true',
        help='Print excessive messages to the screen.')
    parser.add_argument(
        '--include-private',
        dest='private',
        default=False,
        action='store_true',
        help='Include private data in import.')
    parser.add_argument(
        '--travis',
        '-t',
        dest='travis',
        default=False,
        action='store_true',
        help='Run import script in test mode for Travis.')
    parser.add_argument(
        '--clone-depth',
        dest='clone_depth',
        default=0,
        type=int,
        help=('When cloning git repos, only clone out to this depth '
              '(default: 0 = all levels).'))
    parser.add_argument(
        '--purge-outputs',
        dest='purge_outputs',
        default=False,
        action='store_true',
        help=('Purge git outputs after cloning.'))
    parser.add_argument(
        '--log',
        dest='log_filename',
        default=None,
        help='Filename to which to store logging information.')

    # If output files should be written or not
    # ----------------------------------------
    write_group = parser.add_mutually_exclusive_group()
    write_group.add_argument(
        '--write',
        action='store_true',
        dest='write_entries',
        default=True,
        help='Write entries to files [default].')
    write_group.add_argument(
        '--no-write',
        action='store_false',
        dest='write_entries',
        default=True,
        help='do not write entries to file.')

    # If previously cleared output files should be deleted or not
    # -----------------------------------------------------------
    delete_group = parser.add_mutually_exclusive_group()
    delete_group.add_argument(
        '--predelete',
        action='store_true',
        dest='delete_old',
        default=True,
        help='Delete all old event files to begin [default].')
    delete_group.add_argument(
        '--no-predelete',
        action='store_false',
        dest='delete_old',
        default=True,
        help='Do not delete all old event files to start.')

    args, sub_clargs = parser.parse_known_args(args=clargs)
    # Print the help information if no command is given
    if args.command is None:
        parser.print_help()
        return None, None

    return args, sub_clargs


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
    from astrocats.catalog.utils import logger

    # Determine verbosity ('None' means use default)
    log_stream_level = None
    if args.debug:
        log_stream_level = logger.DEBUG
    elif args.verbose:
        log_stream_level = logger.INFO

    # Create log
    log = logger.get_logger(
        stream_level=log_stream_level, tofile=args.log_filename)
    log._verbose = args.verbose
    log._debug = args.debug
    return log


def get_git():
    """Get a string representing the current git status (tag and commit hash).

    Returns
    -------
    git_vers : str
    """
    import subprocess
    git_vers = subprocess.check_output(["git", "describe", "--always"]).strip()
    return git_vers


if __name__ == "__main__":
    main()
