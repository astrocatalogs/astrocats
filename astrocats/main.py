"""
"""


def main():
    from datetime import datetime

    beg_time = datetime.now()
    # Process command-line arguments to determine action
    #    If no subcommand (e.g. 'impoter') is given, returns 'None' --> exit
    args = load_args()
    if args is None:
        return

    # FIX
    # LOAD SUPERNOVAE SPECIFIC STUFF EXPLICITLY FOR NOW.  LATER, CHOOSE BASED
    #    ON ARGS WHAT TO IMPORT AND INITIALIZE

    from astrocats.supernovae.catalog import Catalog
    catalog = Catalog(args)
    git_vers = get_git()
    title_str = "Open Supernova Catalog, version: {}".format(git_vers)
    catalog.log.warning("\n\n{}\n{}\n{}\n".format(
        title_str, '=' * len(title_str), beg_time.ctime()))

    # Choose which submodule to run (note: can also use `set_default` with
    # function)
    if args._name == 'sn-import':
        catalog.log.info("Running `importer`.")
        catalog.import_data()

    end_time = datetime.now()
    catalog.log.warning("All complete at {}, After {}".format(
        end_time, end_time - beg_time))
    return


def load_args(args=None):
    import argparse

    parser = argparse.ArgumentParser(
        description='Generate a catalog JSON file and plot HTML files from '
        'SNE data.')
    # parser.add_argument('--foo', action='store_true', help='foo help')
    subparsers = parser.add_subparsers(
        description='valid subcommands', dest='_name',
        help='sub-command help')

    # Build a 'parent' parser whose settings are inhereted by children parsers
    pars_parent = argparse.ArgumentParser(add_help=False)
    pars_parent.add_argument(
        '--verbose', '-v', dest='verbose', default=False, action='store_true',
        help='Print more messages to the screen.')
    pars_parent.add_argument(
        '--debug', '-d', dest='debug', default=False, action='store_true',
        help='Print excessive messages to the screen.')
    pars_parent.add_argument(
        '--travis', '-t',  dest='travis',  default=False, action='store_true',
        help='Run import script in test mode for Travis.')
    pars_parent.add_argument(
        '--log',  dest='log_filename',  default=None,
        help='Filename to which to store logging information.')

    # If output files should be written or not
    # ----------------------------------------
    write_group = pars_parent.add_mutually_exclusive_group()
    write_group.add_argument(
        '--write', action='store_true', dest='write_entries', default=True,
        help='Write entries to files [default].')
    write_group.add_argument(
        '--no-write', action='store_false', dest='write_entries', default=True,
        help='do not write entries to file.')

    # If previously ceared output files should be deleted or not
    # ----------------------------------------------------------
    delete_group = pars_parent.add_mutually_exclusive_group()
    delete_group.add_argument(
        '--predelete', action='store_true', dest='delete_old', default=True,
        help='Delete all old event files to begin [default].')
    delete_group.add_argument(
        '--no-predelete', action='store_false', dest='delete_old',
        default=True, help='Do not delete all old event files to start.')

    # `importer` submodule --- importing supernova data
    pars_imp = subparsers.add_parser("sn-import", parents=[pars_parent],
                                     help="Generate a catalog JSON file")
    pars_imp.add_argument('--update', '-u', dest='update',
                          default=False, action='store_true',
                          help='Only update catalog using live sources.')
    pars_imp.add_argument('--refresh', '-r', dest='refresh',
                          default=False, action='store_true',
                          help='Ignore most task caches.')
    pars_imp.add_argument('--full-refresh', '-f', dest='full_refresh',
                          default=False, action='store_true',
                          help='Ignore all task caches.')
    pars_imp.add_argument('--archived', '-a', dest='archived',
                          default=False, action='store_true',
                          help='Always use task caches.')
    pars_imp.add_argument(
        '--refresh-list', '-rl', dest='refresh_list', default='', nargs='+',
        help='Space-delimited list of caches to clear.')

    # Control which 'tasks' are executed
    # ----------------------------------
    pars_imp.add_argument(
        '--tasks', dest='args_task_list', nargs='*', default=None,
        help='space delimited list of tasks to perform (others disabled).')
    pars_imp.add_argument(
        '--yes', dest='yes_task_list', nargs='+', default=None,
        help='space delimited list of tasks to turn on.')
    pars_imp.add_argument(
        '--no', dest='no_task_list', nargs='+', default=None,
        help='space delimited list of tasks to turn off.')
    pars_imp.add_argument(
        '--min-task-priority', dest='min_task_priority',
        default=None,
        help='minimum priority for a task to run')
    pars_imp.add_argument(
        '--max-task-priority', dest='max_task_priority',
        default=None,
        help='maximum priority for a task to run')
    pars_imp.add_argument(
        '--task-groups', dest='task_groups',
        default=None,
        help='predefined group(s) of tasks to run.')

    args = parser.parse_args(args=args)
    # Print the help information if no subcommand is given
    # subcommand is required for operation
    if args._name is None:
        parser.print_help()
        return None

    return args


def get_git():
    """Get a string representing the current git status --- i.e. tag and commit
    hash.
    """
    import subprocess
    git_vers = subprocess.getoutput(["git describe --always"]).strip()
    return git_vers
