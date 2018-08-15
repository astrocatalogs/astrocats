"""Handle user arguments when running AstroCats."""

import argparse

import astrocats


def parse_args():

    # Setup argparse instance
    # ---------------------------------------

    # desc = 'Generate catalogs for astronomical data.'
    desc = None

    parser = argparse.ArgumentParser(prog='astrocats', description=desc)

    # Add positional arguments
    parser.add_argument('catalog', help='path to catalog')

    # Add subparsers for particular commands
    subparsers = parser.add_subparsers(dest='command')

    # Add generally applicable arguments
    parser = _add_general(parser)

    # Add subparser arguments for particular commands
    # ---------------------------------------------------------
    import_parser = subparsers.add_parser('import')
    import_parser = _add_import(import_parser)

    produce_parser = subparsers.add_parser('produce')
    produce_parser = _add_produce(produce_parser)

    git_parser = subparsers.add_parser('git')
    git_parser = _add_git(git_parser)

    analyze_parser = subparsers.add_parser('analyze')
    analyze_parser = _add_analyze(analyze_parser)

    # Parse Arguments
    # ----------------------------
    args = parser.parse_args()

    # print help if no command is given
    if args.command is None:
        print("astrocats requires a command to be given.\n")
        parser.print_help()
        return None

    return args


def _add_general(parser):
    version_info = 'AstroCats v{}, SHA: {}'.format(
        astrocats.__version__, astrocats.__git_version__)

    parser.add_argument('--version', action='version', version=version_info)
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
        '--log',
        dest='log_filename',
        default=None,
        help='Filename to which to store logging information.')
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

    return parser


def _add_import(import_pars):
    import_pars.add_argument(
        '--update', '-u', dest='update',
        default=False, action='store_true',
        help='Only update catalog using live sources.')
    import_pars.add_argument(
        '--load-stubs', dest='load_stubs',
        default=False, action='store_true',
        help='Load stubs before running.')
    import_pars.add_argument(
        '--archived', '-a', dest='archived',
        default=False, action='store_true',
        help='Always use task caches.')
    import_pars.add_argument(
        '--include-private',
        dest='private',
        default=False,
        action='store_true',
        help='Include private data in import.')

    # Control which 'tasks' are executed
    # ----------------------------------
    import_pars.add_argument(
        '--tasks', dest='args_task_list', nargs='*', default=None,
        help='space delimited list of tasks to perform.')
    import_pars.add_argument(
        '--yes', dest='yes_task_list', nargs='+', default=None,
        help='space delimited list of tasks to turn on.')
    import_pars.add_argument(
        '--no', dest='no_task_list', nargs='+', default=None,
        help='space delimited list of tasks to turn off.')
    import_pars.add_argument(
        '--min-task-priority', dest='min_task_priority',
        default=None,
        help='minimum priority for a task to run')
    import_pars.add_argument(
        '--max-task-priority', dest='max_task_priority',
        default=None,
        help='maximum priority for a task to run')
    import_pars.add_argument(
        '--task-groups', dest='task_groups',
        default=None,
        help='predefined group(s) of tasks to run.')

    # If previously cleared output files should be deleted or not
    # -----------------------------------------------------------
    delete_group = import_pars.add_mutually_exclusive_group()
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

    return import_pars


def _add_produce(produce_pars):
    produce_pars.add_argument(
        '--event-list',
        '-el',
        dest='event_list',
        help='Process (only) a list of target events',
        default=None,
        type=str,
        nargs='+')

    produce_pars.add_argument(
        '--boneyard',
        '-by',
        dest='boneyard',
        help='Make "boneyard" catalog',
        default=False,
        action='store_true')

    produce_pars.add_argument(
        '--authors',
        dest='authors',
        help='Query ADS for all-authors information (slow!)',
        default=False,
        action='store_true')

    produce_pars.add_argument(
        '--hosts',
        dest='hosts',
        help='Query skyservice or skyview for host-images (slow!)',
        default=False,
        action='store_true')

    return produce_pars


def _add_git(git_pars):

    git_group = git_pars.add_mutually_exclusive_group()

    git_group.add_argument(
        "--clone", dest="git_clone", action='store_true', default=False,
        help="Clone all defined data repositories if they dont exist.")

    git_group.add_argument(
        "--push", dest="git_push", action='store_true', default=False,
        help="Add all files to data repositories, commit, and push.")

    git_group.add_argument(
        "--pull", dest="git_pull", action='store_true', default=False,
        help="'Pull' all data repositories.")

    git_group.add_argument(
        "--reset-local", dest="git_reset_local", action='store_true', default=False,
        help="Hard reset all data repositories using local 'HEAD'.")

    git_group.add_argument(
        "--reset-origin", dest="git_reset_origin", action='store_true', default=False,
        help="Hard reset all data repositories using 'origin/master'.")

    git_group.add_argument(
        "--status", dest="git_status", action='store_true', default=False,
        help="Get the 'git status' of all data repositories.")

    return git_pars


def _add_analyze(analyze_pars):
    analyze_pars.add_argument(
        '--count', '-c', dest='count_flag',
        default=False, action='store_true',
        help='Determine counts of entries, files, etc.')

    analyze_pars.add_argument(
        '--data-tree', dest='data_tree_flag',
        default=True, action='store_true',
        help='Determine what entries, quantities, parameters are stored in the catalog.')

    return analyze_pars
