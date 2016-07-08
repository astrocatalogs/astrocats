"""Entry point for the Tidal Disruption Catalog
"""
import argparse


def main(args, clargs, log):
    from .catalog import Catalog

    # Load tidal disruption-specific command-line argumenets
    # (adding them to existing settings)
    args = load_command_line_args(args=args, clargs=clargs)
    if args is None:
        return

    catalog = Catalog(args, log)

    if args.subcommand == 'import':
        log.info("Running `importer`.")
        catalog.import_data()

    return


def load_command_line_args(args=None, clargs=None):
    """Load and parse command-line arguments.
    """

    parser = argparse.ArgumentParser(prog='catalog',
                                     description='Parent Catalog class '
                                     'for astrocats.')

    subparsers = parser.add_subparsers(
        description='valid subcommands', dest='subcommand')
    # `import` --- importing tidal disruption data
    import_pars = subparsers.add_parser("import",
                                        help="Import data.")
    return add_parser_arguments(parser, import_pars, args, clargs)


def add_parser_arguments(parser, import_pars, args, clargs):
    import_pars.add_argument('--update', '-u', dest='update',
                             default=False, action='store_true',
                             help='Only update catalog using live sources.')
    import_pars.add_argument('--refresh', '-r', dest='refresh',
                             default=False, action='store_true',
                             help='Ignore most task caches.')
    import_pars.add_argument('--full-refresh', '-f', dest='full_refresh',
                             default=False, action='store_true',
                             help='Ignore all task caches.')
    import_pars.add_argument('--archived', '-a', dest='archived',
                             default=False, action='store_true',
                             help='Always use task caches.')
    import_pars.add_argument(
        '--refresh-list', '-rl', dest='refresh_list', default='', nargs='+',
        help='Space-delimited list of caches to clear.')

    # Control which 'tasks' are executed
    # ----------------------------------
    import_pars.add_argument(
        '--tasks', dest='args_task_list', nargs='*', default=None,
        help='space delimited list of tasks to perform (others disabled).')
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

    args = parser.parse_args(args=clargs, namespace=args)
    # Print the help information if no subcommand is given
    # subcommand is required for operation
    if args.subcommand is None:
        parser.print_help()
        args = None

    return args
