"""Handle user arguments when running AstroCats
"""
import argparse


class ArgsHandler:

    def __init__(self, log):
        self.log = log
        parser = self._setup_argparse()
        self.parser = parser
        return

    def run_subcommand(self, args, catalog):
        # Data Import
        # -----------
        if args.subcommand == 'import':
            self.log.info("Running 'import'.")
            catalog.import_data()

        # Git Subcommands
        # ---------------
        elif args.subcommand == 'git-push':
            self.log.info("Running 'git push'.")
            catalog.git_add_commit_push_all_repos()
        elif args.subcommand == 'git-pull':
            self.log.info("Running 'git pull'.")
            catalog.git_add_commit_push_all_repos()
        elif args.subcommand == 'git-reset':
            self.log.info("Running 'git reset'.")
            catalog.git_add_commit_push_all_repos()
        elif args.subcommand == 'git-status':
            self.log.info("Running 'git status'.")
            catalog.git_add_commit_push_all_repos()

        # Analyze Catalogs
        # ----------------
        elif args.subcommand == 'analyze':
            self.log.info("Running 'analyze'.")
            from .analysis import Analysis
            # Create an `Analysis` instance
            lysis = Analysis(catalog, self.log)
            # Pass the command-line arguments to run.
            lysis.analyze(args)

        return

    def load_args(self, args, clargs):
        """Parse arguments and return configuration settings.
        """
        # Parse All Arguments
        args = self.parser.parse_args(args=clargs, namespace=args)

        # Print the help information if no subcommand is given
        # subcommand is required for operation
        if args.subcommand is None:
            self.parser.print_help()
            args = None

        return args

    def _setup_argparse(self):
        """Create `argparse` instance, and setup with appropriate parameters.
        """
        parser = argparse.ArgumentParser(
            prog='catalog', description='Parent Catalog class for astrocats.')

        subparsers = parser.add_subparsers(
            description='valid subcommands', dest='subcommand')

        # Data Import
        # -----------
        # Add the 'import' command, and related arguments
        self._add_parser_arguments_import(subparsers)

        # Git Subcommands
        # ---------------
        # git-push
        self._add_parser_arguments_git_push(subparsers)
        # git-push
        self._add_parser_arguments_git_pull(subparsers)
        # git-push
        self._add_parser_arguments_git_reset(subparsers)
        # git-push
        self._add_parser_arguments_git_status(subparsers)

        # Analyze Catalogs
        # ----------------
        # Add the 'analyze' command, and related arguments
        self._add_parser_arguments_analyze(subparsers)

        return parser

    def _add_parser_arguments_import(self, subparsers):
        """Create parser for 'import' subcommand, and associated arguments.
        """
        import_pars = subparsers.add_parser(
            "import", help="Import data.")

        import_pars.add_argument(
            '--update', '-u', dest='update',
            default=False, action='store_true',
            help='Only update catalog using live sources.')
        import_pars.add_argument(
            '--refresh', '-r', dest='refresh',
            default=False, action='store_true',
            help='Ignore most task caches.')
        import_pars.add_argument(
            '--full-refresh', '-f', dest='full_refresh',
            default=False, action='store_true',
            help='Ignore all task caches.')
        import_pars.add_argument(
            '--archived', '-a', dest='archived',
            default=False, action='store_true',
            help='Always use task caches.')
        import_pars.add_argument(
            '--refresh-list', '-rl', dest='refresh_list',
            default='', nargs='+',
            help='Space-delimited list of caches to clear.')

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

        return import_pars

    def _add_parser_arguments_git_push(self, subparsers):
        """Create a parser for the 'git-push' subcommand.
        """
        push_pars = subparsers.add_parser(
            "git-push",
            help="Add all files to data repositories, commit, and push.")

        return push_pars

    def _add_parser_arguments_git_pull(self, subparsers):
        """Create a parser for the 'git-pull' subcommand.
        """
        pull_pars = subparsers.add_parser(
            "git-pull",
            help="'Pull' all data repositories.")

        return pull_pars

    def _add_parser_arguments_git_reset(self, subparsers):
        """Create a parser for the 'git-reset' subcommand.
        """
        reset_pars = subparsers.add_parser(
            "git-reset",
            help="Hard reset all data repositories.")

        return reset_pars

    def _add_parser_arguments_git_status(self, subparsers):
        """Create a parser for the 'git-status' subcommand.
        """
        status_pars = subparsers.add_parser(
            "git-status",
            help="Get the 'git status' of all data repositories.")

        return status_pars

    def _add_parser_arguments_analyze(self, subparsers):
        """Create a parser for the 'analyze' subcommand.
        """
        lyze_pars = subparsers.add_parser(
            "analyze",
            help="Perform basic analysis on this catalog.")

        lyze_pars.add_argument(
            '--count', '-c', dest='count',
            default=False, action='store_true',
            help='Determine counts of entries, files, etc.')

        return lyze_pars
