"""Analyze AstroCats Catalogs.
"""
import os
from glob import glob
import numpy as np

_IGNORE_FILES = ['LICENSE', 'README.md']


class Analysis:

    def __init__(self, catalog, log):
        """Initialize `Analysis` instance.

        Arguments
        ---------
        catalog : `astrocats.catalog.catalog.Catalog` subclass
        log : `logging.Logger` object

        """
        self.catalog = catalog
        self.log = log
        return

    def analyze(self, args):
        """Run the analysis routines determined from the given `args`.
        """
        self.log.info("Running catalog analysis")

        if args.count:
            self.count()

        return

    def count(self):
        """Analyze the counts of ...things.
        """
        self.log.info("Running 'count'")
        retvals = {}

        num_tasks = self._count_tasks()
        retvals['num_tasks'] = num_tasks

        num_files = self._count_repo_files()
        retvals['num_files'] = num_files

        return retvals

    def _count_tasks(self):
        """Count the number of tasks, both in the json and directory.
        """
        self.log.warning("Tasks:")
        tasks, task_names = self.catalog._load_task_list_from_file()
        # Total number of all tasks
        num_tasks = len(tasks)
        # Number which are active by default
        num_tasks_act = len([tt for tt, vv in tasks.items() if vv.active])
        # Number of python files in the tasks directory
        num_task_files = os.path.join(self.catalog.PATHS.tasks_dir, '*.py')
        num_task_files = len(glob(num_task_files))
        tasks_str = "{} ({} default active) with {} task-files.".format(
            num_tasks, num_tasks_act, num_task_files)
        self.log.warning(tasks_str)
        return num_tasks

    def _count_repo_files(self):
        """
        """
        self.log.warning("Files:")
        num_files = 0
        repos = self.catalog.PATHS.get_all_repo_folders()
        for rep in repos:
            # Get the last portion of the filepath for this repo
            last_path = _get_last_dirs(rep, 2)
            # Get counts for different file types
            num_json = _count_files_by_type(rep, '*.json')
            num_txt = _count_files_by_type(rep, '*.txt')
            num_all = _count_files_by_type(rep, '*')
            num_oth = num_all - num_json - num_txt
            # Get the number of ignored files
            # (total including ignore, minus 'all')
            num_ign = _count_files_by_type(rep, '*', ignored=True)
            num_ign -= num_all
            f_str = "{}: {} ({} json, {} txt, {} other; {} ignored)".format(
                last_path, num_all, num_json, num_txt, num_oth, num_ign)
            self.log.info(f_str)
            num_files += num_all

        return num_files


def _count_files_by_type(path, suffix, ignored=False):
    """Count files in the given path, with the given pattern.

    If `ignored = True` then files in the `_IGNORE_FILES` list are included.

    Returns
    -------
    num_files : int

    """
    # Get all files matching the given path and pattern
    files = glob(os.path.join(path, suffix))
    # Count the files
    files = [ff for ff in files
             if os.path.split(ff)[-1] not in _IGNORE_FILES or ignored]
    num_files = len(files)
    return num_files


def _get_last_dirs(path, num=1):
    head, tail = os.path.split(path)
    last_path = str(tail)
    for ii in range(num):
        head, tail = os.path.split(head)
        last_path = os.path.join(tail, last_path)

    last_path = "..." + last_path
    return last_path
