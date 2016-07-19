"""Analyze AstroCats Catalogs.

To-Do
-----
* Number of entries with ___ (e.g. spectra)
* Number of ___ (e.g. spectra) etc.


"""
import os
from glob import glob
import numpy as np


class Analysis:

    _IGNORE_FILES = ['LICENSE', 'README.md']
    # If no specific types should be counted, make this an empty list
    _COUNT_FILE_TYPES = ['json', 'txt']

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

        Returns
        -------
        retvals : dict
            Dictionary of 'property-name: counts' pairs for further processing

        """
        self.log.info("Running 'count'")
        retvals = {}

        # Numbers of 'tasks'
        num_tasks = self._count_tasks()
        retvals['num_tasks'] = num_tasks

        # Numbers of 'files'
        num_files = self._count_repo_files()
        retvals['num_files'] = num_files

        return retvals

    def _count_tasks(self):
        """Count the number of tasks, both in the json and directory.

        Returns
        -------
        num_tasks : int
            The total number of all tasks included in the `tasks.json` file.

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
        """Count the number of files in the data repositories.

        `_COUNT_FILE_TYPES` are used to determine which file types are checked
        explicitly.
        `_IGNORE_FILES` determine which files are ignored in (most) counts.

        Returns
        -------
        repo_files : int
            Total number of (non-ignored) files in all data repositories.

        """
        self.log.warning("Files:")
        num_files = 0
        repos = self.catalog.PATHS.get_all_repo_folders()
        num_type = np.zeros(len(self._COUNT_FILE_TYPES), dtype=int)
        num_ign = 0
        for rep in repos:
            # Get the last portion of the filepath for this repo
            last_path = _get_last_dirs(rep, 2)
            # Get counts for different file types
            n_all = self._count_files_by_type(rep, '*')
            n_type = np.zeros(len(self._COUNT_FILE_TYPES), dtype=int)
            for ii, ftype in enumerate(self._COUNT_FILE_TYPES):
                n_type[ii] = self._count_files_by_type(rep, '*.' + ftype)
            # Get the number of ignored files
            # (total including ignore, minus 'all')
            n_ign = self._count_files_by_type(rep, '*', ignore=False)
            n_ign -= n_all
            f_str = self._file_nums_str(n_all, n_type, n_ign)
            f_str = "{}: {}".format(last_path, f_str)
            self.log.warning(f_str)
            # Update cumulative counts
            num_files += n_all
            num_type += n_type
            num_ign += n_ign

        f_str = self._file_nums_str(num_files, num_type, num_ign)
        self.log.warning(f_str)
        return num_files

    def _file_nums_str(self, n_all, n_type, n_ign):
        """Construct a string showing the number of different file types.

        Returns
        -------
        f_str : str
        """
        # 'other' is the difference between all and named
        n_oth = n_all - np.sum(n_type)

        f_str = "{} Files".format(n_all) + " ("
        if len(n_type):
            f_str += ", ".join("{} {}".format(name, num) for name, num in
                               zip(self._COUNT_FILE_TYPES, n_type))
            f_str += ", "
        f_str += "other {}; {} ignored)".format(n_oth, n_ign)
        return f_str

    def _count_files_by_type(self, path, pattern, ignore=True):
        """Count files in the given path, with the given pattern.

        If `ignore = True`, skip files in the `_IGNORE_FILES` list.

        Returns
        -------
        num_files : int

        """
        # Get all files matching the given path and pattern
        files = glob(os.path.join(path, pattern))
        # Count the files
        files = [ff for ff in files
                 if os.path.split(ff)[-1] not in self._IGNORE_FILES
                 or not ignore]
        num_files = len(files)
        return num_files


def _get_last_dirs(path, num=1):
    """Get a path including only the trailing `num` directories.

    Returns
    -------
    last_path : str

    """
    head, tail = os.path.split(path)
    last_path = str(tail)
    for ii in range(num):
        head, tail = os.path.split(head)
        last_path = os.path.join(tail, last_path)

    last_path = "..." + last_path
    return last_path
