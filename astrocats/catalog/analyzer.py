"""Analyze AstroCats Catalogs.
"""
import os
from glob import glob
import numpy as np


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
        """
        """
        self.log.info("Running catalog analysis")

        if args.count:
            self.count()

        return

    def count(self):
        self.log.info("Running 'count'")

        self.log.warning("Tasks:")
        num_tasks = self._count_tasks()
        return

    def _count_tasks(self):
        """
        """
        tasks_list, task_names = self.catalog._load_task_list_from_file()
        # Total number of all tasks
        num_tasks_all = len(tasks_list)
        # Number which are active by default
        num_tasks_act = len([tt for tt, vv in tasks_list.items() if vv.active])
        # Number of python files in the tasks directory
        num_task_files = len(glob(os.path.join(self.catalog.PATHS.tasks_dir, '*.py')))
        tasks_str = "{} ({} default active) with {} task-files.".format(
            num_tasks_all, num_tasks_act, num_task_files)
        self.log.warning(tasks_str)
        return num_tasks_all
