"""Data import methods.
"""

import codecs
import importlib
import json
import os
import resource
import warnings
from collections import OrderedDict

from astrocats import FILENAME, main

from . import task
from .catalog import Catalog
from .utils import pbar, repo_file_list


def import_main(catalog=None):
    """Run all of the import tasks.

    This is executed by the 'scripts.main.py' when the module is run as an
    executable. This can also be run as a method, in which case default
    arguments are loaded, but can be overriden using `**kwargs`.
    """

    # If this is called from `scripts.main`, then `args` will contain
    # parameters. If this is being called as an API function, we need to load
    # default parameters which can then be overwritten below
    if catalog is None:
        from ..supernovae.supernova import Supernova
        warnings.warn("`args` not provided, loading new")
        args = main.load_args(args=['importer'])
        catalog = Catalog(args, Supernova)

    log = catalog.log

    tasks_list = load_task_list(catalog.args, catalog.log)
    warnings.filterwarnings(
        'ignore', r'Warning: converting a masked element to nan.')

    if catalog.args.delete_old:
        log.warning("Deleting all old entry files.")
        catalog.delete_old_entry_files()

    prev_priority = 0
    prev_task_name = ''
    # for task, task_obj in tasks_list.items():
    for task_name, task_obj in tasks_list.items():
        if not task_obj.active:
            continue
        log.warning("Task: '{}'".format(task_name))

        nice_name = task_obj.nice_name
        mod_name = task_obj.module
        func_name = task_obj.function
        priority = task_obj.priority

        # Make sure things are running in the correct order
        if priority < prev_priority:
            raise RuntimeError(("Priority for '{}': '{}', less than prev,"
                                "'{}': '{}'.\n{}").format(
                task_name, priority, prev_task_name, prev_priority, task_obj))

        log.debug("\t{}, {}, {}, {}".format(
            nice_name, priority, mod_name, func_name))
        mod = importlib.import_module('.' + mod_name, package='scripts')
        catalog.current_task = task_obj
        getattr(mod, func_name)(catalog)

        num_events, num_stubs = catalog.count()
        log.warning("Task finished.  Events: {},  Stubs: {}".format(
            num_events, num_stubs))
        catalog.journal_entries()
        num_events, num_stubs = catalog.count()
        log.warning("Journal finished.  Events: {}, Stubs: {}".format(
            num_events, num_stubs))

        prev_priority = priority
        prev_task_name = task_name

    def json_dump(adict, fname):
        json_str = json.dumps(adict, indent='\t', separators=(
            ',', ':'), ensure_ascii=False)
        with codecs.open(fname, 'w', encoding='utf8') as jsf:
            jsf.write(json_str)

    json_dump(bibauthor_dict, FILENAME.BIBAUTHORS)
    json_dump(extinctions_dict, FILENAME.EXTINCT)

    print('Memory used (MBs on Mac, GBs on Linux): ' + '{:,}'.format(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024. / 1024.))
    return


def load_task_list(args, log):
    """Load the list of tasks in the `FILENAME.TASK_LIST` json file.

    A `Task` object is created for each entry, with the parameters filled in.
    These are placed in an OrderedDict, sorted by the `priority` parameter,
    with positive values and then negative values, e.g. [0, 2, 10, -10, -1].
    """

    if args.args_task_list is not None:
        if args.yes_task_list is not None or args.no_task_list is not None:
            raise ValueError(
                "If '--tasks' is used, '--yes' and '--no' shouldnt be.")

    def_task_list_filename = FILENAME.TASK_LIST
    log.debug("Loading task-list from '{}'".format(def_task_list_filename))
    data = json.load(open(def_task_list_filename, 'r'))

    # Make sure 'active' modification lists are all valid
    args_lists = [args.args_task_list, args.yes_task_list, args.no_task_list]
    args_names = ['--tasks', '--yes', '--no']
    for arglist, lname in zip(args_lists, args_names):
        if arglist is not None:
            for tname in arglist:
                if tname not in data.keys():
                    raise ValueError(("Value '{}' in '{}' list does not match"
                                      "any tasks").format(tname, lname))

    tasks = {}
    # `defaults` is a dictionary where each `key` is a task name, and values
    # are its properties
    for key, val in data.items():
        tasks[key] = task.Task(name=key, **val)
        # Modify `active` tasks
        # ---------------------
        # If specific list of tasks is given, make only those active
        if args.args_task_list is not None:
            if key in args.args_task_list:
                tasks[key].active = True
            else:
                tasks[key].active = False
        else:
            # Set 'yes' tasks to *active*
            if args.yes_task_list is not None:
                if key in args.yes_task_list:
                    tasks[key].active = True
            # Set 'no' tasks to *inactive*
            if args.no_task_list is not None:
                if key in args.no_task_list:
                    tasks[key].active = False

    # Sort entries as positive values, then negative values
    #    [0, 1, 2, 2, 10, -100, -10, -1]
    # Tuples are sorted by first element (here: '0' if positive), then second
    # (here normal order)
    tasks = OrderedDict(sorted(tasks.items(), key=lambda t: (
        t[1].priority < 0, t[1].priority, t[1].name)))

    names_act = []
    names_inact = []
    for key, val in tasks.items():
        if val.active:
            names_act.append(key)
        else:
            names_inact.append(key)

    log.info("Active Tasks:\n\t" + ", ".join(nn for nn in names_act))
    log.debug("Inactive Tasks:\n\t" + ", ".join(nn for nn in names_inact))
    return tasks
