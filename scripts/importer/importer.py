#!/usr/local/bin/python3.5

import codecs
import importlib
import json
import os
import resource
import warnings
from collections import OrderedDict

from scripts import FILENAME

from . import Events
from ..utils import get_logger, pbar, repo_file_list
from ..Catalog import Catalog
from .constants import TASK, TRAVIS_QUERY_LIMIT
from .funcs import derive_and_sanitize, has_task


def get_old_tasks():
    tasks = OrderedDict([
        ('deleteoldevents',  {'nicename': 'Deleting old events',
                              'update': False}),
        ('internal',         {'nicename': '%pre metadata and photometry',
                              'update': False}),
        ('radio',            {'nicename': '%pre radio data',
                              'update': False}),
        ('xray',             {'nicename': '%pre X-ray data',
                              'update': False}),
        ('simbad',           {'nicename': '%pre SIMBAD',
                              'update': False}),
        ('vizier',           {'nicename': '%pre VizieR',
                              'update': False}),
        ('donations',        {'nicename': '%pre donations',
                              'update': False}),
        ('pessto-dr1',       {'nicename': '%pre PESSTO DR1',
                              'update': False}),
        ('scp',              {'nicename': '%pre SCP',
                              'update': False}),
        ('ascii',            {'nicename': '%pre ASCII',
                              'update': False}),
        ('cccp',             {'nicename': '%pre CCCP',
                              'update': False, 'archived': True}),
        ('suspect_photo',    {'nicename': '%pre SUSPECT',
                              'update': False}),
        ('cfa_photo',        {'nicename': '%pre CfA archive photometry',
                              'update': False}),
        ('ucb_photo',        {'nicename': '%pre UCB photometry',
                              'update': False, 'archived': True}),
        ('sdss',             {'nicename': '%pre SDSS photometry',
                              'update': False}),
        ('csp_photo',        {'nicename': '%pre CSP photometry',
                              'update': False}),
        ('itep',             {'nicename': '%pre ITEP',
                              'update': False}),
        ('asiago_photo',     {'nicename': '%pre Asiago metadata',
                              'update': False}),
        ('tns',              {'nicename': '%pre TNS metadata',
                              'update': True,  'archived': True}),
        # ('rochester',       {'nicename': '%pre Latest Supernovae',
        #                      'update': True,  'archived': False}),
        ('lennarz',          {'nicename': '%pre Lennarz',
                              'update': False}),
        ('fermi',            {'nicename': '%pre Fermi',
                              'update': False}),
        ('gaia',             {'nicename': '%pre GAIA',
                              'update': True,  'archived': False}),
        ('ogle',             {'nicename': '%pre OGLE',
                              'update': True,  'archived': False}),
        ('snls',             {'nicename': '%pre SNLS',
                              'update': False}),
        ('psthreepi',        {'nicename': '%pre Pan-STARRS 3π',
                              'update': True, 'archived': False}),
        ('psmds',            {'nicename': '%pre Pan-STARRS MDS',
                              'update': False}),
        ('crts',             {'nicename': '%pre CRTS',
                              'update': True,  'archived': False}),
        ('snhunt',           {'nicename': '%pre SNhunt',
                              'update': True,  'archived': False}),
        ('nedd',             {'nicename': '%pre NED-D',
                              'update': False}),
        ('cpcs',             {'nicename': '%pre CPCS',
                              'update': True,  'archived': False}),
        ('ptf',              {'nicename': '%pre PTF',
                              'update': False, 'archived': False}),
        ('des',              {'nicename': '%pre DES',
                              'update': False, 'archived': False}),
        ('asassn',           {'nicename': '%pre ASASSN',
                              'update': True}),
        ('asiago_spectra',   {'nicename': '%pre Asiago spectra',
                              'update': True}),
        ('wiserep_spectra',  {'nicename': '%pre WISeREP spectra',
                              'update': False}),
        ('cfa_spectra',      {'nicename': '%pre CfA archive spectra',
                              'update': False}),
        ('snls_spectra',     {'nicename': '%pre SNLS spectra',
                              'update': False}),
        ('csp_spectra',      {'nicename': '%pre CSP spectra',
                              'update': False}),
        ('ucb_spectra',      {'nicename': '%pre UCB spectra',
                              'update': True,  'archived': True}),
        ('suspect_spectra',  {'nicename': '%pre SUSPECT spectra',
                              'update': False}),
        ('snf_spectra',      {'nicename': '%pre SNH spectra',
                              'update': False}),
        ('superfit_spectra', {'nicename': '%pre Superfit spectra',
                              'update': False}),
        ('mergeduplicates',  {'nicename': 'Merging duplicates',
                              'update': False}),
        ('setprefnames',     {'nicename': 'Setting preferred names',
                              'update': False}),
        ('writeevents',      {'nicename': 'Writing events',
                              'update': True})
    ])
    return tasks


def import_main(args=None, **kwargs):
    """Run all of the import tasks.

    This is executed by the 'scripts.main.py' when the module is run as an
    executable. This can also be run as a method, in which case default
    arguments are loaded, but can be overriden using `**kwargs`.
    """
    log = get_logger()

    catalog = Catalog()

    # If this is called from `scripts.main`, then `args` will contain
    # parameters. If this is being called as an API function, we need to load
    # default parameters which can then be overwritten below
    if args is None:
        log.debug("`args` not provided, loading new")
        from .. import main
        args = main.load_args(args=['importer'])

    # If this is called as an API function, overwrite variables in `args` with
    # those passed to the function as keyword arguments.
    for key, val in kwargs.items():
        log.debug("Overriding `args` '{}' = '{}'".format(key, val))
        setattr(args, key, val)
    log.debug("`args` : " + str(vars(args)))

    tasks_list = load_task_list(args, log)
    tasks = get_old_tasks()
    events = OrderedDict()
    warnings.filterwarnings(
        'ignore', r'Warning: converting a masked element to nan.')

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
        # events = getattr(mod, func_name)(events, args, tasks, task_obj)
        events = getattr(mod, func_name)(events, args, tasks, task_obj, log)
        num_events, num_stubs = Events.count(events)
        log.warning("Task finished.  Events: {},  Stubs: {}".format(
            num_events, num_stubs))
        events = Events.journal_events(tasks, args, events, log)
        log.warning("Journal finished.  Events: {}, Stubs: {}".format(
            num_events, num_stubs))

        prev_priority = priority
        prev_task_name = task_name

    return

    files = repo_file_list()

    for ii, fi in enumerate(
            pbar(files, 'Sanitizing and deriving quantities for events')):
        events = OrderedDict()
        name = os.path.basename(os.path.splitext(fi)[0]).replace('.json', '')
        events, name = Events.add_event(
            tasks, args, events, name, log, load_stubs_if_empty=False)
        events, extinctions_dict, bibauthor_dict = derive_and_sanitize(
            catalog, tasks, args, events)
        if has_task(tasks, args, 'writeevents'):
            Events.write_all_events(
                events, args, empty=True, gz=True, bury=True)
        if args.travis and ii > TRAVIS_QUERY_LIMIT:
            break

    def json_dump(adict, fname):
        json_str = json.dumps(adict, indent='\t', separators=(
            ',', ':'), ensure_ascii=False)
        with codecs.open(fname, 'w', encoding='utf8') as jsf:
            jsf.write(json_str)

    BIBAUTHORS_FILENAME = '../bibauthors.json'
    EXTINCTIONS_FILENAME = '../extinctions.json'
    json_dump(bibauthor_dict, BIBAUTHORS_FILENAME)
    json_dump(extinctions_dict, EXTINCTIONS_FILENAME)

    print('Memory used (MBs on Mac, GBs on Linux): ' + '{:,}'.format(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024. / 1024.))
    return


def delete_old_event_files(events, args, tasks, task_obj, log):
    if len(events):
        err_str = "`delete_old_event_files` with `events` not empty!"
        log.error(err_str)
        raise RuntimeError(err_str)
    # Delete all old event JSON files
    repo_files = repo_file_list()
    for rfil in pbar(repo_files, desc='Deleting old events'):
        os.remove(rfil)
    return events


def load_task_list(args, log):
    """Load the list of tasks in the `FILENAME.TASK_LIST` json file.

    A `TASK` object is created for each entry, with the parameters filled in.
    These are placed in an OrderedDict, sorted by the `priority` parameter,
    with positive values and then negative values, e.g. [0, 2, 10, -10, -1].
    """

    # print("refresh_list = ", args.refresh_list)
    # sys.exit(3189752)

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
        tasks[key] = TASK(name=key, **val)
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
