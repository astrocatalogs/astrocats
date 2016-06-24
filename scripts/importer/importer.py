#!/usr/local/bin/python3.5

import codecs
from collections import OrderedDict
import importlib
import json
import os
import resource
# import sys
import warnings

from scripts import FILENAME
from . import Events
from . funcs import derive_and_sanitize, get_bibauthor_dict, get_extinctions_dict, \
    has_task, name_clean
from . constants import TASK, TRAVIS_QUERY_LIMIT
from .. utils import pbar, repo_file_list, get_logger


def get_old_tasks():
    tasks = OrderedDict([
        ('deleteoldevents', {'nicename': 'Deleting old events',          'update': False}),
        ('internal',        {'nicename': '%pre metadata and photometry', 'update': False}),
        ('radio',           {'nicename': '%pre radio data',              'update': False}),
        ('xray',            {'nicename': '%pre X-ray data',              'update': False}),
        ('simbad',          {'nicename': '%pre SIMBAD',                  'update': False}),
        ('vizier',          {'nicename': '%pre VizieR',                  'update': False}),
        ('donations',       {'nicename': '%pre donations',               'update': False}),
        ('pessto-dr1',      {'nicename': '%pre PESSTO DR1',              'update': False}),
        ('scp',             {'nicename': '%pre SCP',                     'update': False}),
        ('ascii',           {'nicename': '%pre ASCII',                   'update': False}),
        ('cccp',            {'nicename': '%pre CCCP',                    'update': False, 'archived': True}),
        ('suspect_photo',   {'nicename': '%pre SUSPECT',                 'update': False}),
        ('cfa_photo',       {'nicename': '%pre CfA archive photometry',  'update': False}),
        ('ucb_photo',       {'nicename': '%pre UCB photometry',          'update': False, 'archived': True}),
        ('sdss',            {'nicename': '%pre SDSS photometry',         'update': False}),
        ('csp_photo',       {'nicename': '%pre CSP photometry',          'update': False}),
        ('itep',            {'nicename': '%pre ITEP',                    'update': False}),
        ('asiago_photo',    {'nicename': '%pre Asiago metadata',         'update': False}),
        ('tns',             {'nicename': '%pre TNS metadata',            'update': True,  'archived': True}),
        # ('rochester',       {'nicename': '%pre Latest Supernovae',       'update': True,  'archived': False}),
        ('lennarz',         {'nicename': '%pre Lennarz',                 'update': False}),
        ('fermi',           {'nicename': '%pre Fermi',                   'update': False}),
        ('gaia',            {'nicename': '%pre GAIA',                    'update': True,  'archived': False}),
        ('ogle',            {'nicename': '%pre OGLE',                    'update': True,  'archived': False}),
        ('snls',            {'nicename': '%pre SNLS',                    'update': False}),
        ('psthreepi',       {'nicename': '%pre Pan-STARRS 3π',           'update': True,  'archived': False}),
        ('psmds',           {'nicename': '%pre Pan-STARRS MDS',          'update': False}),
        ('crts',            {'nicename': '%pre CRTS',                    'update': True,  'archived': False}),
        ('snhunt',          {'nicename': '%pre SNhunt',                  'update': True,  'archived': False}),
        ('nedd',            {'nicename': '%pre NED-D',                   'update': False}),
        ('cpcs',            {'nicename': '%pre CPCS',                    'update': True,  'archived': False}),
        ('ptf',             {'nicename': '%pre PTF',                     'update': False, 'archived': False}),
        ('des',             {'nicename': '%pre DES',                     'update': False, 'archived': False}),
        ('asassn',          {'nicename': '%pre ASASSN',                  'update': True}),
        ('asiago_spectra',  {'nicename': '%pre Asiago spectra',          'update': True}),
        ('wiserep_spectra', {'nicename': '%pre WISeREP spectra',         'update': False}),
        ('cfa_spectra',     {'nicename': '%pre CfA archive spectra',     'update': False}),
        ('snls_spectra',    {'nicename': '%pre SNLS spectra',            'update': False}),
        ('csp_spectra',     {'nicename': '%pre CSP spectra',             'update': False}),
        ('ucb_spectra',     {'nicename': '%pre UCB spectra',             'update': True,  'archived': True}),
        ('suspect_spectra', {'nicename': '%pre SUSPECT spectra',         'update': False}),
        ('snf_spectra',     {'nicename': '%pre SNH spectra',             'update': False}),
        ('superfit_spectra', {'nicename': '%pre Superfit spectra',        'update': False}),
        ('mergeduplicates', {'nicename': 'Merging duplicates',           'update': False}),
        ('setprefnames',    {'nicename': 'Setting preferred names',      'update': False}),
        ('writeevents',     {'nicename': 'Writing events',               'update': True})
    ])
    return tasks


def import_main(args=None, **kwargs):
    """Run all of the import tasks.

    This is executed by the 'scripts.main.py' when the module is run as an executable.
    This can also be run as a method, in which case default arguments are loaded, but can be
    overriden using `**kwargs`.
    """
    log = get_logger()

    # If this is called from `scripts.main`, then `args` will contain parameters.
    #    If this is being called as an API function, we need to load default parameters which can
    #    then be overwritten below
    if args is None:
        log.debug("`args` not provided, loading new")
        from .. import main
        args = main.load_args(args=['importer'])

    # If this is called as an API function, overwrite variables in `args` with those passed to the
    #    function as keyword arguments.
    for key, val in kwargs.items():
        log.debug("Overriding `args` '{}' = '{}'".format(key, val))
        setattr(args, key, val)
    log.debug("`args` : " + str(vars(args)))

    tasks_list = load_task_list(args, log)
    tasks = get_old_tasks()
    events = OrderedDict()
    # FIX: stubs only need to be loaded for `args.update` ??
    stubs = OrderedDict()
    log.exception("WARNING: not loading stubs for testing!!!")
    # stubs = Events.load_stubs(tasks, args, events, log)
    warnings.filterwarnings('ignore', r'Warning: converting a masked element to nan.')

    prev_priority = 0
    prev_task_name = ''
    # for task, task_obj in tasks_list.items():
    for task_name, task_obj in tasks_list.items():
        if not task_obj.active: continue
        log.info("Task: '{}'".format(task_name))

        nice_name = task_obj.nice_name
        mod_name = task_obj.module
        func_name = task_obj.function
        priority = task_obj.priority
        if priority < prev_priority:
            raise RuntimeError("Priority for '{}': '{}', less than prev, '{}': '{}'.\n{}".format(
                task_name, priority, prev_task_name, prev_priority, task_obj))
        log.debug("\t{}, {}, {}, {}".format(nice_name, priority, mod_name, func_name))
        mod = importlib.import_module('.' + mod_name, package='scripts')
        # events = getattr(mod, func_name)(events, args, tasks, task_obj)
        getattr(mod, func_name)(events, stubs, args, tasks, task_obj, log)
        log.debug("Task finished.  Events: {},  Stubs: {}".format(len(events), len(stubs)))
        events, stubs = Events.journal_events(tasks, args, events, stubs, log)
        log.debug("Journal finished.  Events: {}, Stubs: {}".format(len(events), len(stubs)))

        prev_priority = priority
        prev_task_name = task_name

    return

    files = repo_file_list()

    bibauthor_dict = get_bibauthor_dict()
    extinctions_dict = get_extinctions_dict()

    for ii, fi in enumerate(tq(files, 'Sanitizing and deriving quantities for events')):
        events = OrderedDict()
        name = os.path.basename(os.path.splitext(fi)[0]).replace('.json', '')
        events, name = Events.add_event(tasks, args, events, name, log, load_stubs_if_empty=False)
        events, extinctions_dict, bibauthor_dict = derive_and_sanitize(
            tasks, args, events, extinctions_dict, bibauthor_dict, nedd_dict)
        if has_task(tasks, args, 'writeevents'):
            Events.write_all_events(events, args, empty=True, gz=True, bury=True)
        if args.travis and ii > TRAVIS_QUERY_LIMIT:
            break

    def json_dump(adict, fname):
        json_str = json.dumps(adict, indent='\t', separators=(',', ':'), ensure_ascii=False)
        with codecs.open(fname, 'w', encoding='utf8') as jsf:
            jsf.write(json_str)

    BIBAUTHORS_FILENAME = '../bibauthors.json'
    EXTINCTIONS_FILENAME = '../extinctions.json'
    json_dump(bibauthor_dict, BIBAUTHORS_FILENAME)
    json_dump(extinctions_dict, EXTINCTIONS_FILENAME)

    print('Memory used (MBs on Mac, GBs on Linux): ' + '{:,}'.format(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024./1024.))
    return


def do_nedd(events, stubs, args, tasks, task_obj, log):
    import csv
    from html import unescape
    from . constants import PATH
    from . funcs import uniq_cdl
    from .. utils import is_number, Decimal
    nedd_path = os.path.join(PATH.REPO_EXTERNAL, 'NED25.12.1-D-10.4.0-20151123.csv')
    with open(nedd_path, 'r') as f:
        data = sorted(list(csv.reader(f, delimiter=',', quotechar='"'))[13:], key=lambda x: (x[9], x[3]))
        reference = "NED-D"
        refurl = "http://ned.ipac.caltech.edu/Library/Distances/"
        nedddict = OrderedDict()
        olddistname = ''
        for r, row in enumerate(tq(data, currenttask = currenttask)):
            if r <= 12:
                continue
            distname = row[3]
            name = name_clean(distname)
            distmod = row[4]
            moderr = row[5]
            dist = row[6]
            bibcode = unescape(row[8])
            snname = name_clean(row[9])
            redshift = row[10]
            cleanhost = ''

            if name != snname and (name + ' HOST' != snname):
                cleanhost = host_clean(distname)
                if cleanhost.endswith(' HOST'):
                    cleanhost = ''
                if not is_number(dist):
                    print(dist)
                if dist:
                    nedddict.setdefault(cleanhost,[]).append(Decimal(dist))

            if snname and 'HOST' not in snname:
                events, snname, secondarysource = Events.new_event(
                    tasks, args, events, snname, log,
                    refname = reference, url = refurl, secondary = True)
                if bibcode:
                    source = events[snname].add_source(bibcode = bibcode)
                    sources = uniq_cdl([source, secondarysource])
                else:
                    sources = secondarysource

                if name == snname:
                    if redshift:
                        events[snname].add_quantity('redshift', redshift, sources)
                    if dist:
                        events[snname].add_quantity('comovingdist', dist, sources)
                        if not redshift:
                            try:
                                redshift = pretty_num(z_at_value(cosmo.comoving_distance, float(dist) * un.Mpc, zmax = 5.0), sig = get_sig_digits(str(dist)))
                            except (KeyboardInterrupt, SystemExit):
                                raise
                            except:
                                pass
                            else:
                                cosmosource = events[name].add_source(bibcode='2015arXiv150201589P')
                                events[snname].add_quantity('redshift', redshift, uniq_cdl(sources.split(',') + [cosmosource]))

                if cleanhost:
                    events[snname].add_quantity('host', cleanhost, sources)

                if args.update and olddistname != distname:
                    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
            olddistname = distname

        events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    return events

def delete_old_event_files(*args):
    # Delete all old event JSON files
    repo_files = repo_file_list()
    for rfil in pbar(repo_files, desc='Deleting old events'):
        os.remove(rfil)
    return


def load_task_list(args, log):
    """Load the list of tasks in the `FILENAME.TASK_LIST` json file.

    A `TASK` object is created for each entry, with the parameters filled in.
    These are placed in an OrderedDict, sorted by the `priority` parameter, with positive values
    and then negative values, e.g. [0, 2, 10, -10, -1].
    """

    # print("refresh_list = ", args.refresh_list)
    # sys.exit(3189752)

    if args.args_task_list is not None:
        if args.yes_task_list is not None or args.no_task_list is not None:
            raise ValueError("If '--tasks' is used, '--yes' and '--no' shouldnt be.")

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
                    raise ValueError("Value '{}' in '{}' list does not match any tasks".format(
                        tname, lname))

    tasks = {}
    # `defaults` is a dictionary where each `key` is a task name, and values are its properties
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
    #    Tuples are sorted by first element (here: '0' if positive), then second (here normal order)
    tasks = OrderedDict(sorted(tasks.items(), key=lambda t: (t[1].priority < 0, t[1].priority)))

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
