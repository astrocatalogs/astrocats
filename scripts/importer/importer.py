#!/usr/local/bin/python3.5

import codecs
from collections import OrderedDict
import importlib
import json
import os
import resource
import warnings

from .. utils import pbar, repo_file_list
from . funcs import add_event, derive_and_sanitize, \
    do_task, get_bibauthor_dict, get_extinctions_dict, \
    has_task
from scripts import FILENAME
from . constants import TASK


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
    # If this is called from `scripts.main`, then `args` will contain parameters.
    #    If this is being called as an API function, we need to load default parameters which can
    #    then be overwritten below
    if args is None:
        from .. import main
        args = main.load_args(args=['importer'])

    # If this is called as an API function, overwrite variables in `args` with those passed to the
    #    function as keyword arguments.
    for key, val in kwargs.items():
        setattr(args, key, val)
    print(vars(args))

    tasks_list = load_task_list()
    tasks = get_old_tasks()
    events = OrderedDict()
    warnings.filterwarnings('ignore', r'Warning: converting a masked element to nan.')

    prev_priority = 0
    prev_task_name = ''
    # for task, task_obj in tasks_list.items():
    for task_name, task_obj in tasks_list.items():
        print("\n", task_name)

        nice_name = task_obj.nice_name
        mod_name = task_obj.module
        func_name = task_obj.function
        priority = task_obj.priority
        if priority < prev_priority:
            raise RuntimeError("Priority for '{}': '{}', less than prev, '{}': '{}'.\n{}".format(
                task_name, priority, prev_task_name, prev_priority, task_obj))
        print("\t{}, {}, {}, {}".format(nice_name, priority, mod_name, func_name))
        mod = importlib.import_module('.' + mod_name, package='scripts')
        # events = getattr(mod, func_name)(events, args, tasks)
        getattr(mod, func_name)(events, args, tasks)

        prev_priority = priority
        prev_task_name = task_name

        # break
        continue

        if do_task(tasks, args, task, 'deleteoldevents'):
            # Delete `current_task` here and wrap deletion in `pbar_strings` ??
            delete_old_event_files()

        # ========================================================================
        # ========================================================================
        # ========================================================================
        # ========================================================================

        if do_task(tasks, args, task, 'internal'):
            from mtasks.general_data import do_internal
            events = do_internal(events, args, tasks)

        if do_task(tasks, args, task, 'radio'):
            from mtasks.general_data import do_external_radio
            events = do_external_radio(events, args, tasks)

        if do_task(tasks, args, task, 'xray'):
            from mtasks.general_data import do_external_xray
            events = do_external_xray(events, args, tasks)

        # if do_task(tasks, args, task, 'simbad'):
        #     from mtasks import do_simbad
        #     events = do_simbad(events, args, tasks)

        if do_task(tasks, args, task, 'donations'):
            from mtasks.donations import do_donations
            events = do_donations(events, args, tasks)

        if do_task(tasks, args, task, 'ascii'):
            from mtasks.general_data import do_ascii
            events = do_ascii(events, args, tasks)

        if do_task(tasks, args, task, 'cccp'):
            from mtasks.general_data import do_cccp
            events = do_cccp(events, args, tasks)

        # Suspect catalog
        if do_task(tasks, args, task, 'suspect_photo'):
            from mtasks.suspect import do_suspect_photo
            events = do_suspect_photo(events, args, tasks)

        if do_task(tasks, args, task, 'suspect_spectra'):
            from mtasks.suspect import do_suspect_spectra
            events = do_suspect_spectra(events, args, tasks)

        # CfA data
        if do_task(tasks, args, task, 'cfa_photo'):
            from mtasks.cfa import do_cfa_photo
            events = do_cfa_photo(events, args, tasks)

        if do_task(tasks, args, task, 'cfa_spectra'):
            from mtasks.cfa import do_cfa_spectra
            events = do_cfa_spectra(events, args, tasks)

        if do_task(tasks, args, task, 'wiserep_spectra'):
            from mtasks.wiserep import do_wiserep_spectra
            events = do_wiserep_spectra(events, args, tasks)

        # Import primary data sources from Vizier
        if do_task(tasks, args, task, 'vizier'):
            from mtasks.vizier import do_vizier
            events = do_vizier(events, args, tasks)

        if do_task(tasks, args, task, 'lennarz'):
            from mtasks.vizier import do_lennarz
            events = do_lennarz(events, args, tasks)

        if do_task(tasks, args, task, 'snls_spectra'):
            from mtasks.vizier import do_snls_spectra
            events = do_snls_spectra(events, args, tasks)

        if do_task(tasks, args, task, 'pessto-dr1'):
            from mtasks.general_data import do_pessto
            events = do_pessto(events, args, tasks)

        if do_task(tasks, args, task, 'scp'):
            from mtasks.general_data import do_scp
            events = do_scp(events, args, tasks)

        # New UCB import
        if do_task(tasks, args, task, 'ucb_photo'):
            from mtasks.berkeley import do_ucb_photo
            events = do_ucb_photo(events, args, tasks)

        if do_task(tasks, args, task, 'ucb_spectra'):
            from mtasks.berkeley import do_ucb_spectra
            events = do_ucb_spectra(events, args, tasks)

        # Import SDSS
        if do_task(tasks, args, task, 'sdss'):
            from mtasks.general_data import do_sdss
            events = do_sdss(events, args, tasks)

        # Import GAIA
        if do_task(tasks, args, task, 'gaia'):
            from mtasks.general_data import do_gaia
            events = do_gaia(events, args, tasks)

        # Import CSP
        # VizieR catalogs exist for this: J/AJ/139/519, J/AJ/142/156. Should replace eventually.
        if do_task(tasks, args, task, 'csp_photo'):
            from mtasks.carnegie import do_csp_photo
            events = do_csp_photo(events, args, tasks)

        if do_task(tasks, args, task, 'csp_spectra'):
            from mtasks.carnegie import do_csp_spectra
            events = do_csp_spectra(events, args, tasks)

        if do_task(tasks, args, task, 'crts'):
            from mtasks.general_data import do_crts
            events = do_crts(events, args, tasks)

        if do_task(tasks, args, task, 'rochester'):
            from mtasks.rochester import do_rochester
            events = do_rochester(events, args, tasks)

        if do_task(tasks, args, task, 'ogle'):
            from mtasks.ogle import do_ogle
            events = do_ogle(events, args, tasks)

        # Now import the Asiago catalog
        if do_task(tasks, args, task, 'asiago_photo'):
            from mtasks.asiago import do_asiago_photo
            events = do_asiago_photo(events, args, tasks)

        if do_task(tasks, args, task, 'asiago_spectra'):
            from mtasks.asiago import do_asiago_spectra
            events = do_asiago_spectra(events, args, tasks)

        if do_task(tasks, args, task, 'snls'):
            from mtasks.general_data import do_snls
            events = do_snls(events, args, tasks)

        if do_task(tasks, args, task, 'psmds'):
            from mtasks.panstarrs import do_ps_mds
            events = do_ps_mds(events, args, tasks)

        if do_task(tasks, args, task, 'psthreepi'):
            from mtasks.panstarrs import do_ps_threepi
            events = do_ps_threepi(events, args, tasks)

        if do_task(tasks, args, task, 'fermi'):
            from mtasks.general_data import do_fermi
            events = do_fermi(events, args, tasks)

        if do_task(tasks, args, task, 'ptf'):
            from mtasks.palomar import do_ptf
            events = do_ptf(events, args, tasks)

        # Import ITEP
        if do_task(tasks, args, task, 'itep'):
            from mtasks.general_data import do_itep
            events = do_itep(events, args, tasks)

        if do_task(tasks, args, task, 'tns'):
            from mtasks.general_data import do_tns
            events = do_tns(events, args, tasks)

        if do_task(tasks, args, task, 'snhunt'):
            from mtasks.general_data import do_snhunt
            events = do_snhunt(events, args, tasks)

        if do_task(tasks, args, task, 'snf_spectra'):
            from mtasks.snfactory import do_snf_specta
            events = do_snf_specta(events, args, tasks)

        if do_task(tasks, args, task, 'asassn'):
            from mtasks.asassn import do_asassn
            events = do_asassn(events, args, tasks)

        if do_task(tasks, args, task, 'nedd'):
            from mtasks.general_data import do_nedd
            events = do_nedd(events, args, tasks)

        # Import CPCS
        if do_task(tasks, args, task, 'cpcs'):
            from mtasks.general_data import do_cpcs
            events = do_cpcs(events, args, tasks)

        if do_task(tasks, args, task, 'des'):
            from mtasks.general_data import do_des
            events = do_des(events, args, tasks)

        if do_task(tasks, args, task, 'superfit_spectra'):
            from mtasks.general_data import do_superfit_spectra
            events = do_superfit_spectra(events, args, tasks)

        # ========================================================================
        # ========================================================================
        # ========================================================================
        # ========================================================================

        if do_task(tasks, args, task, 'mergeduplicates'):
            from . funcs import merge_duplicates
            events = merge_duplicates(tasks, args, events)

        if do_task(tasks, args, task, 'setprefnames'):
            from . funcs import set_preferred_names
            events = set_preferred_names(tasks, args, events)


    return


    files = repo_file_list()

    bibauthor_dict = get_bibauthor_dict()
    extinctions_dict = get_extinctions_dict()

    for fi in pbar(files, 'Sanitizing and deriving quantities for events'):
        events = OrderedDict()
        name = os.path.basename(os.path.splitext(fi)[0]).replace('.json', '')
        name = add_event(tasks, args, events, args, name, loadifempty=False)
        events, extinctions_dict, bibauthor_dict = derive_and_sanitize(
            tasks, args, events, extinctions_dict, bibauthor_dict, nedd_dict)
        if has_task(tasks, args, 'writeevents'):
            write_all_events(events, args, empty=True, gz=True, bury=True)

    jsonstring = json.dumps(bibauthor_dict, indent='\t', separators=(',', ':'), ensure_ascii=False)
    with codecs.open('../bibauthors.json', 'w', encoding='utf8') as f:
        f.write(jsonstring)
    jsonstring = json.dumps(extinctions_dict, indent='\t', separators=(',', ':'), ensure_ascii=False)
    with codecs.open('../extinctions.json', 'w', encoding='utf8') as f:
        f.write(jsonstring)

    print('Memory used (MBs on Mac, GBs on Linux): ' + '{:,}'.format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024./1024.))
    return


def do_nedd(events, args, tasks):
    import csv
    from html import unescape
    from . constants import PATH
    from . funcs import add_quantity, add_source, journal_events, uniq_cdl
    from .. utils import is_number, Decimal
    nedd_path = os.path.join(PATH.REPO_EXTERNAL, 'NED25.12.1-D-10.4.0-20151123.csv')
    data = csv.reader(open(nedd_path, 'r'), delimiter=',', quotechar='"')
    reference = 'NED-D'
    refurl = 'http://ned.ipac.caltech.edu/Library/Distances/'
    nedd_dict = OrderedDict()
    oldhostname = ''
    for r, row in enumerate(data):
        if r <= 12:
            continue
        hostname = row[3]
        if args.update and oldhostname != hostname:
            events = journal_events(tasks, args, events)
        # distmod = row[4]
        # moderr = row[5]
        dist = row[6]
        bibcode = unescape(row[8])
        name = ''
        if hostname.startswith('SN '):
            if is_number(hostname[3:7]):
                name = 'SN' + hostname[3:]
            else:
                name = hostname[3:]
        elif hostname.startswith('SNLS '):
            name = 'SNLS-' + hostname[5:].split()[0]
        else:
            cleanhost = hostname.replace('MESSIER 0', 'M').replace('MESSIER ', 'M').strip()
            if True in [x in cleanhost for x in ['UGC', 'PGC', 'IC']]:
                cleanhost = ' '.join([x.lstrip('0') for x in cleanhost.split()])
            if 'ESO' in cleanhost:
                cleanhost = cleanhost.replace(' ', '').replace('ESO', 'ESO ')
            nedd_dict.setdefault(cleanhost, []).append(Decimal(dist))

        if name:
            name = add_event(tasks, args, events, name)
            sec_source = add_source(events, name, refname=reference, url=refurl, secondary=True)
            add_quantity(events, name, 'alias', name, sec_source)
            if bibcode:
                source = add_source(events, name, bibcode=bibcode)
                sources = uniq_cdl([source, sec_source])
            else:
                sources = sec_source
            add_quantity(events, name, 'comovingdist', dist, sources)
        oldhostname = hostname

    events = journal_events(tasks, args, events)
    return events


def delete_old_event_files(*args):
    # Delete all old event JSON files
    repo_files = repo_file_list()
    for rfil in pbar(repo_files, desc='Deleting old events'):
        os.remove(rfil)
    return


def load_task_list():
    """Load the list of tasks in the `FILENAME.TASK_LIST` json file.

    A `TASK` object is created for each entry, with the parameters filled in.
    These are placed in an OrderedDict, sorted by the `priority` parameter, with positive values
    and then negative values, e.g. [0, 2, 10, -10, -1].
    """
    def_task_list_filename = FILENAME.TASK_LIST
    data = json.load(open(def_task_list_filename, 'r'))

    tasks = {}
    # `defaults` is a dictionary where each `key` is a task name, and values are its properties
    for key, val in data.items():
        tasks[key] = TASK(name=key, **val)

    # Sort entries as positive values, then negative values
    #    [0, 1, 2, 2, 10, -100, -10, -1]
    #    Tuples are sorted by first element (here: '0' if positive), then second (here normal order)
    tasks = OrderedDict(sorted(tasks.items(), key=lambda t: (t[1].priority < 0, t[1].priority)))
    print("load_task_list")
    for key, val in tasks.items():
        print(key, val.nice_name)

    print("\n\n")
    return tasks

if __name__ == '__main__':
    import_main()
