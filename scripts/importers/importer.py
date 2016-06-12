#!/usr/local/bin/python3.5

import argparse
from astropy.time import Time as astrotime
from bs4 import BeautifulSoup, Tag, NavigableString
import calendar
from cdecimal import Decimal
import codecs
from collections import OrderedDict
import csv
from datetime import timedelta
from glob import glob
from html import unescape
import json
from math import log10, floor, ceil
import numpy as np
import os
import re
import requests
import resource
from string import ascii_letters
import sys
import urllib
import warnings

from .. utils import get_sig_digits, pretty_num, pbar, pbar_strings, is_number, round_sig, tprint, \
    repo_file_list
from . funcs import add_event, add_photometry, add_quantity, add_source, add_spectrum, \
    archived_task, clean_snname, derive_and_sanitize, delete_old_event_files, \
    do_task, event_exists, get_aliases, get_bibauthor_dict, get_extinctions_dict, \
    get_max_light, get_preferred_name, has_task, jd_to_mjd, journal_events, \
    load_cached_url, make_date_string, merge_duplicates, \
    set_preferred_names, uniq_cdl, utf8, write_all_events
from scripts import PATH, FILENAME
from . constants import TRAVIS_QUERY_LIMIT, TASK


def import_main():
    """
    """
    args = load_args()
    current_task = ''
    # eventnames = []
    events = OrderedDict()
    warnings.filterwarnings('ignore', r'Warning: converting a masked element to nan.')

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
        ('asiago_spectra',   {'nicename': '%pre Asiago spectra',          'update': True}),
        ('wiserepspectra',  {'nicename': '%pre WISeREP spectra',         'update': False}),
        ('cfa_spectra',     {'nicename': '%pre CfA archive spectra',     'update': False}),
        ('snlsspectra',     {'nicename': '%pre SNLS spectra',            'update': False}),
        ('csp_spectra',      {'nicename': '%pre CSP spectra',             'update': False}),
        ('ucb_spectra',      {'nicename': '%pre UCB spectra',             'update': True,  'archived': True}),
        ('suspect_spectra', {'nicename': '%pre SUSPECT spectra',         'update': False}),
        ('snfspectra',      {'nicename': '%pre SNH spectra',             'update': False}),
        ('superfitspectra', {'nicename': '%pre Superfit spectra',        'update': False}),
        ('mergeduplicates', {'nicename': 'Merging duplicates',           'update': False}),
        ('setprefnames',    {'nicename': 'Setting preferred names',      'update': False}),
        ('writeevents',     {'nicename': 'Writing events',               'update': True})
    ])

    for task in tasks:
        if do_task(tasks, args, task, 'deleteoldevents'):
            # Delete `current_task` here and wrap deletion in `pbar_strings` ??
            current_task = 'Deleting old events'
            delete_old_event_files()

        # Import data provided directly to OSC, in standard format
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
            from mtasks.suspect import do_suspect_photometry
            events = do_suspect_photometry(events, args, tasks)

        if do_task(tasks, args, task, 'suspect_spectra'):
            from mtasks.suspect import do_suspect_spectra
            events = do_suspect_spectra(events, args, tasks)

        # CfA data
        if do_task(tasks, args, task, 'cfa_photo'):
            from mtasks.cfa import do_cfa_photometry
            events = do_cfa_photometry(events, args, tasks)

        if do_task(tasks, args, task, 'cfa_spectra'):
            from mtasks.cfa import do_cfa_spectra
            events = do_cfa_spectra(events, args, tasks)

        if do_task(tasks, args, task, 'wiserepspectra'):
            from mtasks.wiserep import do_wiserep_spectra
            events = do_wiserep_spectra(events, args, tasks)

        # Import primary data sources from Vizier
        if do_task(tasks, args, task, 'vizier'):
            from mtasks.vizier import do_vizier
            events = do_vizier(events, args, tasks)

        if do_task(tasks, args, task, 'lennarz'):
            from mtasks.vizier import do_lennarz
            events = do_lennarz(events, args, tasks)

        if do_task(tasks, args, task, 'snlsspectra'):
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

        if do_task(tasks, args, task, 'snfspectra'):
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

        if do_task(tasks, args, task, 'superfitspectra'):
            from mtasks.general_data import do_superfit_spectra
            events = do_superfit_spectra(events, args, tasks)

        if do_task(tasks, args, task, 'mergeduplicates'):
            if args.update and not len(events):
                tprint('No sources changed, event files unchanged in update.')
                sys.exit(1)
            merge_duplicates(tasks, args, events)

        if do_task(tasks, args, task, 'setprefnames'):
            set_preferred_names(tasks, args, events)

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
    return tasks


def load_args():
    parser = argparse.ArgumentParser(description='Generate a catalog JSON file and plot HTML files from SNE data.')
    parser.add_argument('--update', '-u',       dest='update',      help='Only update catalog using live sources.',    default=False, action='store_true')
    parser.add_argument('--verbose', '-v',      dest='verbose',     help='Print more messages to the screen.',         default=False, action='store_true')
    parser.add_argument('--refresh', '-r',      dest='refresh',     help='Ignore most task caches.',                   default=False, action='store_true')
    parser.add_argument('--full-refresh', '-f', dest='fullrefresh', help='Ignore all task caches.',                    default=False, action='store_true')
    parser.add_argument('--archived', '-a',     dest='archived',    help='Always use task caches.',                    default=False, action='store_true')
    parser.add_argument('--travis', '-tr',      dest='travis',      help='Run import script in test mode for Travis.', default=False, action='store_true')
    parser.add_argument('--refreshlist', '-rl', dest='refreshlist', help='Comma-delimited list of caches to clear.',   default='')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    import_main()
