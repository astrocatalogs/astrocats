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
        ('ucb',             {'nicename': '%pre UCB photometry',          'update': False, 'archived': True}),
        ('sdss',            {'nicename': '%pre SDSS photometry',         'update': False}),
        ('csp',             {'nicename': '%pre CSP photometry',          'update': False}),
        ('itep',            {'nicename': '%pre ITEP',                    'update': False}),
        ('asiago_photo',          {'nicename': '%pre Asiago metadata',         'update': False}),
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
        ('cspspectra',      {'nicename': '%pre CSP spectra',             'update': False}),
        ('ucbspectra',      {'nicename': '%pre UCB spectra',             'update': True,  'archived': True}),
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
        if do_task(tasks, args, task, 'ucb'):
            from mtasks.general_data import do_ucb
            events = do_ucb(events, args, tasks)

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
        if do_task(tasks, args, task, 'csp'):
            from mtasks.vizier import do_csp
            events = do_csp(events, args, tasks)

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

        # Import ITEP
        if do_task(tasks, args, task, 'itep'):
            itepbadsources = ['2004ApJ...602..571B']

            needsbib = []
            with open(os.path.join(PATH.REPO_EXTERNAL, 'itep-refs.txt'), 'r') as f:
                refrep = f.read().splitlines()
            refrepf = dict(list(zip(refrep[1::2], refrep[::2])))
            f = open(os.path.join(PATH.REPO_EXTERNAL, 'itep-lc-cat-28dec2015.txt'), 'r')
            tsvin = csv.reader(f, delimiter='|', skipinitialspace=True)
            curname = ''
            for r, row in enumerate(pbar(tsvin, current_task)):
                if r <= 1 or len(row) < 7:
                    continue
                name = 'SN' + row[0].strip()
                mjd = str(jd_to_mjd(Decimal(row[1].strip())))
                band = row[2].strip()
                magnitude = row[3].strip()
                e_magnitude = row[4].strip()
                reference = row[6].strip().strip(',')

                if curname != name:
                    curname = name
                    name = add_event(tasks, args, events, name)

                    secondaryreference = 'Sternberg Astronomical Institute Supernova Light Curve Catalogue'
                    secondaryrefurl = 'http://dau.itep.ru/sn/node/72'
                    secondarysource = add_source(events, name, refname=secondaryreference, url=secondaryrefurl, secondary=True)
                    add_quantity(events, name, 'alias', name, secondarysource)

                    year = re.findall(r'\d+', name)[0]
                    add_quantity(events, name, 'discoverdate', year, secondarysource)
                if reference in refrepf:
                    bibcode = unescape(refrepf[reference])
                    source = add_source(events, name, bibcode=bibcode)
                else:
                    needsbib.append(reference)
                    source = add_source(events, name, refname=reference) if reference else ''

                if bibcode not in itepbadsources:
                    add_photometry(events, name, time=mjd, band=band, magnitude=magnitude,
                                   e_magnitude=e_magnitude, source=secondarysource + ',' + source)
            f.close()

            # Write out references that could use a bibcode
            needsbib = list(OrderedDict.fromkeys(needsbib))
            with open('../itep-needsbib.txt', 'w') as f:
                f.writelines(['%s\n' % i for i in needsbib])
            events = journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'tns'):
            session = requests.Session()
            csvtxt = load_cached_url(
                args, 'https://wis-tns.weizmann.ac.il/search?&num_page=1&format=html&sort=desc&order=id&format=csv&page=0',
                os.path.join(PATH.REPO_EXTERNAL, 'TNS/index.csv'))
            if not csvtxt:
                continue
            maxid = csvtxt.splitlines()[1].split(',')[0].strip('"')
            maxpages = ceil(int(maxid)/1000.)

            for page in pbar(range(maxpages), current_task):
                fname = os.path.join(PATH.REPO_EXTERNAL, 'TNS/page-') + str(page).zfill(2) + '.csv'
                if archived_task(tasks, args, 'tns') and os.path.isfile(fname) and page < 7:
                    with open(fname, 'r') as f:
                        csvtxt = f.read()
                else:
                    with open(fname, 'w') as f:
                        session = requests.Session()
                        response = session.get('https://wis-tns.weizmann.ac.il/search?&num_page=1000&format=html&edit[type]=&edit[objname]=&edit[id]=&sort=asc&order=id&display[redshift]=1&display[hostname]=1&display[host_redshift]=1&display[source_group_name]=1&display[programs_name]=1&display[internal_name]=1&display[isTNS_AT]=1&display[public]=1&display[end_pop_period]=0&display[spectra_count]=1&display[discoverymag]=1&display[discmagfilter]=1&display[discoverydate]=1&display[discoverer]=1&display[sources]=1&display[bibcode]=1&format=csv&page=' + str(page))
                        csvtxt = response.text
                        f.write(csvtxt)

                tsvin = csv.reader(csvtxt.splitlines(), delimiter=',')
                for ri, row in enumerate(pbar(tsvin, current_task, leave=False)):
                    if ri == 0:
                        continue
                    if row[4] and 'SN' not in row[4]:
                        continue
                    name = row[1].replace(' ', '')
                    name = add_event(tasks, args, events, name)
                    source = add_source(events, name, refname='Transient Name Server', url='https://wis-tns.weizmann.ac.il')
                    add_quantity(events, name, 'alias', name, source)
                    if row[2] and row[2] != '00:00:00.00':
                        add_quantity(events, name, 'ra', row[2], source)
                    if row[3] and row[3] != '+00:00:00.00':
                        add_quantity(events, name, 'dec', row[3], source)
                    if row[4]:
                        add_quantity(events, name, 'claimedtype', row[4].replace('SN', '').strip(), source)
                    if row[5]:
                        add_quantity(events, name, 'redshift', row[5], source, kind='spectroscopic')
                    if row[6]:
                        add_quantity(events, name, 'host', row[6], source)
                    if row[7]:
                        add_quantity(events, name, 'redshift', row[7], source, kind='host')
                    if row[8]:
                        add_quantity(events, name, 'discoverer', row[8], source)
                    # Currently, all events listing all possible observers. TNS bug?
                    # if row[9]:
                    #    observers = row[9].split(',')
                    #    for observer in observers:
                    #        add_quantity(events, name, 'observer', observer.strip(), source)
                    if row[10]:
                        add_quantity(events, name, 'alias', row[10], source)
                    if row[8] and row[14] and row[15] and row[16]:
                        survey = row[8]
                        magnitude = row[14]
                        band = row[15].split('-')[0]
                        mjd = astrotime(row[16]).mjd
                        add_photometry(events, name, time=mjd, magnitude=magnitude, band=band, survey=survey, source=source)
                    if row[16]:
                        date = row[16].split()[0].replace('-', '/')
                        if date != '0000/00/00':
                            date = date.replace('/00', '')
                            time = row[16].split()[1]
                            if time != '00:00:00':
                                ts = time.split(':')
                                date += pretty_num(timedelta(hours=int(ts[0]), minutes=int(ts[1]), seconds=int(ts[2])).total_seconds()/(24*60*60), sig=6).lstrip('0')
                            add_quantity(events, name, 'discoverdate', date, source)
                    if args.update:
                        events = journal_events(tasks, args, events)
            events = journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'snhunt'):
            html = load_cached_url(args, 'http://nesssi.cacr.caltech.edu/catalina/current.html', os.path.join(PATH.REPO_EXTERNAL, 'SNhunt/current.html'))
            if not html:
                continue
            text = html.splitlines()
            findtable = False
            for ri, row in enumerate(text):
                if 'Supernova Discoveries' in row:
                    findtable = True
                if findtable and '<table' in row:
                    tstart = ri+1
                if findtable and '</table>' in row:
                    tend = ri-1
            tablestr = '<html><body><table>'
            for row in text[tstart:tend]:
                if row[:3] == 'tr>':
                    tablestr = tablestr + '<tr>' + row[3:]
                else:
                    tablestr = tablestr + row
            tablestr = tablestr + '</table></body></html>'
            bs = BeautifulSoup(tablestr, 'html5lib')
            trs = bs.find('table').findAll('tr')
            for tr in pbar(trs, current_task):
                cols = [str(x.text) for x in tr.findAll('td')]
                if not cols:
                    continue
                name = re.sub('<[^<]+?>', '', cols[4]).strip().replace(' ', '').replace('SNHunt', 'SNhunt')
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, refname='Supernova Hunt', url='http://nesssi.cacr.caltech.edu/catalina/current.html')
                add_quantity(events, name, 'alias', name, source)
                host = re.sub('<[^<]+?>', '', cols[1]).strip().replace('_', ' ')
                add_quantity(events, name, 'host', host, source)
                add_quantity(events, name, 'ra', cols[2], source, unit='floatdegrees')
                add_quantity(events, name, 'dec', cols[3], source, unit='floatdegrees')
                dd = cols[0]
                discoverdate = dd[:4] + '/' + dd[4:6] + '/' + dd[6:8]
                add_quantity(events, name, 'discoverdate', discoverdate, source)
                discoverers = cols[5].split('/')
                for discoverer in discoverers:
                    add_quantity(events, name, 'discoverer', 'CRTS', source)
                    add_quantity(events, name, 'discoverer', discoverer, source)
                if args.update:
                    events = journal_events(tasks, args, events)
            events = journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'nedd'):
            f = open(os.path.join(PATH.REPO_EXTERNAL, 'NED25.12.1-D-10.4.0-20151123.csv'), 'r')
            data = csv.reader(f, delimiter=',', quotechar='"')
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
                    secondarysource = add_source(events, name, refname=reference, url=refurl, secondary=True)
                    add_quantity(events, name, 'alias', name, secondarysource)
                    if bibcode:
                        source = add_source(events, name, bibcode=bibcode)
                        sources = uniq_cdl([source, secondarysource])
                    else:
                        sources = secondarysource
                    add_quantity(events, name, 'comovingdist', dist, sources)
                oldhostname = hostname
            events = journal_events(tasks, args, events)

        # Import CPCS
        if do_task(tasks, args, task, 'cpcs'):
            jsontxt = load_cached_url(
                args, 'http://gsaweb.ast.cam.ac.uk/followup/list_of_alerts?format=json&num=100000&published=1&observed_only=1&hashtag=JG_530ad9462a0b8785bfb385614bf178c6',
                os.path.join(PATH.REPO_EXTERNAL, 'CPCS/index.json'))
            if not jsontxt:
                continue
            alertindex = json.loads(jsontxt, object_pairs_hook=OrderedDict)
            ids = [x['id'] for x in alertindex]

            for i, ai in enumerate(pbar(ids, current_task)):
                name = alertindex[i]['ivorn'].split('/')[-1].strip()
                # Skip a few weird entries
                if name == 'ASASSNli':
                    continue
                # Just use a whitelist for now since naming seems inconsistent
                if True in [x in name.upper() for x in ['GAIA', 'OGLE', 'ASASSN', 'MASTER', 'OTJ', 'PS1', 'IPTF']]:
                    name = name.replace('Verif', '').replace('_', ' ')
                    if 'ASASSN' in name and name[6] != '-':
                        name = 'ASASSN-' + name[6:]
                    if 'MASTEROTJ' in name:
                        name = name.replace('MASTEROTJ', 'MASTER OT J')
                    if 'OTJ' in name:
                        name = name.replace('OTJ', 'MASTER OT J')
                    if name.upper().startswith('IPTF'):
                        name = 'iPTF' + name[4:]
                    # Only add events that are classified as SN.
                    if event_exists(events, name):
                        continue
                    name = add_event(tasks, args, events, name)
                else:
                    continue

                secondarysource = add_source(events, name, refname='Cambridge Photometric Calibration Server', url='http://gsaweb.ast.cam.ac.uk/followup/', secondary=True)
                add_quantity(events, name, 'alias', name, secondarysource)
                add_quantity(events, name, 'ra', str(alertindex[i]['ra']), secondarysource, unit='floatdegrees')
                add_quantity(events, name, 'dec', str(alertindex[i]['dec']), secondarysource, unit='floatdegrees')

                alerturl = 'http://gsaweb.ast.cam.ac.uk/followup/get_alert_lc_data?alert_id=' + str(ai)
                source = add_source(events, name, refname='CPCS Alert ' + str(ai), url=alerturl)
                fname = os.path.join(PATH.REPO_EXTERNAL, 'CPCS/alert-') + str(ai).zfill(2) + '.json'
                if archived_task(tasks, args, 'cpcs') and os.path.isfile(fname):
                    with open(fname, 'r') as f:
                        jsonstr = f.read()
                else:
                    session = requests.Session()
                    response = session.get(alerturl + '&hashtag=JG_530ad9462a0b8785bfb385614bf178c6')
                    with open(fname, 'w') as f:
                        jsonstr = response.text
                        f.write(jsonstr)

                try:
                    cpcsalert = json.loads(jsonstr)
                except:
                    continue

                mjds = [round_sig(x, sig=9) for x in cpcsalert['mjd']]
                mags = [round_sig(x, sig=6) for x in cpcsalert['mag']]
                errs = [round_sig(x, sig=6) if (is_number(x) and float(x) > 0.0) else '' for x in cpcsalert['magerr']]
                bnds = cpcsalert['filter']
                obs  = cpcsalert['observatory']
                for mi, mjd in enumerate(mjds):
                    add_photometry(
                        events, name, time=mjd, magnitude=mags[mi], e_magnitude=errs[mi],
                        band=bnds[mi], observatory=obs[mi], source=uniq_cdl([source, secondarysource]))
                if args.update:
                    events = journal_events(tasks, args, events)
            events = journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'ptf'):
            # response = urllib.request.urlopen('http://wiserep.weizmann.ac.il/objects/list')
            # bs = BeautifulSoup(response, 'html5lib')
            # select = bs.find('select', {'name': 'objid'})
            # options = select.findAll('option')
            # for option in options:
            #    print(option.text)
            #    name = option.text
            #    if ((name.startswith('PTF') and is_number(name[3:5])) or
            #        name.startswith('PTFS') or name.startswith('iPTF')):
            #        name = add_event(tasks, args, events, name)

            if archived_task(tasks, args, 'ptf'):
                with open(os.path.join(PATH.REPO_EXTERNAL, 'PTF/update.html'), 'r') as f:
                    html = f.read()
            else:
                session = requests.Session()
                response = session.get('http://wiserep.weizmann.ac.il/spectra/update')
                html = response.text
                with open(os.path.join(PATH.REPO_EXTERNAL, 'PTF/update.html'), 'w') as f:
                    f.write(html)

            bs = BeautifulSoup(html, 'html5lib')
            select = bs.find('select', {'name': 'objid'})
            options = select.findAll('option')
            for option in options:
                name = option.text
                if (((name.startswith('PTF') and is_number(name[3:5])) or
                     name.startswith('PTFS') or name.startswith('iPTF'))):
                    if '(' in name:
                        alias = name.split('(')[0].strip(' ')
                        name = name.split('(')[-1].strip(') ').replace('sn', 'SN')
                        name = add_event(tasks, args, events, name)
                        source = add_source(events, name, bibcode='2012PASP..124..668Y')
                        add_quantity(events, name, 'alias', alias, source)
                    else:
                        name = add_event(tasks, args, events, name)

            with open(os.path.join(PATH.REPO_EXTERNAL, 'PTF/old-ptf-events.csv')) as f:
                for suffix in f.read().splitlines():
                    name = add_event(tasks, args, events, 'PTF' + suffix)
            with open(os.path.join(PATH.REPO_EXTERNAL, 'PTF/perly-2016.csv')) as f:
                for row in f.read().splitlines():
                    cols = [x.strip() for x in row.split(',')]
                    alias = ''
                    if cols[8]:
                        name = cols[8]
                        alias = 'PTF' + cols[0]
                    else:
                        name = 'PTF' + cols[0]
                    name = add_event(tasks, args, events, name)
                    source = add_source(events, name, bibcode='2016arXiv160408207P')
                    add_quantity(events, name, 'alias', name, source)
                    if alias:
                        add_quantity(events, name, 'alias', alias, source)
                    add_quantity(events, name, 'ra', cols[1], source)
                    add_quantity(events, name, 'dec', cols[2], source)
                    add_quantity(events, name, 'claimedtype', 'SLSN-' + cols[3], source)
                    add_quantity(events, name, 'redshift', cols[4], source, kind='spectroscopic')
                    maxdate = cols[6].replace('-', '/')
                    add_quantity(events, name, 'maxdate', maxdate.lstrip('<'), source, upperlimit=maxdate.startswith('<'))
                    add_quantity(events, name, 'ebv', cols[7], source, kind='spectroscopic')
                    name = add_event(tasks, args, events, 'PTF' + suffix)
            events = journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'des'):
            html = load_cached_url(args, 'https://portal.nersc.gov/des-sn/transients/', os.path.join(PATH.REPO_EXTERNAL, 'DES/transients.html'))
            if not html:
                continue
            bs = BeautifulSoup(html, 'html5lib')
            trs = bs.find('tbody').findAll('tr')
            for tri, tr in enumerate(pbar(trs, current_task)):
                name = ''
                source = ''
                if tri == 0:
                    continue
                tds = tr.findAll('td')
                for tdi, td in enumerate(tds):
                    if tdi == 0:
                        name = add_event(tasks, args, events, td.text.strip())
                    if tdi == 1:
                        (ra, dec) = [x.strip() for x in td.text.split('\xa0')]
                    if tdi == 6:
                        atellink = td.find('a')
                        if atellink:
                            atellink = atellink['href']
                        else:
                            atellink = ''

                sources = [add_source(events, name, url='https://portal.nersc.gov/des-sn/', refname='DES Bright Transients',
                                      acknowledgment='http://www.noao.edu/noao/library/NOAO_Publications_Acknowledgments.html#DESdatause')]
                if atellink:
                    sources.append(add_source(events, name, refname='ATel ' + atellink.split('=')[-1], url=atellink))
                sources += [add_source(events, name, bibcode='2012ApJ...753..152B'),
                            add_source(events, name, bibcode='2015AJ....150..150F'),
                            add_source(events, name, bibcode='2015AJ....150...82G'),
                            add_source(events, name, bibcode='2015AJ....150..172K')]
                sources = ','.join(sources)
                add_quantity(events, name, 'alias', name, sources)
                add_quantity(events, name, 'ra', ra, sources)
                add_quantity(events, name, 'dec', dec, sources)

                html2 = load_cached_url(args, 'https://portal.nersc.gov/des-sn/transients/' + name, os.path.join(PATH.REPO_EXTERNAL, 'DES/') + name + '.html')
                if not html2:
                    continue
                lines = html2.splitlines()
                for line in lines:
                    if 'var data = ' in line:
                        jsontxt = json.loads(line.split('=')[-1].rstrip(';'))
                        for i, band in enumerate(jsontxt['band']):
                            add_photometry(
                                events, name, time=jsontxt['mjd'][i], magnitude=jsontxt['mag'][i], e_magnitude=jsontxt['mag_error'][i],
                                band=band, observatory='CTIO', telescope='Blanco 4m', instrument='DECam',
                                upperlimit=True if float(jsontxt['snr'][i]) <= 3.0 else '', source=sources)
            events = journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'asassn'):
            html = load_cached_url(args, 'http://www.astronomy.ohio-state.edu/~assassin/sn_list.html', os.path.join(PATH.REPO_EXTERNAL, 'ASASSN/sn_list.html'))
            if not html:
                continue
            bs = BeautifulSoup(html, 'html5lib')
            trs = bs.find('table').findAll('tr')
            for tri, tr in enumerate(pbar(trs, current_task)):
                name = ''
                source = ''
                ra = ''
                dec = ''
                redshift = ''
                hostoff = ''
                claimedtype = ''
                host = ''
                atellink = ''
                typelink = ''
                if tri == 0:
                    continue
                tds = tr.findAll('td')
                for tdi, td in enumerate(tds):
                    if tdi == 1:
                        name = add_event(tasks, args, events, td.text.strip())
                        atellink = td.find('a')
                        if atellink:
                            atellink = atellink['href']
                        else:
                            atellink = ''
                    if tdi == 2:
                        discdate = td.text.replace('-', '/')
                    if tdi == 3:
                        ra = td.text
                    if tdi == 4:
                        dec = td.text
                    if tdi == 5:
                        redshift = td.text
                    if tdi == 8:
                        hostoff = td.text
                    if tdi == 9:
                        claimedtype = td.text
                        typelink = td.find('a')
                        if typelink:
                            typelink = typelink['href']
                        else:
                            typelink = ''
                    if tdi == 12:
                        host = td.text

                sources = [add_source(events, name, url='http://www.astronomy.ohio-state.edu/~assassin/sn_list.html', refname='ASAS-SN Supernovae')]
                typesources = sources[:]
                if atellink:
                    sources.append(add_source(events, name, refname='ATel ' + atellink.split('=')[-1], url=atellink))
                if typelink:
                    typesources.append(add_source(events, name, refname='ATel ' + typelink.split('=')[-1], url=typelink))
                sources = ','.join(sources)
                typesources = ','.join(typesources)
                add_quantity(events, name, 'alias', name, sources)
                add_quantity(events, name, 'discoverdate', discdate, sources)
                add_quantity(events, name, 'ra', ra, sources, unit='floatdegrees')
                add_quantity(events, name, 'dec', dec, sources, unit='floatdegrees')
                add_quantity(events, name, 'redshift', redshift, sources)
                add_quantity(events, name, 'hostoffset', hostoff, sources, unit='arcseconds')
                for ct in claimedtype.split('/'):
                    if ct != 'Unk':
                        add_quantity(events, name, 'claimedtype', ct, typesources)
                if host != 'Uncatalogued':
                    add_quantity(events, name, 'host', host, sources)
            events = journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'cspspectra'):
            oldname = ''
            file_names = glob(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'CSP/*'))
            for fi, fname in enumerate(pbar_strings(file_names, current_task=current_task)):
                filename = os.path.basename(fname)
                sfile = filename.split('.')
                if sfile[1] == 'txt':
                    continue
                sfile = sfile[0]
                fileparts = sfile.split('_')
                name = 'SN20' + fileparts[0][2:]
                name = get_preferred_name(events, name)
                if oldname and name != oldname:
                    events = journal_events(tasks, args, events)
                oldname = name
                name = add_event(tasks, args, events, name)
                telescope = fileparts[-2]
                instrument = fileparts[-1]
                source = add_source(events, name, bibcode='2013ApJ...773...53F')
                add_quantity(events, name, 'alias', name, source)

                f = open(fname, 'r')
                data = csv.reader(f, delimiter=' ', skipinitialspace=True)
                specdata = []
                for r, row in enumerate(data):
                    if row[0] == '#JDate_of_observation:':
                        jd = row[1].strip()
                        time = str(jd_to_mjd(Decimal(jd)))
                    elif row[0] == '#Redshift:':
                        add_quantity(events, name, 'redshift', row[1].strip(), source)
                    if r < 7:
                        continue
                    specdata.append(list(filter(None, [x.strip(' ') for x in row])))
                specdata = [list(i) for i in zip(*specdata)]
                wavelengths = specdata[0]
                fluxes = specdata[1]

                add_spectrum(
                    name=name, u_time='MJD', time=time, waveunit='Angstrom', fluxunit='erg/s/cm^2/Angstrom', wavelengths=wavelengths,
                    fluxes=fluxes, telescope=telescope, instrument=instrument, source=source, deredshifted=True, filename=filename)
                if args.travis and fi >= TRAVIS_QUERY_LIMIT:
                    break
            events = journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'ucbspectra'):
            secondaryreference = 'UCB Filippenko Group\'s Supernova Database (SNDB)'
            secondaryrefurl = 'http://heracles.astro.berkeley.edu/sndb/info'
            secondaryrefbib = '2012MNRAS.425.1789S'
            ucbspectracnt = 0

            jsontxt = load_cached_url(
                args, 'http://heracles.astro.berkeley.edu/sndb/download?id=allpubspec',
                os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'UCB/allpub.json'))
            if not jsontxt:
                continue

            spectra = json.loads(jsontxt)
            spectra = sorted(spectra, key=lambda k: k['ObjName'])
            oldname = ''
            for spectrum in pbar(spectra, desc=current_task):
                name = spectrum['ObjName']
                if oldname and name != oldname:
                    events = journal_events(tasks, args, events)
                oldname = name
                name = add_event(tasks, args, events, name)

                secondarysource = add_source(events, name, refname=secondaryreference, url=secondaryrefurl, bibcode=secondaryrefbib, secondary=True)
                add_quantity(events, name, 'alias', name, secondarysource)
                sources = [secondarysource]
                if spectrum['Reference']:
                    sources += [add_source(events, name, bibcode=spectrum['Reference'])]
                sources = uniq_cdl(sources)

                if spectrum['Type'] and spectrum['Type'].strip() != 'NoMatch':
                    for ct in spectrum['Type'].strip().split(','):
                        add_quantity(events, name, 'claimedtype', ct.replace('-norm', '').strip(), sources)
                if spectrum['DiscDate']:
                    add_quantity(events, name, 'discoverdate', spectrum['DiscDate'].replace('-', '/'), sources)
                if spectrum['HostName']:
                    add_quantity(events, name, 'host', urllib.parse.unquote(spectrum['HostName']).replace('*', ''), sources)
                if spectrum['UT_Date']:
                    epoch = str(spectrum['UT_Date'])
                    year = epoch[:4]
                    month = epoch[4:6]
                    day = epoch[6:]
                    sig = get_sig_digits(day) + 5
                    mjd = pretty_num(astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day)), sig=sig)
                filename = spectrum['Filename'] if spectrum['Filename'] else ''
                instrument = spectrum['Instrument'] if spectrum['Instrument'] else ''
                reducer = spectrum['Reducer'] if spectrum['Reducer'] else ''
                observer = spectrum['Observer'] if spectrum['Observer'] else ''
                snr = str(spectrum['SNR']) if spectrum['SNR'] else ''

                if not filename:
                    raise(ValueError('Filename not found for SNDB spectrum!'))
                if not spectrum['SpecID']:
                    raise(ValueError('ID not found for SNDB spectrum!'))

                filepath = os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'UCB/') + filename
                if archived_task('ucbspectra', tasks) and os.path.isfile(filepath):
                    with open(filepath, 'r') as f:
                        spectxt = f.read()
                else:
                    session = requests.Session()
                    response = session.get('http://heracles.astro.berkeley.edu/sndb/download?id=ds:' + str(spectrum['SpecID']))
                    spectxt = response.text
                    with open(filepath, 'w') as f:
                        f.write(spectxt)

                specdata = list(csv.reader(spectxt.splitlines(), delimiter=' ', skipinitialspace=True))
                startrow = 0
                for row in specdata:
                    if row[0][0] == '#':
                        startrow += 1
                    else:
                        break
                specdata = specdata[startrow:]

                haserrors = len(specdata[0]) == 3 and specdata[0][2] and specdata[0][2] != 'NaN'
                specdata = [list(i) for i in zip(*specdata)]

                wavelengths = specdata[0]
                fluxes = specdata[1]
                errors = ''
                if haserrors:
                    errors = specdata[2]

                if not list(filter(None, errors)):
                    errors = ''

                add_spectrum(
                    name=name, u_time='MJD', time=mjd, waveunit='Angstrom', fluxunit='Uncalibrated',
                    wavelengths=wavelengths, filename=filename, fluxes=fluxes, errors=errors, errorunit='Uncalibrated',
                    instrument=instrument, source=sources, snr=snr, observer=observer, reducer=reducer,
                    deredshifted=('-noz' in filename))
                ucbspectracnt = ucbspectracnt + 1
                if args.travis and ucbspectracnt >= TRAVIS_QUERY_LIMIT:
                    break
            events = journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'snfspectra'):
            eventfolders = next(os.walk(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'SNFactory')))[1]
            bibcodes = {'SN2005gj': '2006ApJ...650..510A', 'SN2006D': '2007ApJ...654L..53T', 'SN2007if': '2010ApJ...713.1073S', 'SN2011fe': '2013A&A...554A..27P'}
            oldname = ''
            snfcnt = 0
            for eventfolder in eventfolders:
                name = eventfolder
                name = get_preferred_name(events, name)
                if oldname and name != oldname:
                    events = journal_events(tasks, args, events)
                oldname = name
                name = add_event(tasks, args, events, name)
                secondaryreference = 'Nearby Supernova Factory'
                secondaryrefurl = 'http://snfactory.lbl.gov/'
                secondarybibcode = '2002SPIE.4836...61A'
                secondarysource = add_source(events, name, refname=secondaryreference, url=secondaryrefurl, bibcode=secondarybibcode, secondary=True)
                add_quantity(events, name, 'alias', name, secondarysource)
                bibcode = bibcodes[name]
                source = add_source(events, name, bibcode=bibcode)
                sources = uniq_cdl([source, secondarysource])
                eventspectra = glob(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'SNFactory/')+eventfolder+'/*.dat')
                for spectrum in eventspectra:
                    filename = os.path.basename(spectrum)
                    with open(spectrum) as f:
                        specdata = list(csv.reader(f, delimiter=' ', skipinitialspace=True))
                    specdata = list(filter(None, specdata))
                    newspec = []
                    time = ''
                    telescope = ''
                    instrument = ''
                    observer = ''
                    observatory = ''
                    if 'Keck_20060202_R' in spectrum:
                        time = '53768.23469'
                    elif 'Spectrum05_276' in spectrum:
                        time = pretty_num(astrotime('2005-10-03').mjd, sig=5)
                    elif 'Spectrum05_329' in spectrum:
                        time = pretty_num(astrotime('2005-11-25').mjd, sig=5)
                    elif 'Spectrum05_336' in spectrum:
                        time = pretty_num(astrotime('2005-12-02').mjd, sig=5)
                    for row in specdata:
                        if row[0][0] == '#':
                            joinrow = (' '.join(row)).split('=')
                            if len(joinrow) < 2:
                                continue
                            field = joinrow[0].strip('# ')
                            value = joinrow[1].split('/')[0].strip('\' ')
                            if not time:
                                if field == 'JD':
                                    time = str(jd_to_mjd(Decimal(value)))
                                elif field == 'MJD':
                                    time = value
                                elif field == 'MJD-OBS':
                                    time = value
                            if field == 'OBSERVER':
                                observer = value.capitalize()
                            if field == 'OBSERVAT':
                                observatory = value.capitalize()
                            if field == 'TELESCOP':
                                telescope = value.capitalize()
                            if field == 'INSTRUME':
                                instrument = value.capitalize()
                        else:
                            newspec.append(row)
                    if not time:
                        raise(ValueError('Time missing from spectrum.'))
                    specdata = newspec
                    haserrors = len(specdata[0]) == 3 and specdata[0][2] and specdata[0][2] != 'NaN'
                    specdata = [list(i) for i in zip(*specdata)]

                    wavelengths = specdata[0]
                    fluxes = specdata[1]
                    errors = ''
                    if haserrors:
                        errors = specdata[2]

                    add_spectrum(
                        name=name, u_time='MJD', time=time, waveunit='Angstrom', fluxunit='erg/s/cm^2/Angstrom',
                        wavelengths=wavelengths, fluxes=fluxes, errors=errors, observer=observer, observatory=observatory,
                        telescope=telescope, instrument=instrument,
                        errorunit=('Variance' if name == 'SN2011fe' else 'erg/s/cm^2/Angstrom'), source=sources, filename=filename)
                    snfcnt = snfcnt + 1
                    if args.travis and snfcnt % TRAVIS_QUERY_LIMIT == 0:
                        break
            events = journal_events(tasks, args, events)

        if do_task(tasks, args, task, 'superfitspectra'):
            sfdirs = glob(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'superfit/*'))
            for sfdir in pbar(sfdirs, desc=current_task):
                sffiles = sorted(glob(sfdir + '/*.dat'))
                lastname = ''
                oldname = ''
                for sffile in pbar(sffiles, desc=current_task):
                    basename = os.path.basename(sffile)
                    name = basename.split('.')[0]
                    if name.startswith('sn'):
                        name = 'SN' + name[2:]
                        if len(name) == 7:
                            name = name[:6] + name[6].upper()
                    elif name.startswith('ptf'):
                        name = 'PTF' + name[3:]

                    if 'theory' in name:
                        continue
                    if event_exists(events, name):
                        prefname = get_preferred_name(events, name)
                        if 'spectra' in events[prefname] and lastname != prefname:
                            continue
                    if oldname and name != oldname:
                        events = journal_events(tasks, args, events)
                    oldname = name
                    name = add_event(tasks, args, events, name)
                    epoch = basename.split('.')[1]
                    (mldt, mlmag, mlband, mlsource) = get_max_light(events, name)
                    if mldt:
                        epoff = Decimal(0.0) if epoch == 'max' else (Decimal(epoch[1:]) if epoch[0] == 'p' else -Decimal(epoch[1:]))
                    else:
                        epoff = ''

                    source = add_source(events, name, refname='Superfit', url='http://www.dahowell.com/superfit.html', secondary=True)
                    add_quantity(events, name, 'alias', name, source)

                    with open(sffile) as f:
                        rows = f.read().splitlines()
                    specdata = []
                    for row in rows:
                        if row.strip():
                            specdata.append(list(filter(None, re.split('\t+|\s+', row, maxsplit=0))))
                    specdata = [[x.replace('D', 'E') for x in list(i)] for i in zip(*specdata)]
                    wavelengths = specdata[0]
                    fluxes = specdata[1]

                    mlmjd = str(Decimal(astrotime('-'.join([str(mldt.year), str(mldt.month), str(mldt.day)])).mjd) + epoff) if (epoff != '') else ''
                    add_spectrum(
                        name, u_time='MJD' if mlmjd else '', time=mlmjd, waveunit='Angstrom', fluxunit='Uncalibrated',
                        wavelengths=wavelengths, fluxes=fluxes, source=source)

                    lastname = name
                events = journal_events(tasks, args, events)

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

    sys.exit(0)


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
