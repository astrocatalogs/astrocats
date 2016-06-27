"""General data import tasks.
"""
import csv
import os
import re
from collections import OrderedDict
from datetime import datetime
from glob import glob
from html import unescape
from math import ceil, log10

import requests
from astropy.time import Time as astrotime
from bs4 import BeautifulSoup

from cdecimal import Decimal
from scripts import PATH
from scripts.utils import (is_number, pbar, pbar_strings, pretty_num,
                           round_sig, single_spaces)

from .. import Events
from ..constants import TRAVIS_QUERY_LIMIT
from ..Events import load_event_from_file
from ..funcs import (add_photometry, add_spectrum, event_exists, jd_to_mjd,
                     load_cached_url, make_date_string, uniq_cdl)


def do_external_radio(events, stubs, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    path_pattern = os.path.join(PATH.REPO_EXTERNAL_RADIO, '*.txt')
    for datafile in pbar_strings(glob(path_pattern), desc=current_task):
        oldname = os.path.basename(datafile).split('.')[0]
        events, name = Events.add_event(tasks, args, events, oldname, log)
        radiosourcedict = OrderedDict()
        with open(datafile, 'r') as ff:
            for li, line in enumerate([xx.strip() for xx in ff.read().splitlines()]):
                if line.startswith('(') and li <= len(radiosourcedict):
                    key = line.split()[0]
                    bibc = line.split()[-1]
                    radiosourcedict[key] = events[
                        name].add_source(bibcode=bibc)
                elif li in [xx + len(radiosourcedict) for xx in range(3)]:
                    continue
                else:
                    cols = list(filter(None, line.split()))
                    source = radiosourcedict[cols[6]]
                    add_photometry(
                        events, name, time=cols[0], frequency=cols[
                            2], u_frequency='GHz',
                        fluxdensity=cols[3], e_fluxdensity=cols[
                            4], u_fluxdensity='µJy',
                        instrument=cols[5], source=source)
                    events[name].add_quantity('alias', oldname, source)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_external_xray(events, stubs, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    path_pattern = os.path.join(PATH.REPO_EXTERNAL_XRAY, '*.txt')
    for datafile in pbar_strings(glob(path_pattern), desc=current_task):
        oldname = os.path.basename(datafile).split('.')[0]
        events, name = Events.add_event(tasks, args, events, oldname, log)
        with open(datafile, 'r') as ff:
            for li, line in enumerate(ff.read().splitlines()):
                if li == 0:
                    source = events[name].add_source(bibcode=line.split()[-1])
                elif li in [1, 2, 3]:
                    continue
                else:
                    cols = list(filter(None, line.split()))
                    add_photometry(
                        events, name, time=cols[:2],
                        energy=cols[2:4], u_energy='keV', counts=cols[4], flux=cols[6],
                        unabsorbedflux=cols[8], u_flux='ergs/ss/cm^2',
                        photonindex=cols[15], instrument=cols[
                            17], nhmw=cols[11],
                        upperlimit=(float(cols[5]) < 0), source=source)
                    events[name].add_quantity('alias', oldname, source)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_internal(events, stubs, args, tasks, task_obj, log):
    """Load events from files in the 'internal' repository, and save them.
    """
    current_task = task_obj.current_task(args)
    path_pattern = os.path.join(PATH.REPO_INTERNAL, '*.json')
    files = glob(path_pattern)
    log.debug("found {} files matching '{}'".format(len(files), path_pattern))
    log.error("`do_internal` 'update' section is disabled")
    for datafile in pbar_strings(files, desc=current_task):
        # FIX: do we still need this difference?
        '''
        if args.update:
            if not load_event_from_file(events, args, tasks, path=datafile,
                                        clean=True, delete=False, append=True):
                raise IOError('Failed to find specified file.')
        else:
            if not load_event_from_file(events, args, tasks, path=datafile,
                                        clean=True, delete=False):
                raise IOError('Failed to find specified file.')
        '''
        new_event = load_event_from_file(events, args, tasks, log, path=datafile,
                                         clean=True, delete=False)
        events.update({new_event.name: new_event})

    return events


def do_simbad(events, stubs, args, tasks, task_obj, log):
    # Simbad.list_votable_fields()
    # Some coordinates that SIMBAD claims belong to the SNe actually belong to
    # the host.
    simbadmirrors = ['http://simbad.harvard.edu/simbad/sim-script',
                     'http://simbad.u-strasbg.fr/simbad/sim-script']
    simbadbadcoordbib = ['2013ApJ...770..107C']
    simbadbadnamebib = ['2004AJ....127.2809W', '2005MNRAS.364.1419Z',
                        '2015A&A...574A.112D', '2011MNRAS.417..916G', '2002ApJ...566..880G']
    simbadbannedcats = ['[TBV2008]', 'OGLE-MBR']
    customSimbad = Simbad()
    customSimbad.ROW_LIMIT = -1
    customSimbad.TIMEOUT = 120
    customSimbad.add_votable_fields('otype', 'sptype', 'sp_bibcode', 'id')
    for mirror in simbadmirrors:
        customSimbad.SIMBAD_URL = mirror
        try:
            table = customSimbad.query_criteria('maintype=SN | maintype="SN?"')
        except:
            continue
        else:
            break

    # 2000A&AS..143....9W
    for brow in pbar(table, current_task=current_task):
        row = {x: re.sub(r'b\'(.*)\'', r'\1',
                         str(brow[x])) for x in brow.colnames}
        # Skip items with no bibliographic info aside from SIMBAD, too
        # error-prone
        if row['OTYPE'] == 'Candidate_SN*' and not row['SP_TYPE']:
            continue
        if not row['COO_BIBCODE'] and not row['SP_BIBCODE'] and not row['SP_BIBCODE_2']:
            continue
        if any([x in row['MAIN_ID'] for x in simbadbannedcats]):
            continue
        if row['COO_BIBCODE'] and row['COO_BIBCODE'] in simbadbadnamebib:
            continue
        name = single_spaces(re.sub(r'\[[^)]*\]', '', row['MAIN_ID']).strip())
        if name == 'SN':
            continue
        if is_number(name):
            continue
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(srcname='SIMBAD astronomical database', bibcode="2000A&AS..143....9W",
                                         url="http://simbad.u-strasbg.fr/", secondary=True)
        aliases = row['ID'].split(',')
        for alias in aliases:
            if any([x in alias for x in simbadbannedcats]):
                continue
            ali = single_spaces(re.sub(r'\[[^)]*\]', '', alias).strip())
            if is_number(ali):
                continue
            ali = name_clean(ali)
            events[name].add_quantity('alias', ali, source)
        if row['COO_BIBCODE'] and row['COO_BIBCODE'] not in simbadbadcoordbib:
            csources = ','.join(
                [source, events[name].add_source(bibcode=row['COO_BIBCODE'])])
            events[name].add_quantity('ra', row['RA'], csources)
            events[name].add_quantity('dec', row['DEC'], csources)
        if row['SP_BIBCODE']:
            ssources = uniq_cdl([source, events[name].add_source(bibcode=row['SP_BIBCODE'])] +
                                ([events[name].add_source(bibcode=row['SP_BIBCODE_2'])] if row['SP_BIBCODE_2'] else []))
            events[name].add_quantity('claimedtype', row['SP_TYPE'].replace(
                'SN.', '').replace('SN', '').replace('(~)', '').strip(': '), ssources)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events
