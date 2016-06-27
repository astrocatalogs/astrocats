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


def do_snhunt(events, stubs, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    snh_url = 'http://nesssi.cacr.caltech.edu/catalina/current.html'
    html = load_cached_url(args, current_task, snh_url, os.path.join(
        PATH.REPO_EXTERNAL, 'SNhunt/current.html'))
    if not html:
        return events
    text = html.splitlines()
    findtable = False
    for ri, row in enumerate(text):
        if 'Supernova Discoveries' in row:
            findtable = True
        if findtable and '<table' in row:
            tstart = ri + 1
        if findtable and '</table>' in row:
            tend = ri - 1
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
        cols = [str(xx.text) for xx in tr.findAll('td')]
        if not cols:
            continue
        name = re.sub('<[^<]+?>', '', cols[4]
                      ).strip().replace(' ', '').replace('SNHunt', 'SNhunt')
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(srcname='Supernova Hunt', url=snh_url)
        events[name].add_quantity('alias', name, source)
        host = re.sub('<[^<]+?>', '', cols[1]).strip().replace('_', ' ')
        events[name].add_quantity('host', host, source)
        events[name].add_quantity('ra', cols[2], source, unit='floatdegrees')
        events[name].add_quantity('dec', cols[3], source, unit='floatdegrees')
        dd = cols[0]
        discoverdate = dd[:4] + '/' + dd[4:6] + '/' + dd[6:8]
        events[name].add_quantity('discoverdate', discoverdate, source)
        discoverers = cols[5].split('/')
        for discoverer in discoverers:
            events[name].add_quantity('discoverer', 'CRTS', source)
            events[name].add_quantity('discoverer', discoverer, source)
        if args.update:
            events, stubs = Events.journal_events(
                tasks, args, events, stubs, log)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_superfit_spectra(events, stubs, args, tasks, task_obj, log):
    from .. funcs import get_max_light, get_preferred_name
    superfit_url = 'http://www.dahowell.com/superfit.html'
    current_task = task_obj.current_task(args)
    sfdirs = list(glob(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'superfit/*')))
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
                events, stubs = Events.journal_events(
                    tasks, args, events, stubs, log)
            oldname = name
            events, name = Events.add_event(tasks, args, events, name, log)
            epoch = basename.split('.')[1]
            (mldt, mlmag, mlband, mlsource) = get_max_light(events, name)
            if mldt:
                if epoch == 'max':
                    epoff = Decimal(0.0)
                elif epoch[0] == 'p':
                    epoff = Decimal(epoch[1:])
                else:
                    epoff = -Decimal(epoch[1:])
            else:
                epoff = ''

            source = events[name].add_source(
                srcname='Superfit', url=superfit_url, secondary=True)
            events[name].add_quantity('alias', oldname, source)

            with open(sffile) as ff:
                rows = ff.read().splitlines()
            specdata = []
            for row in rows:
                if row.strip():
                    specdata.append(
                        list(filter(None, re.split('\tt+|\ss+', row, maxsplit=0))))
            specdata = [[xx.replace('D', 'E') for xx in list(ii)]
                        for ii in zip(*specdata)]
            wavelengths = specdata[0]
            fluxes = specdata[1]

            if epoff != '':
                mlmjd = astrotime(
                    '-'.join([str(mldt.year), str(mldt.month), str(mldt.day)])).mjd
                mlmjd = str(Decimal(mlmjd) + epoff)
            else:
                mlmjd = ''
            add_spectrum(
                events, name, 'Angstrom', 'Uncalibrated', u_time='MJD' if mlmjd else '', time=mlmjd,
                wavelengths=wavelengths, fluxes=fluxes, source=source)

            lastname = name

        events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_tns(events, stubs, args, tasks, task_obj, log):
    from datetime import timedelta
    session = requests.Session()
    current_task = task_obj.current_task(args)
    tns_url = 'https://wis-tns.weizmann.ac.il/'
    search_url = tns_url + \
        'search?&num_page=1&format=html&sort=desc&order=id&format=csv&page=0'
    csvtxt = load_cached_url(args, current_task, search_url, os.path.join(
        PATH.REPO_EXTERNAL, 'TNS/index.csv'))
    if not csvtxt:
        return events
    maxid = csvtxt.splitlines()[1].split(',')[0].strip('"')
    maxpages = ceil(int(maxid) / 1000.)

    for page in pbar(range(maxpages), current_task):
        fname = os.path.join(PATH.REPO_EXTERNAL, 'TNS/page-') + \
            str(page).zfill(2) + '.csv'
        if task_obj.load_archive(args) and os.path.isfile(fname) and page < 7:
            with open(fname, 'r') as tns_file:
                csvtxt = tns_file.read()
        else:
            with open(fname, 'w') as tns_file:
                session = requests.Session()
                ses_url = (tns_url + 'search?&num_page=1000&format=html&edit'
                           '[type]=&edit[objname]=&edit[id]=&sort=asc&order=id&display[redshift]=1'
                           '&display[hostname]=1&display[host_redshift]=1'
                           '&display[source_group_name]=1&display[programs_name]=1'
                           '&display[internal_name]=1&display[isTNS_AT]=1&display[public]=1'
                           '&display[end_pop_period]=0&display[spectra_count]=1'
                           '&display[discoverymag]=1&display[discmagfilter]=1'
                           '&display[discoverydate]=1&display[discoverer]=1&display[sources]=1'
                           '&display[bibcode]=1&format=csv&page=' + str(page))
                response = session.get(ses_url)
                csvtxt = response.text
                tns_file.write(csvtxt)

        tsvin = list(csv.reader(csvtxt.splitlines(), delimiter=','))
        for ri, row in enumerate(pbar(tsvin, current_task, leave=False)):
            if ri == 0:
                continue
            if row[4] and 'SN' not in row[4]:
                continue
            name = row[1].replace(' ', '')
            events, name = Events.add_event(tasks, args, events, name, log)
            source = events[name].add_source(
                srcname='Transient Name Server', url=tns_url)
            events[name].add_quantity('alias', name, source)
            if row[2] and row[2] != '00:00:00.00':
                events[name].add_quantity('ra', row[2], source)
            if row[3] and row[3] != '+00:00:00.00':
                events[name].add_quantity('dec', row[3], source)
            if row[4]:
                events[name].add_quantity(
                    'claimedtype', row[4].replace('SN', '').strip(), source)
            if row[5]:
                events[name].add_quantity(
                    'redshift', row[5], source, kind='spectroscopic')
            if row[6]:
                events[name].add_quantity('host', row[6], source)
            if row[7]:
                events[name].add_quantity(
                    'redshift', row[7], source, kind='host')
            if row[8]:
                events[name].add_quantity('discoverer', row[8], source)
            # Currently, all events listing all possible observers. TNS bug?
            # if row[9]:
            #    observers = row[9].split(',')
            #    for observer in observers:
            #        events[name].add_quantity('observer', observer.strip(), source)
            if row[10]:
                events[name].add_quantity('alias', row[10], source)
            if row[8] and row[14] and row[15] and row[16]:
                survey = row[8]
                magnitude = row[14]
                band = row[15].split('-')[0]
                mjd = astrotime(row[16]).mjd
                add_photometry(events, name, time=mjd, magnitude=magnitude, band=band,
                               survey=survey, source=source)
            if row[16]:
                date = row[16].split()[0].replace('-', '/')
                if date != '0000/00/00':
                    date = date.replace('/00', '')
                    time = row[16].split()[1]
                    if time != '00:00:00':
                        ts = time.split(':')
                        dt = timedelta(hours=int(ts[0]), minutes=int(
                            ts[1]), seconds=int(ts[2]))
                        date += pretty_num(dt.total_seconds() /
                                           (24 * 60 * 60), sig=6).lstrip('0')
                    events[name].add_quantity('discoverdate', date, source)
            if args.update:
                events, stubs = Events.journal_events(
                    tasks, args, events, stubs, log)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
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
