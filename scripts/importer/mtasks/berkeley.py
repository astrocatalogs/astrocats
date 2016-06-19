"""Import data from UC Berkeley group.
"""
from astropy.time import Time as astrotime
import csv
import json
from math import floor
import os
import requests
import urllib

from scripts import PATH
from .. constants import TRAVIS_QUERY_LIMIT
from .. funcs import add_event, add_photometry, add_source, add_quantity, add_spectrum, \
    journal_events, load_cached_url, uniq_cdl
from ... utils import get_sig_digits, pbar, pretty_num


def do_ucb_photo(events, args, tasks, task_obj):
    current_task = task_obj.current_task(args)
    sec_ref = 'UCB Filippenko Group\'s Supernova Database (SNDB)'
    sec_refurl = 'http://heracles.astro.berkeley.edu/sndb/info'
    sec_refbib = '2012MNRAS.425.1789S'

    jsontxt = load_cached_url(
        args, 'http://heracles.astro.berkeley.edu/sndb/download?id=allpubphot',
        os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'UCB/allpub.json'))
    if not jsontxt:
        return events

    photom = json.loads(jsontxt)
    photom = sorted(photom, key=lambda kk: kk['ObjName'])
    for phot in pbar(photom, desc=current_task):
        name = phot['ObjName']
        name = add_event(tasks, args, events, name)

        sec_source = add_source(events, name, srcname=sec_ref, url=sec_refurl, bibcode=sec_refbib,
                                secondary=True)
        add_quantity(events, name, 'alias', name, sec_source)
        sources = [sec_source]
        if phot['Reference']:
            sources += [add_source(events, name, bibcode=phot['Reference'])]
        sources = uniq_cdl(sources)

        if phot['Type'] and phot['Type'].strip() != 'NoMatch':
            for ct in phot['Type'].strip().split(','):
                add_quantity(events, name, 'claimedtype', ct.replace('-norm', '').strip(), sources)
        if phot['DiscDate']:
            add_quantity(events, name, 'discoverdate', phot['DiscDate'].replace('-', '/'), sources)
        if phot['HostName']:
            host = urllib.parse.unquote(phot['HostName']).replace('*', '')
            add_quantity(events, name, 'host', host, sources)
        filename = phot['Filename'] if phot['Filename'] else ''

        if not filename:
            raise ValueError('Filename not found for SNDB phot!')
        if not phot['PhotID']:
            raise ValueError('ID not found for SNDB phot!')

        filepath = os.path.join(PATH.REPO_EXTERNAL, 'SNDB/') + filename
        if task_obj.load_archive(args) and os.path.isfile(filepath):
            with open(filepath, 'r') as ff:
                phottxt = ff.read()
        else:
            session = requests.Session()
            response = session.get('http://heracles.astro.berkeley.edu/sndb/download?id=dp:' +
                                   str(phot['PhotID']))
            phottxt = response.text
            with open(filepath, 'w') as ff:
                ff.write(phottxt)

        tsvin = csv.reader(phottxt.splitlines(), delimiter=' ', skipinitialspace=True)

        for rr, row in enumerate(tsvin):
            if len(row) > 0 and row[0] == "#":
                continue
            mjd = row[0]
            magnitude = row[1]
            if magnitude and float(magnitude) > 99.0:
                continue
            e_mag = row[2]
            band = row[4]
            telescope = row[5]
            add_photometry(
                events, name, time=mjd, telescope=telescope, band=band, magnitude=magnitude,
                e_magnitude=e_mag, source=sources)

    events = journal_events(tasks, args, events)
    return events


def do_ucb_spectra(events, args, tasks, task_obj):
    current_task = task_obj.current_task(args)
    sec_reference = 'UCB Filippenko Group\'s Supernova Database (SNDB)'
    sec_refurl = 'http://heracles.astro.berkeley.edu/sndb/info'
    sec_refbib = '2012MNRAS.425.1789S'
    ucbspectracnt = 0

    jsontxt = load_cached_url(
        args, 'http://heracles.astro.berkeley.edu/sndb/download?id=allpubspec',
        os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'UCB/allpub.json'))
    if not jsontxt:
        return events

    spectra = json.loads(jsontxt)
    spectra = sorted(spectra, key=lambda kk: kk['ObjName'])
    oldname = ''
    for spectrum in pbar(spectra, desc=current_task):
        name = spectrum['ObjName']
        if oldname and name != oldname:
            events = journal_events(tasks, args, events)
        oldname = name
        name = add_event(tasks, args, events, name)

        sec_source = add_source(
            events, name, refname=sec_reference, url=sec_refurl, bibcode=sec_refbib, secondary=True)
        add_quantity(events, name, 'alias', name, sec_source)
        sources = [sec_source]
        if spectrum['Reference']:
            sources += [add_source(events, name, bibcode=spectrum['Reference'])]
        sources = uniq_cdl(sources)

        if spectrum['Type'] and spectrum['Type'].strip() != 'NoMatch':
            for ct in spectrum['Type'].strip().split(','):
                add_quantity(events, name, 'claimedtype', ct.replace('-norm', '').strip(), sources)
        if spectrum['DiscDate']:
            ddate = spectrum['DiscDate'].replace('-', '/')
            add_quantity(events, name, 'discoverdate', ddate, sources)
        if spectrum['HostName']:
            host = urllib.parse.unquote(spectrum['HostName']).replace('*', '')
            add_quantity(events, name, 'host', host, sources)
        if spectrum['UT_Date']:
            epoch = str(spectrum['UT_Date'])
            year = epoch[:4]
            month = epoch[4:6]
            day = epoch[6:]
            sig = get_sig_digits(day) + 5
            mjd = astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd
            mjd = pretty_num(mjd + float(day) - floor(float(day)), sig=sig)
        filename = spectrum['Filename'] if spectrum['Filename'] else ''
        instrument = spectrum['Instrument'] if spectrum['Instrument'] else ''
        reducer = spectrum['Reducer'] if spectrum['Reducer'] else ''
        observer = spectrum['Observer'] if spectrum['Observer'] else ''
        snr = str(spectrum['SNR']) if spectrum['SNR'] else ''

        if not filename:
            raise ValueError('Filename not found for SNDB spectrum!')
        if not spectrum['SpecID']:
            raise ValueError('ID not found for SNDB spectrum!')

        filepath = os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'UCB/') + filename
        if task_obj.load_archive(args) and os.path.isfile(filepath):
            with open(filepath, 'r') as ff:
                spectxt = ff.read()
        else:
            session = requests.Session()
            response = session.get('http://heracles.astro.berkeley.edu/sndb/download?id=ds:' +
                                   str(spectrum['SpecID']))
            spectxt = response.text
            with open(filepath, 'w') as ff:
                ff.write(spectxt)

        specdata = list(csv.reader(spectxt.splitlines(), delimiter=' ', skipinitialspace=True))
        startrow = 0
        for row in specdata:
            if row[0][0] == '#':
                startrow += 1
            else:
                break
        specdata = specdata[startrow:]

        haserrors = len(specdata[0]) == 3 and specdata[0][2] and specdata[0][2] != 'NaN'
        specdata = [list(ii) for ii in zip(*specdata)]

        wavelengths = specdata[0]
        fluxes = specdata[1]
        errors = ''
        if haserrors:
            errors = specdata[2]

        if not list(filter(None, errors)):
            errors = ''

        units = 'Uncalibrated'
        add_spectrum(
            name=name, u_time='MJD', time=mjd, waveunit='Angstrom', fluxunit=units,
            wavelengths=wavelengths, filename=filename, fluxes=fluxes, errors=errors,
            errorunit=units, instrument=instrument, source=sources, snr=snr, observer=observer,
            reducer=reducer, deredshifted=('-noz' in filename))
        ucbspectracnt = ucbspectracnt + 1
        if args.travis and ucbspectracnt >= TRAVIS_QUERY_LIMIT:
            break

    events = journal_events(tasks, args, events)
    return events