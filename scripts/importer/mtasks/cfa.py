"""Imports from the 'CfA' sources.
"""
import csv
import os
from glob import glob
from math import floor

from astropy.time import Time as astrotime

from cdecimal import Decimal
from scripts import PATH

from .. import Events
from ...utils import is_number, pbar, pbar_strings
from ..constants import ACKN_CFA, TRAVIS_QUERY_LIMIT
from ..funcs import (add_photometry, add_spectrum, clean_snname,
                     get_preferred_name, jd_to_mjd, uniq_cdl)


def do_cfa_photo(catalog):
    from html import unescape
    import re
    current_task = task_obj.current_task(args)
    file_names = glob(os.path.join(PATH.REPO_EXTERNAL, 'cfa-input/*.dat'))
    for fname in pbar_strings(file_names, desc=current_task):
        f = open(fname, 'r')
        tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)
        csv_data = []
        for r, row in enumerate(tsvin):
            new = []
            for item in row:
                new.extend(item.split('\t'))
            csv_data.append(new)

        for r, row in enumerate(csv_data):
            for c, col in enumerate(row):
                csv_data[r][c] = col.strip()
            csv_data[r] = [_f for _f in csv_data[r] if _f]

        eventname = os.path.basename(os.path.splitext(fname)[0])

        eventparts = eventname.split('_')

        name = clean_snname(eventparts[0])
        name = catalog.add_event(name)
        secondaryname = 'CfA Supernova Archive'
        secondaryurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
        secondarysource = catalog.events[name].add_source(srcname=secondaryname,
                                                  url=secondaryurl,
                                                  secondary=True,
                                                  acknowledgment=ACKN_CFA)
        catalog.events[name].add_quantity('alias', name, secondarysource)

        year = re.findall(r'\d+', name)[0]
        catalog.events[name].add_quantity('discoverdate', year, secondarysource)

        eventbands = list(eventparts[1])

        tu = 'MJD'
        jdoffset = Decimal(0.)
        for rc, row in enumerate(csv_data):
            if len(row) > 0 and row[0][0] == "#":
                if len(row[0]) > 2 and row[0][:3] == '#JD':
                    tu = 'JD'
                    rowparts = row[0].split('-')
                    jdoffset = Decimal(rowparts[1])
                elif len(row[0]) > 6 and row[0][:7] == '#Julian':
                    tu = 'JD'
                    jdoffset = Decimal(0.)
                elif len(row) > 1 and row[1].lower() == 'photometry':
                    for ci, col in enumerate(row[2:]):
                        if col[0] == "(":
                            refstr = ' '.join(row[2 + ci:])
                            refstr = refstr.replace('(', '').replace(')', '')
                            bibcode = unescape(refstr)
                            source = catalog.events[name].add_source(bibcode=bibcode)
                elif len(row) > 1 and row[1] == 'HJD':
                    tu = 'HJD'
                continue

            elif len(row) > 0:
                mjd = row[0]
                for v, val in enumerate(row):
                    if v == 0:
                        if tu == 'JD':
                            mjd = str(jd_to_mjd(Decimal(val) + jdoffset))
                            tuout = 'MJD'
                        elif tu == 'HJD':
                            mjd = str(jd_to_mjd(Decimal(val)))
                            tuout = 'MJD'
                        else:
                            mjd = val
                            tuout = tu
                    elif v % 2 != 0:
                        if float(row[v]) < 90.0:
                            src = secondarysource + ',' + source
                            add_photometry(
                                events, name, u_time=tuout, time=mjd,
                                band=eventbands[(v - 1) // 2],
                                magnitude=row[v], e_magnitude=row[v + 1],
                                source=src)
        f.close()

    # Hicken 2012
    with open(os.path.join(PATH.REPO_EXTERNAL,
                           'hicken-2012-standard.dat'), 'r') as infile:
        tsvin = csv.reader(infile, delimiter='|', skipinitialspace=True)
        for r, row in enumerate(pbar(tsvin, current_task)):
            if r <= 47:
                continue

            if row[0][:2] != 'sn':
                name = 'SN' + row[0].strip()
            else:
                name = row[0].strip()

            name = catalog.add_event(name)

            source = catalog.events[name].add_source(bibcode='2012ApJS..200...12H')
            catalog.events[name].add_quantity('alias', name, source)
            catalog.events[name].add_quantity('claimedtype', 'Ia', source)
            add_photometry(
                events, name, u_time='MJD', time=row[2].strip(),
                band=row[1].strip(),
                magnitude=row[6].strip(), e_magnitude=row[7].strip(),
                source=source)

        # Bianco 2014
        tsvin = open(os.path.join(PATH.REPO_EXTERNAL,
                                  'bianco-2014-standard.dat'), 'r')
        tsvin = csv.reader(tsvin, delimiter=' ', skipinitialspace=True)
        for row in pbar(tsvin, current_task):
            name = 'SN' + row[0]
            name = catalog.add_event(name)

            source = catalog.events[name].add_source(bibcode='2014ApJS..213...19B')
            catalog.events[name].add_quantity('alias', name, source)
            add_photometry(
                events, name, u_time='MJD', time=row[2], band=row[1],
                magnitude=row[3],
                e_magnitude=row[4], telescope=row[5], system='Standard',
                source=source)

    catalog.journal_events()
    return


def do_cfa_spectra(catalog):
    current_task = task_obj.current_task(args)
    # Ia spectra
    oldname = ''
    file_names = next(os.walk(os.path.join(
        PATH.REPO_EXTERNAL_SPECTRA, 'CfA_SNIa')))[1]
    for name in pbar_strings(file_names, current_task):
        fullpath = os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'CfA_SNIa/') + name
        origname = name
        if name.startswith('sn') and is_number(name[2:6]):
            name = 'SN' + name[2:]
        if name.startswith('snf') and is_number(name[3:7]):
            name = 'SNF' + name[3:]
        name = get_preferred_name(events, name)
        if oldname and name != oldname:
            events = Events.journal_events(
                tasks, args, events, log)
        oldname = name
        name = catalog.add_event(name)
        reference = 'CfA Supernova Archive'
        refurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
        source = catalog.events[name].add_source(
            srcname=reference, url=refurl, secondary=True,
            acknowledgment=ACKN_CFA)
        catalog.events[name].add_quantity('alias', name, source)
        for fi, fname in enumerate(sorted(glob(fullpath + '/*'),
                                          key=lambda s: s.lower())):
            filename = os.path.basename(fname)
            fileparts = filename.split('-')
            if origname.startswith('sn') and is_number(origname[2:6]):
                year = fileparts[1][:4]
                month = fileparts[1][4:6]
                day = fileparts[1][6:]
                instrument = fileparts[2].split('.')[0]
            else:
                year = fileparts[2][:4]
                month = fileparts[2][4:6]
                day = fileparts[2][6:]
                instrument = fileparts[3].split('.')[0]
            time = str(astrotime(year + '-' + month + '-' +
                                 str(floor(float(day))).zfill(2)).mjd +
                       float(day) - floor(float(day)))
            f = open(fname, 'r')
            data = csv.reader(f, delimiter=' ', skipinitialspace=True)
            data = [list(i) for i in zip(*data)]
            wavelengths = data[0]
            fluxes = data[1]
            errors = data[2]
            sources = uniq_cdl([source,
                                (events[name]
                                 .add_source(bibcode='2012AJ....143..126B')),
                                (events[name]
                                 .add_source(bibcode='2008AJ....135.1598M'))])
            add_spectrum(
                events, name, 'Angstrom', 'erg/s/cm^2/Angstrom',
                filename=filename,
                wavelengths=wavelengths, fluxes=fluxes, u_time='MJD' if time
                else '', time=time, instrument=instrument,
                errorunit='ergs/s/cm^2/Angstrom', errors=errors,
                source=sources, dereddened=False, deredshifted=False)
            if args.travis and fi >= TRAVIS_QUERY_LIMIT:
                break
    catalog.journal_events()

    # Ibc spectra
    oldname = ''
    file_names = next(os.walk(os.path.join(
        PATH.REPO_EXTERNAL_SPECTRA, 'CfA_SNIbc')))[1]
    for name in pbar(file_names, current_task):
        fullpath = os.path.join(
            PATH.REPO_EXTERNAL_SPECTRA, 'CfA_SNIbc/') + name
        if name.startswith('sn') and is_number(name[2:6]):
            name = 'SN' + name[2:]
        name = get_preferred_name(events, name)
        if oldname and name != oldname:
            events = Events.journal_events(
                tasks, args, events, log)
        oldname = name
        name = catalog.add_event(name)
        reference = 'CfA Supernova Archive'
        refurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
        source = catalog.events[name].add_source(
            srcname=reference, url=refurl, secondary=True,
            acknowledgment=ACKN_CFA)
        catalog.events[name].add_quantity('alias', name, source)
        for fi, fname in enumerate(sorted(glob(fullpath + '/*'),
                                          key=lambda s: s.lower())):
            filename = os.path.basename(fname)
            fileparts = filename.split('-')
            instrument = ''
            year = fileparts[1][:4]
            month = fileparts[1][4:6]
            day = fileparts[1][6:].split('.')[0]
            if len(fileparts) > 2:
                instrument = fileparts[-1].split('.')[0]
            time = str(astrotime(year + '-' + month + '-' +
                                 str(floor(float(day))).zfill(2)).mjd +
                       float(day) - floor(float(day)))
            f = open(fname, 'r')
            data = csv.reader(f, delimiter=' ', skipinitialspace=True)
            data = [list(i) for i in zip(*data)]
            wavelengths = data[0]
            fluxes = data[1]
            sources = uniq_cdl(
                [source,
                 catalog.events[name].add_source(bibcode='2014AJ....147...99M')])
            add_spectrum(
                events, name, 'Angstrom', 'erg/s/cm^2/Angstrom',
                wavelengths=wavelengths, filename=filename,
                fluxes=fluxes, u_time='MJD' if time else '', time=time,
                instrument=instrument, source=sources,
                dereddened=False, deredshifted=False)
            if args.travis and fi >= TRAVIS_QUERY_LIMIT:
                break
    catalog.journal_events()

    # Other spectra
    oldname = ''
    file_names = next(os.walk(os.path.join(
        PATH.REPO_EXTERNAL_SPECTRA, 'CfA_Extra')))[1]
    for name in pbar_strings(file_names, current_task):
        fullpath = os.path.join(
            PATH.REPO_EXTERNAL_SPECTRA, 'CfA_Extra/') + name
        if name.startswith('sn') and is_number(name[2:6]):
            name = 'SN' + name[2:]
        name = get_preferred_name(events, name)
        if oldname and name != oldname:
            events = Events.journal_events(
                tasks, args, events, log)
        oldname = name
        name = catalog.add_event(name)
        reference = 'CfA Supernova Archive'
        refurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
        source = catalog.events[name].add_source(
            srcname=reference, url=refurl, secondary=True,
            acknowledgment=ACKN_CFA)
        catalog.events[name].add_quantity('alias', name, source)
        for fi, fname in enumerate(sorted(glob(fullpath + '/*'), key=lambda s:
                                          s.lower())):
            if not os.path.isfile(fname):
                continue
            filename = os.path.basename(fname)
            if ((not filename.startswith('sn') or
                 not filename.endswith('flm') or
                 any(x in filename for x in
                     ['-interp', '-z', '-dered', '-obj', '-gal']))):
                continue
            fileparts = filename.split('.')[0].split('-')
            instrument = ''
            time = ''
            if len(fileparts) > 1:
                year = fileparts[1][:4]
                month = fileparts[1][4:6]
                day = fileparts[1][6:]
                if is_number(year) and is_number(month) and is_number(day):
                    if len(fileparts) > 2:
                        instrument = fileparts[-1]
                    time = str(astrotime(year + '-' + month + '-' +
                                         str(floor(float(day))).zfill(2)).mjd +
                               float(day) - floor(float(day)))
            f = open(fname, 'r')
            data = csv.reader(f, delimiter=' ', skipinitialspace=True)
            data = [list(i) for i in zip(*data)]
            wavelengths = data[0]
            fluxes = [str(Decimal(x) * Decimal(1.0e-15)) for x in data[1]]
            add_spectrum(
                events, name, 'Angstrom', 'erg/s/cm^2/Angstrom',
                wavelengths=wavelengths, filename=filename,
                fluxes=fluxes, u_time='MJD' if time else '', time=time,
                instrument=instrument, source=source,
                dereddened=False, deredshifted=False)
            if args.travis and fi >= TRAVIS_QUERY_LIMIT:
                break

    catalog.journal_events()
    return
