"""Data import tasks for individual data 'donations'.
"""
import csv
import json
import os
from glob import glob
from math import isnan

from scripts import PATH
from ... utils import pbar, is_number, pbar_strings
from .. funcs import add_event, add_photometry, add_source, add_quantity, \
    journal_events


def do_donations(events, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    # Nicholl 04-01-16 donation
    with open(os.path.join(PATH.REPO_EXTERNAL, 'Nicholl-04-01-16/bibcodes.json'), 'r') as f:
        bcs = json.loads(f.read())

    file_names = glob(os.path.join(PATH.REPO_EXTERNAL, 'Nicholl-04-01-16/*.txt'))
    for datafile in pbar_strings(file_names, current_task + ': Nicholl-04-01-16'):
        name = os.path.basename(datafile).split('_')[0]
        events, name = add_event(tasks, args, events, name, log)
        bibcode = ''
        for bc in bcs:
            if name in bcs[bc]:
                bibcode = bc
        if not bibcode:
            raise ValueError('Bibcode not found!')
        source = add_source(events, name, bibcode=bibcode)
        add_quantity(events, name, 'alias', name, source)
        with open(datafile, 'r') as f:
            tsvin = csv.reader(f, delimiter='\t', skipinitialspace=True)
            for r, rrow in enumerate(tsvin):
                row = list(filter(None, rrow))
                if not row:
                    continue
                if row[0][0] == '#' and row[0] != '#MJD':
                    continue
                if row[0] == '#MJD':
                    bands = [x for x in row[1:] if x and 'err' not in x]
                    continue
                mjd = row[0]
                if not is_number(mjd):
                    continue
                for v, val in enumerate(row[1::2]):
                    upperlimit = ''
                    if '>' in val:
                        upperlimit = True
                    mag = val.strip('>')
                    if not is_number(mag) or isnan(float(mag)) or float(mag) > 90.0:
                        continue
                    err = ''
                    if is_number(row[2*v+2]) and not isnan(float(row[2*v+2])):
                        err = row[2*v+2]
                    add_photometry(
                        events, name, time=mjd, band=bands[v], magnitude=mag,
                        e_magnitude=err, upperlimit=upperlimit, source=source)
    events = journal_events(tasks, args, events)

    # Maggi 04-11-16 donation (MC SNRs)
    with open(os.path.join(PATH.REPO_EXTERNAL, 'Maggi-04-11-16/LMCSNRs_OpenSNe.csv')) as f:
        tsvin = csv.reader(f, delimiter=',')
        for row in pbar(tsvin, current_task + ': Maggi-04-11-16/LMCSNRs'):
            name = 'MCSNR ' + row[0]
            events, name = add_event(tasks, args, events, name, log)
            source = add_source(events, name, bibcode='2016A&A...585A.162M')
            add_quantity(events, name, 'alias', name, source)
            if row[1] != 'noname':
                add_quantity(events, name, 'alias', row[1], source)
            add_quantity(events, name, 'ra', row[2], source)
            add_quantity(events, name, 'dec', row[3], source)
            add_quantity(events, name, 'host', 'LMC', source)
            if row[4] == '1':
                add_quantity(events, name, 'claimedtype', 'Ia', source)
            elif row[4] == '2':
                add_quantity(events, name, 'claimedtype', 'CC', source)
    with open(os.path.join(PATH.REPO_EXTERNAL, 'Maggi-04-11-16/SMCSNRs_OpenSNe.csv')) as f:
        tsvin = csv.reader(f, delimiter=',')
        for row in pbar(tsvin, current_task + ': Maggi-04-11-16/SMCSNRs'):
            name = 'MCSNR ' + row[0]
            events, name = add_event(tasks, args, events, name, log)
            source = add_source(events, name, srcname='Pierre Maggi')
            add_quantity(events, name, 'alias', name, source)
            add_quantity(events, name, 'alias', row[1], source)
            add_quantity(events, name, 'alias', row[2], source)
            add_quantity(events, name, 'ra', row[3], source)
            add_quantity(events, name, 'dec', row[4], source)
            add_quantity(events, name, 'host', 'SMC', source)
    events = journal_events(tasks, args, events)

    # Galbany 04-18-16 donation
    folders = next(os.walk(os.path.join(PATH.REPO_EXTERNAL, 'galbany-04-18-16/')))[1]
    bibcode = '2016AJ....151...33G'
    for folder in folders:
        infofiles = glob(os.path.join(PATH.REPO_EXTERNAL, 'galbany-04-18-16/') + folder + '/*.info')
        photfiles = glob(os.path.join(PATH.REPO_EXTERNAL, 'galbany-04-18-16/') + folder + '/*.out*')

        zhel = ''
        zcmb = ''
        zerr = ''
        for path in infofiles:
            with open(path, 'r') as f:
                lines = f.read().splitlines()
                for line in lines:
                    splitline = line.split(':')
                    field = splitline[0].strip().lower()
                    value = splitline[1].strip()
                    if field == 'name':
                        name = value[:6].upper()
                        name += (value[6].upper() if len(value) == 7 else value[6:])
                        events, name = add_event(tasks, args, events, name, log)
                        source = add_source(events, name, bibcode=bibcode)
                        add_quantity(events, name, 'alias', name, source)
                    elif field == 'type':
                        claimedtype = value.replace('SN', '')
                        add_quantity(events, name, 'claimedtype', claimedtype, source)
                    elif field == 'zhel':
                        zhel = value
                    elif field == 'redshift_error':
                        zerr = value
                    elif field == 'zcmb':
                        zcmb = value
                    elif field == 'ra':
                        add_quantity(events, name, 'ra', value, source, unit='floatdegrees')
                    elif field == 'dec':
                        add_quantity(events, name, 'dec', value, source, unit='floatdegrees')
                    elif field == 'host':
                        value = value.replace('- ', '-').replace('G ', 'G')
                        add_quantity(events, name, 'host', value, source)
                    elif field == 'e(b-v)_mw':
                        add_quantity(events, name, 'ebv', value, source)

        add_quantity(events, name, 'redshift', zhel, source, error=zerr, kind='heliocentric')
        add_quantity(events, name, 'redshift', zcmb, source, error=zerr, kind='cmb')

        for path in photfiles:
            with open(path, 'r') as f:
                band = ''
                lines = f.read().splitlines()
                for li, line in enumerate(lines):
                    if li in [0, 2, 3]:
                        continue
                    if li == 1:
                        band = line.split(':')[-1].strip()
                    else:
                        cols = list(filter(None, line.split()))
                        if not cols:
                            continue
                        add_photometry(
                            events, name, time=cols[0], magnitude=cols[1], e_magnitude=cols[2],
                            band=band, system=cols[3], telescope=cols[4], source=source)
    events = journal_events(tasks, args, events)

    # Brown 05-14-16
    files = glob(os.path.join(PATH.REPO_EXTERNAL, 'brown-05-14-16/*.dat'))
    for fi in pbar(files, current_task):
        name = os.path.basename(fi).split('_')[0]
        events, name = add_event(tasks, args, events, name, log)
        source = add_source(
            events, name, refname='Swift Supernovae', bibcode='2014Ap&SS.354...89B',
            url='http://people.physics.tamu.edu/pbrown/SwiftSN/swift_sn.html')
        add_quantity(events, name, 'alias', name, source)
        with open(fi, 'r') as f:
            lines = f.read().splitlines()
            for line in lines:
                if not line or line[0] == '#':
                    continue
                cols = list(filter(None, line.split()))
                band = cols[0]
                mjd = cols[1]
                # Skip lower limit entries for now
                if cols[2] == 'NULL' and cols[6] == 'NULL':
                    continue
                isupp = cols[2] == 'NULL' and cols[6] != 'NULL'
                mag = cols[2] if not isupp else cols[4]
                e_mag = cols[3] if not isupp else ''
                upp = '' if not isupp else True
                add_photometry(events, name, time=mjd, magnitude=mag, e_magnitude=e_mag,
                               upperlimit=upp, band=band, source=source,
                               telescope='Swift', instrument='UVOT', system='Vega')
    events = journal_events(tasks, args, events)

    # Nicholl 05-03-16
    files = glob(os.path.join(PATH.REPO_EXTERNAL, 'nicholl-05-03-16/*.txt'))
    events, name = add_event(tasks, args, events, 'SN2015bn', log)
    source = add_source(events, name, bibcode='2016arXiv160304748N')
    add_quantity(events, name, 'alias', name, source)
    add_quantity(events, name, 'alias', 'PS15ae', source)
    for fi in pbar(files, current_task):
        telescope = os.path.basename(fi).split('_')[1]
        with open(fi, 'r') as f:
            lines = f.read().splitlines()
            for li, line in enumerate(lines):
                if not line or (line[0] == '#' and li != 0):
                    continue
                cols = list(filter(None, line.split()))
                if not cols:
                    continue
                if li == 0:
                    bands = cols[1:]
                    continue

                mjd = cols[0]
                for ci, col in enumerate(cols[1::2]):
                    if not is_number(col):
                        continue

                    emag = cols[2*ci+2]
                    upp = ''
                    if not is_number(emag):
                        emag = ''
                        upp = True
                    instrument = 'UVOT' if telescope == 'Swift' else ''
                    add_photometry(
                        events, name, time=mjd, magnitude=col, e_magnitude=emag, upperlimit=upp,
                        band=bands[ci], source=source, telescope=telescope, instrument=instrument)

    events = journal_events(tasks, args, events)
    return events
