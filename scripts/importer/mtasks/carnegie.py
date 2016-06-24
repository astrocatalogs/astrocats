"""Imports from the 'Vizier' catalog.
"""
from cdecimal import Decimal
import csv
from glob import glob
import os

from scripts import PATH
from .. import Events
from .. funcs import add_photometry, add_spectrum, \
    clean_snname, get_preferred_name, jd_to_mjd
from .. constants import TRAVIS_QUERY_LIMIT
from ... utils import pbar_strings


def do_csp_photo(events, stubs, args, tasks, task_obj, log):
    import re
    cspbands = ['u', 'B', 'V', 'g', 'r', 'i', 'Y', 'J', 'H', 'K']
    file_names = glob(os.path.join(PATH.REPO_EXTERNAL, 'CSP/*.dat'))
    current_task = task_obj.current_task(args)
    for fname in pbar_strings(file_names, desc=current_task):
        tsvin = csv.reader(open(fname, 'r'), delimiter='\t', skipinitialspace=True)
        eventname = os.path.basename(os.path.splitext(fname)[0])
        eventparts = eventname.split('opt+')
        name = clean_snname(eventparts[0])
        events, name = Events.add_event(tasks, args, events, name, log)

        reference = 'Carnegie Supernova Project'
        refbib = '2010AJ....139..519C'
        refurl = 'http://csp.obs.carnegiescience.edu/data'
        source = events[name].add_source(bibcode=refbib, srcname=reference, url=refurl)
        events[name].add_quantity('alias', name, source)

        year = re.findall(r'\d+', name)[0]
        events[name].add_quantity('discoverdate', year, source)

        for r, row in enumerate(tsvin):
            if len(row) > 0 and row[0][0] == "#":
                if r == 2:
                    redz = row[0].split(' ')[-1]
                    events[name].add_quantity('redshift', redz, source, kind='cmb')
                    events[name].add_quantity('ra', row[1].split(' ')[-1], source)
                    events[name].add_quantity('dec', row[2].split(' ')[-1], source)
                continue
            for v, val in enumerate(row):
                if v == 0:
                    mjd = val
                elif v % 2 != 0:
                    if float(row[v]) < 90.0:
                        add_photometry(
                            events, name, time=mjd, observatory='LCO', band=cspbands[(v-1)//2],
                            system='CSP', magnitude=row[v], e_magnitude=row[v+1], source=source)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_csp_spectra(events, stubs, args, tasks, task_obj, log):
    oldname = ''
    current_task = task_obj.current_task(args)
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
            events, stubs = Events.journal_events(tasks, args, events, stubs, log)
        oldname = name
        events, name = Events.add_event(tasks, args, events, name, log)
        telescope = fileparts[-2]
        instrument = fileparts[-1]
        source = events[name].add_source(bibcode='2013ApJ...773...53F')
        events[name].add_quantity('alias', name, source)

        data = csv.reader(open(fname, 'r'), delimiter=' ', skipinitialspace=True)
        specdata = []
        for r, row in enumerate(data):
            if row[0] == '#JDate_of_observation:':
                jd = row[1].strip()
                time = str(jd_to_mjd(Decimal(jd)))
            elif row[0] == '#Redshift:':
                events[name].add_quantity('redshift', row[1].strip(), source)
            if r < 7:
                continue
            specdata.append(list(filter(None, [x.strip(' ') for x in row])))
        specdata = [list(i) for i in zip(*specdata)]
        wavelengths = specdata[0]
        fluxes = specdata[1]

        add_spectrum(
            name=name, u_time='MJD', time=time, waveunit='Angstrom', fluxunit='erg/s/cm^2/Angstrom',
            wavelengths=wavelengths, fluxes=fluxes, telescope=telescope, instrument=instrument,
            source=source, deredshifted=True, filename=filename)
        if args.travis and fi >= TRAVIS_QUERY_LIMIT:
            break

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events
