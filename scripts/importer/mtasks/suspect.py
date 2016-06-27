"""General data import tasks.
"""
import csv
import json
import os
import re
import urllib
from glob import glob
from html import unescape
from math import floor

from astropy.time import Time as astrotime
from bs4 import BeautifulSoup

from cdecimal import Decimal
from scripts import PATH

from .. import Events
from ...utils import get_sig_digits, is_number, pbar, pbar_strings, pretty_num
from ..constants import TRAVIS_QUERY_LIMIT
from ..funcs import (add_photometry, add_spectrum, get_preferred_name,
                     jd_to_mjd, uniq_cdl)


def do_suspect_photo(catalog):
    current_task = task_obj.current_task(args)
    with open(os.path.join(PATH.REPO_EXTERNAL,
                           'suspectreferences.csv'), 'r') as f:
        tsvin = csv.reader(f, delimiter=',', skipinitialspace=True)
        suspectrefdict = {}
        for row in tsvin:
            suspectrefdict[row[0]] = row[1]

    file_names = glob(os.path.join(PATH.REPO_EXTERNAL, 'SUSPECT/*.html'))
    for datafile in pbar_strings(file_names, desc=current_task):
        basename = os.path.basename(datafile)
        basesplit = basename.split('-')
        oldname = basesplit[1]
        name = catalog.add_event(oldname)
        if name.startswith('SN') and is_number(name[2:]):
            name = name + 'A'
        band = basesplit[3].split('.')[0]
        ei = int(basesplit[2])
        bandlink = 'file://' + os.path.abspath(datafile)
        bandresp = urllib.request.urlopen(bandlink)
        bandsoup = BeautifulSoup(bandresp, 'html5lib')
        bandtable = bandsoup.find('table')

        names = bandsoup.body.findAll(text=re.compile('Name'))
        reference = ''
        for link in bandsoup.body.findAll('a'):
            if 'adsabs' in link['href']:
                reference = str(link).replace('"', "'")

        bibcode = unescape(suspectrefdict[reference])
        source = catalog.events[name].add_source(bibcode=bibcode)

        sec_ref = 'SUSPECT'
        sec_refurl = 'https://www.nhn.ou.edu/~suspect/'
        sec_source = catalog.events[name].add_source(
            srcname=sec_ref, url=sec_refurl, secondary=True)
        catalog.events[name].add_quantity('alias', oldname, sec_source)

        if ei == 1:
            year = re.findall(r'\d+', name)[0]
            catalog.events[name].add_quantity('discoverdate', year, sec_source)
            catalog.events[name].add_quantity('host', names[1].split(':')[
                                      1].strip(), sec_source)

            redshifts = bandsoup.body.findAll(text=re.compile('Redshift'))
            if redshifts:
                catalog.events[name].add_quantity(
                    'redshift', redshifts[0].split(':')[1].strip(),
                    sec_source, kind='heliocentric')
            # hvels = bandsoup.body.findAll(text=re.compile('Heliocentric
            # Velocity'))
            # if hvels:
            #     vel = hvels[0].split(':')[1].strip().split(' ')[0]
            #     catalog.events[name].add_quantity('velocity', vel, sec_source,
            # kind='heliocentric')
            types = bandsoup.body.findAll(text=re.compile('Type'))

            catalog.events[name].add_quantity(
                'claimedtype', types[0].split(':')[1].strip().split(' ')[0],
                sec_source)

        for r, row in enumerate(bandtable.findAll('tr')):
            if r == 0:
                continue
            col = row.findAll('td')
            mjd = str(jd_to_mjd(Decimal(col[0].contents[0])))
            mag = col[3].contents[0]
            if mag.isspace():
                mag = ''
            else:
                mag = str(mag)
            e_magnitude = col[4].contents[0]
            if e_magnitude.isspace():
                e_magnitude = ''
            else:
                e_magnitude = str(e_magnitude)
            add_photometry(events, name, time=mjd, band=band, magnitude=mag,
                           e_magnitude=e_magnitude,
                           source=sec_source + ',' + source)

    events = Events.journal_events(tasks, args, events, log)
    return events


def do_suspect_spectra(catalog):
    current_task = task_obj.current_task(args)
    with open(os.path.join(PATH.REPO_EXTERNAL_SPECTRA,
                           'Suspect/sources.json'), 'r') as f:
        sourcedict = json.loads(f.read())

    with open(os.path.join(PATH.REPO_EXTERNAL_SPECTRA,
                           'Suspect/filename-changes.txt'), 'r') as f:
        rows = f.readlines()
        changedict = {}
        for row in rows:
            if not row.strip() or row[0] == "#":
                continue
            items = row.strip().split(' ')
            changedict[items[1]] = items[0]

    suspectcnt = 0
    folders = next(os.walk(os.path.join(
        PATH.REPO_EXTERNAL_SPECTRA, 'Suspect')))[1]
    for folder in pbar(folders, current_task):
        eventfolders = next(os.walk(os.path.join(
            PATH.REPO_EXTERNAL_SPECTRA, 'Suspect/') + folder))[1]
        oldname = ''
        for eventfolder in pbar(eventfolders, current_task):
            name = eventfolder
            if is_number(name[:4]):
                name = 'SN' + name
            name = get_preferred_name(events, name)
            if oldname and name != oldname:
                events = Events.journal_events(
                    tasks, args, events, log)
            oldname = name
            name = catalog.add_event(name)
            sec_ref = 'SUSPECT'
            sec_refurl = 'https://www.nhn.ou.edu/~suspect/'
            sec_bibc = '2001AAS...199.8408R'
            sec_source = catalog.events[name].add_source(
                srcname=sec_ref, url=sec_refurl, bibcode=sec_bibc,
                secondary=True)
            catalog.events[name].add_quantity('alias', name, sec_source)
            fpath = os.path.join(PATH.REPO_EXTERNAL_SPECTRA,
                                 'Suspect', folder, eventfolder)
            eventspectra = next(os.walk(fpath))[2]
            for spectrum in eventspectra:
                sources = [sec_source]
                bibcode = ''
                if spectrum in changedict:
                    specalias = changedict[spectrum]
                else:
                    specalias = spectrum
                if specalias in sourcedict:
                    bibcode = sourcedict[specalias]
                elif name in sourcedict:
                    bibcode = sourcedict[name]
                if bibcode:
                    source = catalog.events[name].add_source(bibcode=unescape(bibcode))
                    sources += [source]
                sources = uniq_cdl(sources)

                date = spectrum.split('_')[1]
                year = date[:4]
                month = date[4:6]
                day = date[6:]
                sig = get_sig_digits(day) + 5
                day_fmt = str(floor(float(day))).zfill(2)
                time = astrotime(year + '-' + month + '-' + day_fmt).mjd
                time = time + float(day) - floor(float(day))
                time = pretty_num(time, sig=sig)

                fpath = os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'Suspect',
                                     folder,
                                     eventfolder, spectrum)
                with open() as f:
                    specdata = list(csv.reader(
                        f, delimiter=' ', skipinitialspace=True))
                    specdata = list(filter(None, specdata))
                    newspec = []
                    oldval = ''
                    for row in specdata:
                        if row[1] == oldval:
                            continue
                        newspec.append(row)
                        oldval = row[1]
                    specdata = newspec
                haserrors = len(specdata[0]) == 3 and specdata[
                    0][2] and specdata[0][2] != 'NaN'
                specdata = [list(i) for i in zip(*specdata)]

                wavelengths = specdata[0]
                fluxes = specdata[1]
                errors = ''
                if haserrors:
                    errors = specdata[2]

                add_spectrum(
                    events, name, 'Angstrom', 'Uncalibrated', u_time='MJD',
                    time=time,
                    wavelengths=wavelengths, fluxes=fluxes, errors=errors,
                    errorunit='Uncalibrated',
                    source=sources, filename=spectrum)
                suspectcnt = suspectcnt + 1
                if args.travis and suspectcnt % TRAVIS_QUERY_LIMIT == 0:
                    break

    events = Events.journal_events(tasks, args, events, log)
    return events
