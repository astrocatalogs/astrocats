"""General data import tasks.
"""
from astropy.time import Time as astrotime
from cdecimal import Decimal
import csv
from glob import glob
import os

from scripts import PATH, TRAVIS_QUERY_LIMIT
from .. import Events
from .. funcs import add_spectrum, get_preferred_name, jd_to_mjd, uniq_cdl
from ... utils import pretty_num


def do_snf_aliases(events, stubs, args, tasks, task_obj, log):
    with open('../sne-external/SNF/snf-aliases.csv') as f:
        for row in [x.split(',') for x in f.read().splitlines()]:
            events, name, source = Events.new_event(tasks, args, events, row[0], log,
                                                    bibcode = oscbibcode, refname = oscname,
                                                    url = oscurl, secondary = True)
            events[name].add_quantity('alias', row[1], source)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_snf_specta(events, stubs, args, tasks, task_obj, log):
    bibcodes = {'SN2005gj': '2006ApJ...650..510A', 'SN2006D': '2007ApJ...654L..53T',
                'SN2007if': '2010ApJ...713.1073S', 'SN2011fe': '2013A&A...554A..27P'}
    oldname = ''
    snfcnt = 0
    eventfolders = next(os.walk(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'SNFactory')))[1]
    for eventfolder in eventfolders:
        name = eventfolder
        name = get_preferred_name(events, name)
        if oldname and name != oldname:
            events, stubs = Events.journal_events(tasks, args, events, stubs, log)
        oldname = name
        events, name = Events.add_event(tasks, args, events, name, log)
        sec_reference = 'Nearby Supernova Factory'
        sec_refurl = 'http://snfactory.lbl.gov/'
        sec_bibcode = '2002SPIE.4836...61A'
        sec_source = events[name].add_source(
            refname=sec_reference, url=sec_refurl, bibcode=sec_bibcode, secondary=True)
        events[name].add_quantity('alias', name, sec_source)
        bibcode = bibcodes[name]
        source = events[name].add_source(bibcode=bibcode)
        sources = uniq_cdl([source, sec_source])
        use_path = os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'SNFactory', eventfolder, '*.dat')
        eventspectra = glob(use_path)
        for spectrum in eventspectra:
            filename = os.path.basename(spectrum)
            with open(spectrum) as spec_file:
                specdata = list(csv.reader(spec_file, delimiter=' ', skipinitialspace=True))
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
                raise ValueError('Time missing from spectrum.')
            specdata = newspec
            haserrors = len(specdata[0]) == 3 and specdata[0][2] and specdata[0][2] != 'NaN'
            specdata = [list(i) for i in zip(*specdata)]

            wavelengths = specdata[0]
            fluxes = specdata[1]
            errors = ''
            if haserrors:
                errors = specdata[2]

            unit_err = 'Variance' if name == 'SN2011fe' else 'erg/s/cm^2/Angstrom'
            unit_flx = 'erg/s/cm^2/Angstrom'
            add_spectrum(
                name=name, u_time='MJD', time=time, waveunit='Angstrom', fluxunit=unit_flx,
                wavelengths=wavelengths, fluxes=fluxes, errors=errors, observer=observer,
                observatory=observatory, telescope=telescope, instrument=instrument,
                errorunit=unit_err, source=sources, filename=filename)
            snfcnt = snfcnt + 1
            if args.travis and snfcnt % TRAVIS_QUERY_LIMIT == 0:
                break

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events
