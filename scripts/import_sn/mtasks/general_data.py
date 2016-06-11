"""General data import tasks.
"""
from astropy.time import Time as astrotime
from bs4 import BeautifulSoup
from collections import OrderedDict
import csv
from datetime import datetime
from glob import glob
from html import unescape
import os
import re
import requests
import urllib

from .. scripts import PATH
from ... utils import Decimal, is_number, pbar, pbar_strings
from .. funcs import add_event, add_photometry, add_source, add_quantity, archived_task, \
    jd_to_mjd, journal_events, load_event_from_file, make_date_string


def do_ascii(events, args, tasks):
    current_task = 'ASCII'

    # 2006ApJ...645..841N
    with open(os.path.join(PATH.REPO_EXTERNAL, '2006ApJ...645..841N-table3.csv'), 'r') as f:
        tsvin = csv.reader(f, delimiter=',')
        for ri, row in enumerate(pbar(tsvin, current_task)):
            name = 'SNLS-' + row[0]
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2006ApJ...645..841N')
            add_quantity(events, name, 'alias', name, source)
            add_quantity(events, name, 'redshift', row[1], source, kind='spectroscopic')
            astrot = astrotime(float(row[4]) + 2450000., format='jd').datetime
            date_str = make_date_string(astrot.year, astrot.month, astrot.day)
            add_quantity(events, name, 'discoverdate', date_str, source)
    events = journal_events(tasks, args, events)

    # Anderson 2014
    file_names = glob(os.path.join(PATH.REPO_EXTERNAL, 'SNII_anderson2014/*.dat'))
    for datafile in pbar_strings(file_names, desc=current_task):
        basename = os.path.basename(datafile)
        if not is_number(basename[:2]):
            continue
        if basename == '0210_V.dat':
            name = 'SN0210'
        else:
            name = ('SN20' if int(basename[:2]) < 50 else 'SN19') + basename.split('_')[0]
        name = add_event(tasks, args, events, name)
        source = add_source(events, name, bibcode='2014ApJ...786...67A')
        add_quantity(events, name, 'alias', name, source)

        if name in ['SN1999ca', 'SN2003dq', 'SN2008aw']:
            system = 'Swope'
        else:
            system = 'Landolt'

        with open(datafile, 'r') as f:
            tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)
            for row in tsvin:
                if not row[0]:
                    continue
                time = str(jd_to_mjd(Decimal(row[0])))
                add_photometry(
                    events, name, time=time, band='V', magnitude=row[1], e_magnitude=row[2],
                    system=system, source=source)
    events = journal_events(tasks, args, events)

    # stromlo
    stromlobands = ['B', 'V', 'R', 'I', 'VM', 'RM']
    with open(os.path.join(PATH.REPO_EXTERNAL, 'J_A+A_415_863-1/photometry.csv'), 'r') as f:
        tsvin = csv.reader(f, delimiter=',')
        for row in pbar(tsvin, current_task):
            name = row[0]
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2004A&A...415..863G')
            add_quantity(events, name, 'alias', name, source)
            mjd = str(jd_to_mjd(Decimal(row[1])))
            for ri, ci in enumerate(range(2, len(row), 3)):
                if not row[ci]:
                    continue
                band = stromlobands[ri]
                upperlimit = True if (not row[ci+1] and row[ci+2]) else False
                e_upper_magnitude = str(abs(Decimal(row[ci+1]))) if row[ci+1] else ''
                e_lower_magnitude = str(abs(Decimal(row[ci+2]))) if row[ci+2] else ''
                teles = 'MSSSO 1.3m' if band in ['VM', 'RM'] else 'CTIO'
                instr = 'MaCHO' if band in ['VM', 'RM'] else ''
                add_photometry(
                    events, name, time=mjd, band=band, magnitude=row[ci],
                    e_upper_magnitude=e_upper_magnitude, e_lower_magnitude=e_lower_magnitude,
                    upperlimit=upperlimit, telescope=teles, instrument=instr, source=source)
    events = journal_events(tasks, args, events)

    # 2015MNRAS.449..451W
    with open(os.path.join(PATH.REPO_EXTERNAL, '2015MNRAS.449..451W.dat'), 'r') as f:
        data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace=True)
        for r, row in enumerate(pbar(data, current_task)):
            if r == 0:
                continue
            namesplit = row[0].split('/')
            name = namesplit[-1]
            if name.startswith('SN'):
                name = name.replace(' ', '')
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2015MNRAS.449..451W')
            add_quantity(events, name, 'alias', name, source)
            if len(namesplit) > 1:
                add_quantity(events, name, 'alias', namesplit[0], source)
            add_quantity(events, name, 'claimedtype', row[1], source)
            add_photometry(events, name, time=row[2], band=row[4], magnitude=row[3], source=source)
    events = journal_events(tasks, args, events)

    # 2016MNRAS.459.1039T
    with open(os.path.join(PATH.REPO_EXTERNAL, '2016MNRAS.459.1039T.tsv'), 'r') as f:
        data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace=True)
        name = add_event(tasks, args, events, 'LSQ13zm')
        source = add_source(events, name, bibcode='2016MNRAS.459.1039T')
        add_quantity(events, name, 'alias', name, source)
        for r, row in enumerate(pbar(data, current_task)):
            if row[0][0] == '#':
                bands = [x.replace('(err)', '') for x in row[3:-1]]
                continue
            mjd = row[1]
            mags = [re.sub(r'\([^)]*\)', '', x) for x in row[3:-1]]
            upps = [True if '>' in x else '' for x in mags]
            mags = [x.replace('>', '') for x in mags]
            errs = [x[x.find('(')+1:x.find(')')] if '(' in x else '' for x in row[3:-1]]
            for mi, mag in enumerate(mags):
                if not is_number(mag):
                    continue
                add_photometry(
                    events, name, time=mjd, band=bands[mi], magnitude=mag, e_magnitude=errs[mi],
                    instrument=row[-1], upperlimit=upps[mi], source=source)
    events = journal_events(tasks, args, events)

    # 2015ApJ...804...28G
    with open(os.path.join(PATH.REPO_EXTERNAL, '2015ApJ...804...28G.tsv'), 'r') as f:
        data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace=True)
        name = add_event(tasks, args, events, 'PS1-13arp')
        source = add_source(events, name, bibcode='2015ApJ...804...28G')
        add_quantity(events, name, 'alias', name, source)
        for r, row in enumerate(pbar(data, current_task)):
            if r == 0:
                continue
            mjd = row[1]
            mag = row[3]
            upp = True if '<' in mag else ''
            mag = mag.replace('<', '')
            err = row[4] if is_number(row[4]) else ''
            ins = row[5]
            add_photometry(
                events, name, time=mjd, band=row[0], magnitude=mag, e_magnitude=err,
                instrument=ins, upperlimit=upp, source=source)
    events = journal_events(tasks, args, events)

    # 2016ApJ...819...35A
    with open(os.path.join(PATH.REPO_EXTERNAL, '2016ApJ...819...35A.tsv'), 'r') as f:
        data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace=True)
        for r, row in enumerate(pbar(data, current_task)):
            if row[0][0] == '#':
                continue
            name = add_event(tasks, args, events, row[0])
            source = add_source(events, name, bibcode='2016ApJ...819...35A')
            add_quantity(events, name, 'alias', name, source)
            add_quantity(events, name, 'ra', row[1], source)
            add_quantity(events, name, 'dec', row[2], source)
            add_quantity(events, name, 'redshift', row[3], source)
            disc_date = datetime.strptime(row[4], '%Y %b %d').isoformat()
            disc_date = disc_date.split('T')[0].replace('-', '/')
            add_quantity(events, name, 'discoverdate', disc_date, source)
    events = journal_events(tasks, args, events)

    # 2014ApJ...784..105W
    with open(os.path.join(PATH.REPO_EXTERNAL, '2014ApJ...784..105W.tsv'), 'r') as f:
        data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace=True)
        for r, row in enumerate(pbar(data, current_task)):
            if row[0][0] == '#':
                continue
            name = add_event(tasks, args, events, row[0])
            source = add_source(events, name, bibcode='2014ApJ...784..105W')
            add_quantity(events, name, 'alias', name, source)
            mjd = row[1]
            band = row[2]
            mag = row[3]
            err = row[4]
            add_photometry(
                events, name, time=mjd, band=row[2], magnitude=mag, e_magnitude=err,
                instrument='WHIRC', telescope='WIYN 3.5 m', observatory='NOAO',
                system='WHIRC', source=source)
    events = journal_events(tasks, args, events)

    # 2012MNRAS.425.1007B
    with open(os.path.join(PATH.REPO_EXTERNAL, '2012MNRAS.425.1007B.tsv'), 'r') as f:
        data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace=True)
        for r, row in enumerate(pbar(data, current_task)):
            if row[0][0] == '#':
                bands = row[2:]
                continue
            name = add_event(tasks, args, events, row[0])
            source = add_source(events, name, bibcode='2012MNRAS.425.1007B')
            add_quantity(events, name, 'alias', name, source)
            mjd = row[1]
            mags = [x.split('±')[0].strip() for x in row[2:]]
            errs = [x.split('±')[1].strip() if '±' in x else '' for x in row[2:]]
            if row[0] == 'PTF09dlc':
                ins = 'HAWK-I'
                tel = 'VLT 8.1m'
                obs = 'ESO'
            else:
                ins = 'NIRI'
                tel = 'Gemini North 8.2m'
                obs = 'Gemini'

            for mi, mag in enumerate(mags):
                if not is_number(mag):
                    continue
                add_photometry(
                    events, name, time=mjd, band=bands[mi], magnitude=mag, e_magnitude=errs[mi],
                    instrument=ins, telescope=tel, observatory=obs,
                    system='Natural', source=source)

    events = journal_events(tasks, args, events)
    return events


def do_cccp(events, args, tasks):
    current_task = 'CCCP'
    cccpbands = ['B', 'V', 'R', 'I']
    file_names = glob(os.path.join(PATH.REPO_EXTERNAL, 'CCCP/apj407397*.txt'))
    for datafile in pbar_strings(file_names, current_task + ': apj407397...'):
        with open(datafile, 'r') as f:
            tsvin = csv.reader(f, delimiter='\t', skipinitialspace=True)
            for r, row in enumerate(tsvin):
                if r == 0:
                    continue
                elif r == 1:
                    name = 'SN' + row[0].split('SN ')[-1]
                    name = add_event(tasks, args, events, name)
                    source = add_source(events, name, bibcode='2012ApJ...744...10K')
                    add_quantity(events, name, 'alias', name, source)
                elif r >= 5:
                    mjd = str(Decimal(row[0]) + 53000)
                    for b, band in enumerate(cccpbands):
                        if row[2*b + 1]:
                            mag = row[2*b + 1].strip('>')
                            upl = (not row[2*b + 2])
                            add_photometry(
                                events, name, time=mjd, band=band, magnitude=mag,
                                e_magnitude=row[2*b + 2], upperlimit=upl, source=source)

    if archived_task(tasks, args, 'cccp'):
        with open(os.path.join(PATH.REPO_EXTERNAL, 'CCCP/sc_cccp.html'), 'r') as f:
            html = f.read()
    else:
        session = requests.Session()
        response = session.get('https://webhome.weizmann.ac.il/home/iair/sc_cccp.html')
        html = response.text
        with open(os.path.join(PATH.REPO_EXTERNAL, 'CCCP/sc_cccp.html'), 'w') as f:
            f.write(html)

    soup = BeautifulSoup(html, 'html5lib')
    links = soup.body.findAll("a")
    for link in pbar(links, current_task + ': links'):
        if 'sc_sn' in link['href']:
            name = add_event(tasks, args, events, link.text.replace(' ', ''))
            source = add_source(events, name, refname='CCCP',
                                url='https://webhome.weizmann.ac.il/home/iair/sc_cccp.html')
            add_quantity(events, name, 'alias', name, source)

            if archived_task(tasks, args, 'cccp'):
                fname = os.path.join(PATH.REPO_EXTERNAL, 'CCCP/') + link['href'].split('/')[-1]
                with open(fname, 'r') as f:
                    html2 = f.read()
            else:
                response2 = session.get('https://webhome.weizmann.ac.il/home/iair/' + link['href'])
                html2 = response2.text
                fname = os.path.join(PATH.REPO_EXTERNAL, 'CCCP/') + link['href'].split('/')[-1]
                with open(fname, 'w') as f:
                    f.write(html2)

            soup2 = BeautifulSoup(html2, 'html5lib')
            links2 = soup2.body.findAll("a")
            for link2 in links2:
                if '.txt' in link2['href'] and '_' in link2['href']:
                    band = link2['href'].split('_')[1].split('.')[0].upper()
                    if archived_task(tasks, args, 'cccp'):
                        fname = os.path.join(PATH.REPO_EXTERNAL, 'CCCP/')
                        fname += link2['href'].split('/')[-1]
                        if not os.path.isfile(fname):
                            continue
                        with open(fname, 'r') as f:
                            html3 = f.read()
                    else:
                        response3 = session.get('https://webhome.weizmann.ac.il/home/iair/cccp/' +
                                                link2['href'])
                        if response3.status_code == 404:
                            continue
                        html3 = response3.text
                        fname = os.path.join(PATH.REPO_EXTERNAL, 'CCCP/')
                        fname += link2['href'].split('/')[-1]
                        with open(fname, 'w') as f:
                            f.write(html3)
                    table = [[str(Decimal(y.strip())).rstrip('0') for y in x.split(',')]
                             for x in list(filter(None, html3.split('\n')))]
                    for row in table:
                        add_photometry(
                            events, name, time=str(Decimal(row[0]) + 53000), band=band,
                            magnitude=row[1], e_magnitude=row[2], source=source)

    events = journal_events(tasks, args, events)
    return events


def do_external_radio(events, args, tasks):
    current_task = 'External Radio'
    path_pattern = os.path.join(PATH.REPO_EXTERNAL_RADIO, '*.txt')
    for datafile in pbar_strings(glob(path_pattern), desc=current_task):
        name = add_event(tasks, args, events, os.path.basename(datafile).split('.')[0])
        radiosourcedict = OrderedDict()
        with open(datafile, 'r') as f:
            for li, line in enumerate([x.strip() for x in f.read().splitlines()]):
                if line.startswith('(') and li <= len(radiosourcedict):
                    key = line.split()[0]
                    bibc = line.split()[-1]
                    radiosourcedict[key] = add_source(events, name, bibcode=bibc)
                elif li in [x + len(radiosourcedict) for x in range(3)]:
                    continue
                else:
                    cols = list(filter(None, line.split()))
                    source = radiosourcedict[cols[6]]
                    add_photometry(
                        events, name, time=cols[0], frequency=cols[2], u_frequency='GHz',
                        fluxdensity=cols[3], e_fluxdensity=cols[4], u_fluxdensity='µJy',
                        instrument=cols[5], source=source)
                    add_quantity(events, name, 'alias', name, source)

    events = journal_events(tasks, args, events)
    return events


def do_external_suspect_photometry(events, args, tasks):
    current_task = 'External: Suspect - Photometry'
    with open(os.path.join(PATH.REPO_EXTERNAL, 'suspectreferences.csv'), 'r') as f:
        tsvin = csv.reader(f, delimiter=',', skipinitialspace=True)
        suspectrefdict = {}
        for row in tsvin:
            suspectrefdict[row[0]] = row[1]

    file_names = glob(os.path.join(PATH.REPO_EXTERNAL, 'SUSPECT/*.html'))
    for datafile in pbar_strings(file_names, desc=current_task):
        basename = os.path.basename(datafile)
        basesplit = basename.split('-')
        name = basesplit[1]
        name = add_event(tasks, args, events, name)
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
        source = add_source(events, name, bibcode=bibcode)

        secondaryreference = 'SUSPECT'
        secondaryrefurl = 'https://www.nhn.ou.edu/~suspect/'
        secondarysource = add_source(events, name, refname=secondaryreference, url=secondaryrefurl,
                                     secondary=True)
        add_quantity(events, name, 'alias', name, secondarysource)

        if ei == 1:
            year = re.findall(r'\d+', name)[0]
            add_quantity(events, name, 'discoverdate', year, secondarysource)
            add_quantity(events, name, 'host', names[1].split(':')[1].strip(), secondarysource)

            redshifts = bandsoup.body.findAll(text=re.compile('Redshift'))
            if redshifts:
                add_quantity(
                    events, name, 'redshift', redshifts[0].split(':')[1].strip(),
                    secondarysource, kind='heliocentric')
            # hvels = bandsoup.body.findAll(text=re.compile('Heliocentric Velocity'))
            # if hvels:
            #    add_quantity(events, name, 'velocity', hvels[0].split(':')[1].strip().split(' ')[0],
            #        secondarysource, kind='heliocentric')
            types = bandsoup.body.findAll(text=re.compile('Type'))

            add_quantity(
                events, name, 'claimedtype', types[0].split(':')[1].strip().split(' ')[0],
                secondarysource)

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
            add_photometry(
                events, name, time=mjd, band=band, magnitude=mag, e_magnitude=e_magnitude,
                source=secondarysource + ',' + source)

    events = journal_events(tasks, args, events)
    return events


def do_external_suspect_spectra(events, args, tasks):
    import json
    from .. constants import TRAVIS_QUERY_LIMIT
    from .. funcs import add_spectrum, get_preferred_name, uniq_cdl
    from math import floor
    from ... utils import get_sig_digits, pretty_num
    current_task = 'External: Suspect - Spectra'
    with open(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'Suspect/sources.json'), 'r') as f:
        sourcedict = json.loads(f.read())

    with open(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'Suspect/filename-changes.txt'), 'r') as f:
        rows = f.readlines()
        changedict = {}
        for row in rows:
            if not row.strip() or row[0] == "#":
                continue
            items = row.strip().split(' ')
            changedict[items[1]] = items[0]

    suspectcnt = 0
    folders = next(os.walk(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'Suspect')))[1]
    for folder in pbar(folders, current_task):
        eventfolders = next(os.walk(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'Suspect/')+folder))[1]
        oldname = ''
        for eventfolder in pbar(eventfolders, current_task):
            name = eventfolder
            if is_number(name[:4]):
                name = 'SN' + name
            name = get_preferred_name(events, name)
            if oldname and name != oldname:
                events = journal_events(tasks, args, events)
            oldname = name
            name = add_event(tasks, args, events, name)
            sec_ref = 'SUSPECT'
            sec_refurl = 'https://www.nhn.ou.edu/~suspect/'
            sec_bibc = '2001AAS...199.8408R'
            sec_source = add_source(
                events, name, refname=sec_ref, url=sec_refurl, bibcode=sec_bibc, secondary=True)
            add_quantity(events, name, 'alias', name, sec_source)
            fpath = os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'Suspect', folder, eventfolder)
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
                    source = add_source(events, name, bibcode=unescape(bibcode))
                    sources += source
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

                fpath = os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'Suspect', folder,
                                     eventfolder, spectrum)
                with open() as f:
                    specdata = list(csv.reader(f, delimiter=' ', skipinitialspace=True))
                    specdata = list(filter(None, specdata))
                    newspec = []
                    oldval = ''
                    for row in specdata:
                        if row[1] == oldval:
                            continue
                        newspec.append(row)
                        oldval = row[1]
                    specdata = newspec
                haserrors = len(specdata[0]) == 3 and specdata[0][2] and specdata[0][2] != 'NaN'
                specdata = [list(i) for i in zip(*specdata)]

                wavelengths = specdata[0]
                fluxes = specdata[1]
                errors = ''
                if haserrors:
                    errors = specdata[2]

                add_spectrum(
                    name=name, u_time='MJD', time=time, waveunit='Angstrom',
                    fluxunit='Uncalibrated', wavelengths=wavelengths,
                    fluxes=fluxes, errors=errors, errorunit='Uncalibrated',
                    source=sources, filename=spectrum)
                suspectcnt = suspectcnt + 1
                if args.travis and suspectcnt % TRAVIS_QUERY_LIMIT == 0:
                    break

    events = journal_events(tasks, args, events)
    return events


def do_external_xray(events, args, tasks):
    current_task = 'External X-ray'
    path_pattern = os.path.join(PATH.REPO_EXTERNAL_XRAY, '*.txt')
    for datafile in pbar_strings(glob(path_pattern), desc=current_task):
        name = add_event(tasks, args, events, os.path.basename(datafile).split('.')[0])
        with open(datafile, 'r') as f:
            for li, line in enumerate(f.read().splitlines()):
                if li == 0:
                    source = add_source(events, name, bibcode=line.split()[-1])
                elif li in [1, 2, 3]:
                    continue
                else:
                    cols = list(filter(None, line.split()))
                    add_photometry(
                        events, name, time=cols[:2],
                        energy=cols[2:4], u_energy='keV', counts=cols[4], flux=cols[6],
                        unabsorbedflux=cols[8], u_flux='ergs/s/cm^2',
                        photonindex=cols[15], instrument=cols[17], nhmw=cols[11],
                        upperlimit=(float(cols[5]) < 0), source=source)
                    add_quantity(events, name, 'alias', name, source)

    events = journal_events(tasks, args, events)
    return events


def do_internal(events, args, tasks):
    """Load events from files in the 'internal' repository, and save them.
    """
    current_task = 'Internal'
    path_pattern = os.path.join(PATH.REPO_INTERNAL, '*.json')
    files = glob(path_pattern)
    for datafile in pbar_strings(files, desc=current_task):
        if args.update:
            if not load_event_from_file(events, args, tasks, path=datafile,
                                        clean=True, delete=False, append=True):
                raise IOError('Failed to find specified file.')
        else:
            if not load_event_from_file(events, args, tasks, path=datafile,
                                        clean=True, delete=False):
                raise IOError('Failed to find specified file.')

    events = journal_events(tasks, args, events)
    return events


'''
def do_simbad(events, args, tasks):
    Simbad.list_votable_fields()
    customSimbad = Simbad()
    customSimbad.add_votable_fields('otype', 'id(opt)')
    result = customSimbad.query_object('SN 20[0-9][0-9]*', wildcard=True)
    for r, row in enumerate(result):
        if row['OTYPE'].decode() != 'SN':
            continue
        name = row['MAIN_ID'].decode()
        aliases = Simbad.query_objectids(name)
        print(aliases)
        if name[:3] == 'SN ':
            name = 'SN' + name[3:]
        if name[:2] == 'SN' and is_number(name[2:]):
            name = name + 'A'
        name = add_event(tasks, args, events, name)
    events = journal_events(tasks, args, events)
    return events
'''
