"""General data import tasks.
"""
import csv
import os
import re
from datetime import datetime
from glob import glob

from astropy.time import Time as astrotime

from cdecimal import Decimal
from scripts import PATH
from scripts.utils import is_number, pbar, pbar_strings

from .. import Events
from ..funcs import add_photometry, jd_to_mjd, make_date_string


def do_ascii(catalog):
    current_task = catalog.current_task

    # 2006ApJ...645..841N
    file_path = os.path.join(
        PATH.REPO_EXTERNAL, '2006ApJ...645..841N-table3.csv')
    tsvin = list(csv.reader(open(file_path, 'r'), delimiter=','))
    for ri, row in enumerate(pbar(tsvin, current_task)):
        name = 'SNLS-' + row[0]
        name = catalog.add_event(name)
        source = catalog.events[name].add_source(bibcode='2006ApJ...645..841N')
        catalog.events[name].add_quantity('alias', name, source)
        catalog.events[name].add_quantity(
            'redshift', row[1], source, kind='spectroscopic')
        astrot = astrotime(float(row[4]) + 2450000., format='jd').datetime
        date_str = make_date_string(astrot.year, astrot.month, astrot.day)
        catalog.events[name].add_quantity('discoverdate', date_str, source)
    catalog.journal_events()

    # Anderson 2014
    file_names = list(
        glob(os.path.join(PATH.REPO_EXTERNAL, 'SNII_anderson2014/*.dat')))
    for datafile in pbar_strings(file_names, desc=current_task):
        basename = os.path.basename(datafile)
        if not is_number(basename[:2]):
            continue
        if basename == '0210_V.dat':
            name = 'SN0210'
        else:
            name = ('SN20' if int(basename[:2]) <
                    50 else 'SN19') + basename.split('_')[0]
        name = catalog.add_event(name)
        source = catalog.events[name].add_source(bibcode='2014ApJ...786...67A')
        catalog.events[name].add_quantity('alias', name, source)

        if name in ['SN1999ca', 'SN2003dq', 'SN2008aw']:
            system = 'Swope'
        else:
            system = 'Landolt'

        with open(datafile, 'r') as ff:
            tsvin = csv.reader(ff, delimiter=' ', skipinitialspace=True)
            for row in tsvin:
                if not row[0]:
                    continue
                time = str(jd_to_mjd(Decimal(row[0])))
                add_photometry(catalog.events, name, time=time, band='V',
                               magnitude=row[1], e_magnitude=row[2],
                               system=system, source=source)
    catalog.journal_events()

    # stromlo
    stromlobands = ['B', 'V', 'R', 'I', 'VM', 'RM']
    file_path = os.path.join(
        PATH.REPO_EXTERNAL, 'J_A+A_415_863-1/photometry.csv')
    tsvin = list(csv.reader(open(file_path, 'r'), delimiter=','))
    for row in pbar(tsvin, current_task):
        name = row[0]
        name = catalog.add_event(name)
        source = catalog.events[name].add_source(bibcode='2004A&A...415..863G')
        catalog.events[name].add_quantity('alias', name, source)
        mjd = str(jd_to_mjd(Decimal(row[1])))
        for ri, ci in enumerate(range(2, len(row), 3)):
            if not row[ci]:
                continue
            band = stromlobands[ri]
            upperlimit = True if (not row[ci + 1] and row[ci + 2]) else False
            e_upper_magnitude = str(
                abs(Decimal(row[ci + 1]))) if row[ci + 1] else ''
            e_lower_magnitude = str(
                abs(Decimal(row[ci + 2]))) if row[ci + 2] else ''
            teles = 'MSSSO 1.3m' if band in ['VM', 'RM'] else 'CTIO'
            instr = 'MaCHO' if band in ['VM', 'RM'] else ''
            add_photometry(catalog.events, name, time=mjd, band=band,
                           magnitude=row[ci],
                           e_upper_magnitude=e_upper_magnitude,
                           e_lower_magnitude=e_lower_magnitude,
                           upperlimit=upperlimit, telescope=teles,
                           instrument=instr, source=source)
    catalog.journal_events()

    # 2015MNRAS.449..451W
    file_path = os.path.join(PATH.REPO_EXTERNAL, '2015MNRAS.449..451W.dat')
    data = list(csv.reader(open(file_path, 'r'), delimiter='\t',
                           quotechar='"', skipinitialspace=True))
    for rr, row in enumerate(pbar(data, current_task)):
        if rr == 0:
            continue
        namesplit = row[0].split('/')
        name = namesplit[-1]
        if name.startswith('SN'):
            name = name.replace(' ', '')
        name = catalog.add_event(name)
        source = catalog.events[name].add_source(bibcode='2015MNRAS.449..451W')
        catalog.events[name].add_quantity('alias', name, source)
        if len(namesplit) > 1:
            catalog.events[name].add_quantity('alias', namesplit[0], source)
        catalog.events[name].add_quantity('claimedtype', row[1], source)
        add_photometry(catalog.events, name, time=row[2], band=row[
                       4], magnitude=row[3], source=source)
    catalog.journal_events()

    # 2016MNRAS.459.1039T
    file_path = os.path.join(PATH.REPO_EXTERNAL, '2016MNRAS.459.1039T.tsv')
    data = list(csv.reader(open(file_path, 'r'), delimiter='\t',
                           quotechar='"', skipinitialspace=True))
    name = catalog.add_event('LSQ13zm')
    source = catalog.events[name].add_source(bibcode='2016MNRAS.459.1039T')
    catalog.events[name].add_quantity('alias', name, source)
    for rr, row in enumerate(pbar(data, current_task)):
        if row[0][0] == '#':
            bands = [xx.replace('(err)', '') for xx in row[3:-1]]
            continue
        mjd = row[1]
        mags = [re.sub(r'\([^)]*\)', '', xx) for xx in row[3:-1]]
        upps = [True if '>' in xx else '' for xx in mags]
        mags = [xx.replace('>', '') for xx in mags]
        errs = [xx[xx.find('(') + 1:xx.find(')')]
                if '(' in xx else '' for xx in row[3:-1]]
        for mi, mag in enumerate(mags):
            if not is_number(mag):
                continue
            add_photometry(
                catalog.events, name, time=mjd, band=bands[
                    mi], magnitude=mag, e_magnitude=errs[mi],
                instrument=row[-1], upperlimit=upps[mi], source=source)
    catalog.journal_events()

    # 2015ApJ...804...28G
    file_path = os.path.join(PATH.REPO_EXTERNAL, '2015ApJ...804...28G.tsv')
    data = list(csv.reader(open(file_path, 'r'), delimiter='\t',
                           quotechar='"', skipinitialspace=True))
    name = catalog.add_event('PS1-13arp')
    source = catalog.events[name].add_source(bibcode='2015ApJ...804...28G')
    catalog.events[name].add_quantity('alias', name, source)
    for rr, row in enumerate(pbar(data, current_task)):
        if rr == 0:
            continue
        mjd = row[1]
        mag = row[3]
        upp = True if '<' in mag else ''
        mag = mag.replace('<', '')
        err = row[4] if is_number(row[4]) else ''
        ins = row[5]
        add_photometry(
            catalog.events, name, time=mjd, band=row[
                0], magnitude=mag, e_magnitude=err,
            instrument=ins, upperlimit=upp, source=source)
    catalog.journal_events()

    # 2016ApJ...819...35A
    file_path = os.path.join(PATH.REPO_EXTERNAL, '2016ApJ...819...35A.tsv')
    data = list(csv.reader(open(file_path, 'r'), delimiter='\t',
                           quotechar='"', skipinitialspace=True))
    for rr, row in enumerate(pbar(data, current_task)):
        if row[0][0] == '#':
            continue
        name = catalog.add_event(row[0])
        source = catalog.events[name].add_source(bibcode='2016ApJ...819...35A')
        catalog.events[name].add_quantity('alias', name, source)
        catalog.events[name].add_quantity('ra', row[1], source)
        catalog.events[name].add_quantity('dec', row[2], source)
        catalog.events[name].add_quantity('redshift', row[3], source)
        disc_date = datetime.strptime(row[4], '%Y %b %d').isoformat()
        disc_date = disc_date.split('T')[0].replace('-', '/')
        catalog.events[name].add_quantity('discoverdate', disc_date, source)
    catalog.journal_events()

    # 2014ApJ...784..105W
    file_path = os.path.join(PATH.REPO_EXTERNAL, '2014ApJ...784..105W.tsv')
    data = list(csv.reader(open(file_path, 'r'), delimiter='\t',
                           quotechar='"', skipinitialspace=True))
    for rr, row in enumerate(pbar(data, current_task)):
        if row[0][0] == '#':
            continue
        name = catalog.add_event(row[0])
        source = catalog.events[name].add_source(bibcode='2014ApJ...784..105W')
        catalog.events[name].add_quantity('alias', name, source)
        mjd = row[1]
        band = row[2]
        mag = row[3]
        err = row[4]
        add_photometry(
            catalog.events, name, time=mjd, band=row[
                2], magnitude=mag, e_magnitude=err,
            instrument='WHIRC', telescope='WIYN 3.5 m', observatory='NOAO',
            system='WHIRC', source=source)
    catalog.journal_events()

    # 2012MNRAS.425.1007B
    file_path = os.path.join(PATH.REPO_EXTERNAL, '2012MNRAS.425.1007B.tsv')
    data = list(csv.reader(open(file_path, 'r'), delimiter='\t',
                           quotechar='"', skipinitialspace=True))
    for rr, row in enumerate(pbar(data, current_task)):
        if row[0][0] == '#':
            bands = row[2:]
            continue
        name = catalog.add_event(row[0])
        source = catalog.events[name].add_source(bibcode='2012MNRAS.425.1007B')
        catalog.events[name].add_quantity('alias', name, source)
        mjd = row[1]
        mags = [xx.split('±')[0].strip() for xx in row[2:]]
        errs = [xx.split('±')[1].strip()
                if '±' in xx else '' for xx in row[2:]]
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
                catalog.events, name, time=mjd, band=bands[
                    mi], magnitude=mag, e_magnitude=errs[mi],
                instrument=ins, telescope=tel, observatory=obs,
                system='Natural', source=source)

        catalog.journal_events()

    # 2014ApJ...783...28G
    file_path = os.path.join(PATH.REPO_EXTERNAL, 'apj490105t2_ascii.txt')
    with open(file_path, 'r') as f:
        data = list(csv.reader(f, delimiter='\t',
                               quotechar='"', skipinitialspace=True))
        for r, row in enumerate(pbar(data, current_task)):
            if row[0][0] == '#':
                continue
            (name,
             source) = catalog.new_event(row[0],
                                        bibcode='2014ApJ...783...28G')
            catalog.events[name].add_quantity('alias', row[1], source)
            catalog.events[name].add_quantity(
                'discoverdate', '20' + row[0][3:5], source)
            catalog.events[name].add_quantity('ra', row[2], source)
            catalog.events[name].add_quantity('dec', row[3], source)
            catalog.events[name].add_quantity(
                'redshift', row[13] if is_number(row[13]) else row[10], source)
    catalog.journal_events()

    # 2005ApJ...634.1190H
    file_path = os.path.join(PATH.REPO_EXTERNAL, '2005ApJ...634.1190H.tsv')
    with open(file_path, 'r') as f:
        data = list(csv.reader(f, delimiter='\t',
                               quotechar='"', skipinitialspace=True))
        for r, row in enumerate(pbar(data, current_task)):
            (name,
             source) = catalog.new_event('SNLS-' + row[0],
                                         bibcode='2005ApJ...634.1190H')
            catalog.events[name].add_quantity(
                'discoverdate', '20' + row[0][:2], source)
            catalog.events[name].add_quantity('ra', row[1], source)
            catalog.events[name].add_quantity('dec', row[2], source)
            catalog.events[name].add_quantity('redshift', row[5].replace(
                '?', ''), source, error=row[6], kind='host')
            catalog.events[name].add_quantity(
                'claimedtype', row[7].replace('SN', '').strip(':* '), source)
    catalog.journal_events()

    # 2014MNRAS.444.2133S
    file_path = os.path.join(PATH.REPO_EXTERNAL, '2014MNRAS.444.2133S.tsv')
    with open(file_path, 'r') as f:
        data = list(csv.reader(f, delimiter='\t',
                               quotechar='"', skipinitialspace=True))
        for r, row in enumerate(pbar(data, current_task)):
            if row[0][0] == '#':
                continue
            name = row[0]
            if is_number(name[:4]):
                name = 'SN' + name
            events, name, source = Events.new_event(
                tasks, args, events, name, log, bibcode='2014MNRAS.444.2133S')
            catalog.events[name].add_quantity('ra', row[1], source)
            catalog.events[name].add_quantity('dec', row[2], source)
            catalog.events[name].add_quantity('redshift', row[3], source, kind='host')
    catalog.journal_events()

    # 2009MNRAS.398.1041B
    file_path = os.path.join(PATH.REPO_EXTERNAL, '2009MNRAS.398.1041B.tsv')
    with open(file_path, 'r') as f:
        data = list(csv.reader(f, delimiter='\t',
                               quotechar='"', skipinitialspace=True))
        for r, row in enumerate(pbar(data, current_task)):
            if row[0][0] == '#':
                bands = row[2:-1]
                continue
            (catalog.events,
             name,
             source) = catalog.new_event('SN2008S',
                                        bibcode='2009MNRAS.398.1041B')
            mjd = str(jd_to_mjd(Decimal(row[0])))
            mags = [x.split('±')[0].strip() for x in row[2:]]
            upps = [('<' in x.split('±')[0]) for x in row[2:]]
            errs = [x.split('±')[1].strip()
                    if '±' in x else '' for x in row[2:]]

            instrument = row[-1]

            for mi, mag in enumerate(mags):
                if not is_number(mag):
                    continue
                add_photometry(catalog.events, name, time=mjd, band=bands[mi],
                               magnitude=mag, e_magnitude=errs[mi],
                               instrument=instrument, source=source)
    catalog.journal_events()

    # 2010arXiv1007.0011P
    file_path = os.path.join(PATH.REPO_EXTERNAL, '2010arXiv1007.0011P.tsv')
    with open(file_path, 'r') as f:
        data = list(csv.reader(f, delimiter='\t',
                               quotechar='"', skipinitialspace=True))
        for r, row in enumerate(pbar(data, current_task)):
            if row[0][0] == '#':
                bands = row[1:]
                continue
            (catalog.events,
             name,
             source) = catalog.new_event('SN2008S',
                                        bibcode='2010arXiv1007.0011P')
            mjd = row[0]
            mags = [x.split('±')[0].strip() for x in row[1:]]
            errs = [x.split('±')[1].strip()
                    if '±' in x else '' for x in row[1:]]

            for mi, mag in enumerate(mags):
                if not is_number(mag):
                    continue
                add_photometry(catalog.events, name, time=mjd, band=bands[mi],
                               magnitude=mag, e_magnitude=errs[mi],
                               instrument='LBT', source=source)
    catalog.journal_events()

    # 2000ApJ...533..320G
    file_path = os.path.join(PATH.REPO_EXTERNAL, '2000ApJ...533..320G.tsv')
    with open(file_path, 'r') as f:
        data = list(csv.reader(f, delimiter='\t',
                               quotechar='"', skipinitialspace=True))
        (catalog.events,
         name,
         source) = catalog.new_event('SN1997cy',
                                    bibcode='2000ApJ...533..320G')
        for r, row in enumerate(pbar(data, current_task)):
            if row[0][0] == '#':
                bands = row[1:-1]
                continue
            mjd = str(jd_to_mjd(Decimal(row[0])))
            mags = row[1:len(bands)]
            for mi, mag in enumerate(mags):
                if not is_number(mag):
                    continue
                add_photometry(catalog.events, name, time=mjd, band=bands[mi],
                               magnitude=mag,
                               observatory='Mount Stromlo', telescope='MSSSO',
                               source=source, kcorrected=True)

    catalog.journal_events()
    return
