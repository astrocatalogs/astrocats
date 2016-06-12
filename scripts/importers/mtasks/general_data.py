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


def do_crts(events, args, tasks):
    crtsnameerrors = ['2011ax']

    folders = ['catalina', 'MLS', 'SSS']
    for fold in pbar(folders, current_task):
        html = load_cached_url(args, 'http://nesssi.cacr.caltech.edu/' + fold + '/AllSN.html',
                               os.path.join(PATH.REPO_EXTERNAL, 'CRTS', fold + '.html'))
        if not html:
            continue
        bs = BeautifulSoup(html, 'html5lib')
        trs = bs.findAll('tr')
        for tr in pbar(trs, current_task):
            tds = tr.findAll('td')
            if not tds:
                continue
            refs = []
            aliases = []
            ttype = ''
            ctype = ''
            for tdi, td in enumerate(tds):
                if tdi == 0:
                    crtsname = td.contents[0].text.strip()
                elif tdi == 1:
                    ra = td.contents[0]
                elif tdi == 2:
                    dec = td.contents[0]
                elif tdi == 11:
                    lclink = td.find('a')['onclick']
                    lclink = lclink.split("'")[1]
                elif tdi == 13:
                    aliases = re.sub('[()]', '', re.sub('<[^<]+?>', '', td.contents[-1].strip()))
                    aliases = [x.strip('; ') for x in list(filter(None, aliases.split(' ')))]

            name = ''
            hostmag = ''
            hostupper = False
            validaliases = []
            for ai, alias in enumerate(aliases):
                if alias in ['SN', 'SDSS']:
                    continue
                if alias in crtsnameerrors:
                    continue
                if alias == 'mag':
                    if ai < len(aliases) - 1:
                        ind = ai+1
                        if aliases[ai+1] in ['SDSS']:
                            ind = ai+2
                        elif aliases[ai+1] in ['gal', 'obj', 'object', 'source']:
                            ind = ai-1
                        if '>' in aliases[ind]:
                            hostupper = True
                        hostmag = aliases[ind].strip('>~').replace(',', '.')
                    continue
                if is_number(alias[:4]) and alias[:2] == '20' and len(alias) > 4:
                    name = 'SN' + alias
                # lalias = alias.lower()
                if ((('asassn' in alias and len(alias) > 6) or
                     ('ptf' in alias and len(alias) > 3) or
                     ('ps1' in alias and len(alias) > 3) or 'snhunt' in alias or
                     ('mls' in alias and len(alias) > 3) or 'gaia' in alias or
                     ('lsq' in alias and len(alias) > 3))):
                    alias = alias.replace('SNHunt', 'SNhunt')
                    validaliases.append(alias)

            if not name:
                name = crtsname
            name = add_event(tasks, args, events, name)
            source = add_source(
                events, name, refname='Catalina Sky Survey', bibcode='2009ApJ...696..870D',
                url='http://nesssi.cacr.caltech.edu/catalina/AllSN.html')
            add_quantity(events, name, 'alias', name, source)
            for alias in validaliases:
                add_quantity(events, name, 'alias', alias, source)
            add_quantity(events, name, 'ra', ra, source, unit='floatdegrees')
            add_quantity(events, name, 'dec', dec, source, unit='floatdegrees')

            if hostmag:
                # 1.0 magnitude error based on Drake 2009 assertion that SN are only considered
                #    real if they are 2 mags brighter than host.
                add_photometry(
                    events, name, band='C', magnitude=hostmag, e_magnitude=1.0, source=source,
                    host=True, telescope='Catalina Schmidt', upperlimit=hostupper)

            fname2 = (PATH.REPO_EXTERNAL + '/' + fold + '/' +
                      lclink.split('.')[-2].rstrip('p').split('/')[-1] + '.html')
            if ((not args.fullrefresh and archived_task(tasks, args, 'crts') and
                 os.path.isfile(fname2))):
                with open(fname2, 'r') as f:
                    html2 = f.read()
            else:
                with open(fname2, 'w') as f:
                    response2 = urllib.request.urlopen(lclink)
                    html2 = response2.read().decode('utf-8')
                    f.write(html2)

            lines = html2.splitlines()
            teles = 'Catalina Schmidt'
            for line in lines:
                if 'javascript:showx' in line:
                    mjdstr = re.search("showx\('(.*?)'\)", line).group(1).split('(')[0].strip()
                    if not is_number(mjdstr):
                        continue
                    mjd = str(Decimal(mjdstr) + Decimal(53249.0))
                else:
                    continue
                if 'javascript:showy' in line:
                    mag = re.search("showy\('(.*?)'\)", line).group(1)
                if 'javascript:showz' in line:
                    err = re.search("showz\('(.*?)'\)", line).group(1)
                e_mag = err if float(err) > 0.0 else ''
                upl = (float(err) == 0.0)
                add_photometry(
                    events, name, time=mjd, band='C', magnitude=mag, source=source,
                    includeshost=True, telescope=teles, e_magnitude=e_mag, upperlimit=upl)
            if args.update:
                events = journal_events(tasks, args, events)

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


def do_fermi(events, args, tasks):
    with open(os.path.join(PATH.REPO_EXTERNAL, '1SC_catalog_v01.asc'), 'r') as f:
        tsvin = csv.reader(f, delimiter=',')
        for ri, row in enumerate(pbar(tsvin, current_task)):
            if row[0].startswith('#'):
                if len(row) > 1 and 'UPPER_LIMITS' in row[1]:
                    break
                continue
            if 'Classified' not in row[1]:
                continue
            name = row[0].replace('SNR', 'G')
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2016ApJS..224....8A')
            add_quantity(events, name, 'alias', name, source)
            add_quantity(events, name, 'alias', row[0].replace('SNR', 'MWSNR'), source)
            add_quantity(events, name, 'ra', row[2], source, unit='floatdegrees')
            add_quantity(events, name, 'dec', row[3], source, unit='floatdegrees')
    events = journal_events(tasks, args, events)
    return events


def do_gaia(events, args, tasks):
    fname = os.path.join(PATH.REPO_EXTERNAL, 'GAIA/alerts.csv')
    csvtxt = load_cached_url(args, 'http://gsaweb.ast.cam.ac.uk/alerts/alerts.csv', fname)
    if not csvtxt:
        continue
    tsvin = csv.reader(csvtxt.splitlines(), delimiter=',', skipinitialspace=True)
    reference = 'Gaia Photometric Science Alerts'
    refurl = 'http://gsaweb.ast.cam.ac.uk/alerts/alertsindex'
    for ri, row in enumerate(pbar(tsvin, current_task)):
        if ri == 0 or not row:
            continue
        name = add_event(tasks, args, events, row[0])
        source = add_source(events, name, refname=reference, url=refurl)
        add_quantity(events, name, 'alias', name, source)
        year = '20' + re.findall(r'\d+', row[0])[0]
        add_quantity(events, name, 'discoverdate', year, source)
        add_quantity(events, name, 'ra', row[2], source, unit='floatdegrees')
        add_quantity(events, name, 'dec', row[3], source, unit='floatdegrees')
        if row[7] and row[7] != 'unknown':
            type = row[7].replace('SNe', '').replace('SN', '').strip()
            add_quantity(events, name, 'claimedtype', type, source)
        elif any([x in row[9].upper() for x in ['SN CANDIATE', 'CANDIDATE SN', 'HOSTLESS SN']]):
            add_quantity(events, name, 'claimedtype', 'Candidate', source)

        if 'aka' in row[9].replace('gakaxy', 'galaxy').lower() and 'AKARI' not in row[9]:
            commentsplit = (row[9].replace('_', ' ').replace('MLS ', 'MLS').replace('CSS ', 'CSS').
                            replace('SN iPTF', 'iPTF').replace('SN ', 'SN').replace('AT ', 'AT'))
            commentsplit = commentsplit.split()
            for csi, cs in enumerate(commentsplit):
                if 'aka' in cs.lower() and csi < len(commentsplit) - 1:
                    alias = commentsplit[csi+1].strip('(),:.').replace('PSNJ', 'PSN J')
                    if alias[:6] == 'ASASSN' and alias[6] != '-':
                        alias = 'ASASSN-' + alias[6:]
                    add_quantity(events, name, 'alias', alias, source)
                    break

        fname = os.path.join(PATH.REPO_EXTERNAL, 'GAIA/') + row[0] + '.csv'
        if not args.fullrefresh and archived_task(tasks, args, 'gaia') and os.path.isfile(fname):
            with open(fname, 'r') as f:
                csvtxt = f.read()
        else:
            response = urllib.request.urlopen('http://gsaweb.ast.cam.ac.uk/alerts/alert/' +
                                              row[0] + '/lightcurve.csv')
            with open(fname, 'w') as f:
                csvtxt = response.read().decode('utf-8')
                f.write(csvtxt)

        tsvin2 = csv.reader(csvtxt.splitlines())
        for ri2, row2 in enumerate(tsvin2):
            if ri2 <= 1 or not row2:
                continue
            mjd = str(jd_to_mjd(Decimal(row2[1].strip())))
            magnitude = row2[2].strip()
            if magnitude == 'null':
                continue
            e_mag = 0.
            telescope = 'GAIA'
            band = 'G'
            add_photometry(
                events, name, time=mjd, telescope=telescope, band=band, magnitude=magnitude,
                e_magnitude=e_mag, source=source)
        if args.update:
            events = journal_events(tasks, args, events)
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


def do_pessto(events, args, tasks):
    with open(os.path.join(PATH.REPO_EXTERNAL, 'PESSTO_MPHOT.csv'), 'r') as f:
        tsvin = csv.reader(f, delimiter=',')
        for ri, row in enumerate(tsvin):
            if ri == 0:
                bands = [x.split('_')[0] for x in row[3::2]]
                systems = [x.split('_')[1].capitalize().replace('Ab', 'AB') for x in row[3::2]]
                continue
            name = row[1]
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2015A&A...579A..40S')
            add_quantity(events, name, 'alias', name, source)
            for hi, ci in enumerate(range(3, len(row)-1, 2)):
                if not row[ci]:
                    continue
                add_photometry(
                    events, name, time=row[2], magnitude=row[ci], e_magnitude=row[ci+1],
                    band=bands[hi], system=systems[hi], telescope='Swift' if systems[hi] == 'Swift' else '',
                    source=source)

    events = journal_events(tasks, args, events)
    return events


def do_scp(events, args, tasks):
    with open(os.path.join(PATH.REPO_EXTERNAL, 'SCP09.csv'), 'r') as f:
        tsvin = csv.reader(f, delimiter=',')
        for ri, row in enumerate(pbar(tsvin, current_task)):
            if ri == 0:
                continue
            name = row[0].replace('SCP', 'SCP-')
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, refname='Supernova Cosmology Project',
                                url='http://supernova.lbl.gov/2009ClusterSurvey/')
            add_quantity(events, name, 'alias', name, source)
            if row[1]:
                add_quantity(events, name, 'alias', row[1], source)
            if row[2]:
                kind = 'spectroscopic' if row[3] == 'sn' else 'host'
                add_quantity(events, name, 'redshift', row[2], source, kind=kind)
            if row[4]:
                add_quantity(events, name, 'redshift', row[2], source, kind='cluster')
            if row[6]:
                claimedtype = row[6].replace('SN ', '')
                kind = ('spectroscopic/light curve' if 'a' in row[7] and 'c' in row[7] else
                        'spectroscopic' if 'a' in row[7] else
                        'light curve' if 'c' in row[7]
                        else '')
                if claimedtype != '?':
                    add_quantity(events, name, 'claimedtype', claimedtype, source, kind=kind)
    events = journal_events(tasks, args, events)
    return events


def do_sdss(events, args, tasks):
    with open(os.path.join(PATH.REPO_EXTERNAL, 'SDSS/2010ApJ...708..661D.txt'), 'r') as f:
        bibcodes2010 = f.read().split('\n')
    sdssbands = ['u', 'g', 'r', 'i', 'z']
    file_names = glob(os.path.join(PATH.REPO_EXTERNAL, 'SDSS/*.sum'))
    for fname in pbar_strings(file_names, desc=current_task):
        f = open(fname, 'r')
        tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)

        basename = os.path.basename(fname)
        if basename in bibcodes2010:
            bibcode = '2010ApJ...708..661D'
        else:
            bibcode = '2008AJ....136.2306H'

        for r, row in enumerate(tsvin):
            if r == 0:
                if row[5] == 'RA:':
                    name = 'SDSS-II ' + row[3]
                else:
                    name = 'SN' + row[5]
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode=bibcode)
                add_quantity(events, name, 'alias', name, source)
                add_quantity(events, name, 'alias', 'SDSS-II ' + row[3], source)

                if row[5] != 'RA:':
                    year = re.findall(r'\d+', name)[0]
                    add_quantity(events, name, 'discoverdate', year, source)

                add_quantity(events, name, 'ra', row[-4], source, unit='floatdegrees')
                add_quantity(events, name, 'dec', row[-2], source, unit='floatdegrees')
            if r == 1:
                error = row[4] if float(row[4]) >= 0.0 else ''
                add_quantity(events, name, 'redshift', row[2], source, error=error,
                             kind='heliocentric')
            if r >= 19:
                # Skip bad measurements
                if int(row[0]) > 1024:
                    continue

                mjd = row[1]
                band = sdssbands[int(row[2])]
                magnitude = row[3]
                e_mag = row[4]
                telescope = 'SDSS'
                add_photometry(
                    events, name, time=mjd, telescope=telescope, band=band, magnitude=magnitude,
                    e_magnitude=e_mag, source=source, system='SDSS')
        f.close()
    events = journal_events(tasks, args, events)
    return events


def do_snls(events, args, tasks):
    with open(os.path.join(PATH.REPO_EXTERNAL, 'SNLS-ugriz.dat'), 'r') as f:
        data = csv.reader(f, delimiter=' ', quotechar='"', skipinitialspace=True)
        for row in data:
            flux = row[3]
            err = row[4]
            # Being extra strict here with the flux constraint, see note below.
            if float(flux) < 3.0*float(err):
                continue
            name = 'SNLS-' + row[0]
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2010A&A...523A...7G')
            add_quantity(events, name, 'alias', name, source)
            band = row[1]
            mjd = row[2]
            sig = get_sig_digits(flux.split('E')[0])+1
            # Conversion comes from SNLS-Readme
            # NOTE: Datafiles avail for download suggest diff zeropoints than 30, need to inquire.
            magnitude = pretty_num(30.0-2.5*log10(float(flux)), sig=sig)
            e_mag = pretty_num(2.5*log10(1.0 + float(err)/float(flux)), sig=sig)
            # e_mag = pretty_num(2.5*(log10(float(flux) + float(err)) - log10(float(flux))), sig=sig)
            add_photometry(
                events, name, time=mjd, band=band, magnitude=magnitude, e_magnitude=e_mag, counts=flux,
                e_counts=err, source=source)

    events = journal_events(tasks, args, events)
    return events


def do_ucb(events, args, tasks):
    sec_ref = 'UCB Filippenko Group\'s Supernova Database (SNDB)'
    sec_refurl = 'http://heracles.astro.berkeley.edu/sndb/info'
    sec_refbib = '2012MNRAS.425.1789S'

    jsontxt = load_cached_url(
        args, 'http://heracles.astro.berkeley.edu/sndb/download?id=allpubphot',
        os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'UCB/allpub.json'))
    if not jsontxt:
        continue

    photom = json.loads(jsontxt)
    photom = sorted(photom, key=lambda k: k['ObjName'])
    for phot in pbar(photom, desc=current_task):
        name = phot['ObjName']
        name = add_event(tasks, args, events, name)

        sec_source = add_source(events, name, refname=sec_ref, url=sec_refurl, bibcode=sec_refbib,
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
            raise(ValueError('Filename not found for SNDB phot!'))
        if not phot['PhotID']:
            raise(ValueError('ID not found for SNDB phot!'))

        filepath = os.path.join(PATH.REPO_EXTERNAL, 'SNDB/') + filename
        if archived_task(tasks, args, 'ucb') and os.path.isfile(filepath):
            with open(filepath, 'r') as f:
                phottxt = f.read()
        else:
            session = requests.Session()
            response = session.get('http://heracles.astro.berkeley.edu/sndb/download?id=dp:' +
                                   str(phot['PhotID']))
            phottxt = response.text
            with open(filepath, 'w') as f:
                f.write(phottxt)

        tsvin = csv.reader(phottxt.splitlines(), delimiter=' ', skipinitialspace=True)

        for r, row in enumerate(tsvin):
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
