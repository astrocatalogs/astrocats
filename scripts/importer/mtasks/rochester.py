"""General data import tasks.
"""
from astropy.time import Time as astrotime
from bs4 import BeautifulSoup
import csv
from math import floor
import os
import re
from string import ascii_letters

from scripts import PATH
from .. import Events
from .. funcs import add_photometry, get_aliases, load_cached_url, make_date_string, uniq_cdl
from ... utils import is_number, pbar


def do_rochester(events, stubs, args, tasks, task_obj, log):
    rochestermirrors = ['http://www.rochesterastronomy.org/', 'http://www.supernova.thistlethwaites.com/']
    rochesterpaths = ['snimages/snredshiftall.html', 'sn2016/snredshift.html', 'snimages/snredboneyard.html']
    rochesterupdate = [False, True, True]
    current_task = task_obj.current_task(args)

    for pp, path in enumerate(pbar(rochesterpaths, current_task)):
        if args.update and not rochesterupdate[pp]:
            continue

        filepath = os.path.join(PATH.REPO_EXTERNAL, 'rochester/') + os.path.basename(path)
        for mirror in rochestermirrors:
            html = load_cached_url(mirror + path, filepath, failhard = (mirror != rochestermirrors[-1]))
            if html:
                break

        if not html:
            continue

        soup = BeautifulSoup(html, 'html5lib')
        rows = soup.findAll('tr')
        sec_ref = 'Latest Supernovae'
        sec_refurl = 'http://www.rochesterastronomy.org/snimages/snredshiftall.html'
        for rr, row in enumerate(pbar(rows, current_task)):
            if rr == 0:
                continue
            cols = row.findAll('td')
            if not len(cols):
                continue

            name = ''
            if cols[14].contents:
                aka = str(cols[14].contents[0]).strip()
                if is_number(aka.strip('?')):
                    aka = 'SN' + aka.strip('?') + 'A'
                    oldname = aka
                    events, name = Events.add_event(tasks, args, events, aka, log)
                elif len(aka) = 4 and is_number(aka[:4]):
                    aka = 'SN' + aka
                    oldname = aka
                    events, name = Events.add_event(tasks, args, events, aka, log)

            ra = str(cols[3].contents[0]).strip()
            dec = str(cols[4].contents[0]).strip()

            sn = re.sub('<[^<]+?>', '', str(cols[0].contents[0])).strip()
            if is_number(sn.strip('?')):
                sn = 'SN' + sn.strip('?') + 'A'
            elif len(sn) == 4 and is_number(sn[:4]):
                sn = 'SN' + sn
            if not name:
                if not sn:
                    continue
                if sn[:8] == 'MASTER J':
                    sn = sn.replace('MASTER J', 'MASTER OT J').replace('SNHunt', 'SNhunt')
                if 'POSSIBLE' in sn.upper() and ra and dec:
                    sn = 'PSN J' + ra.replace(':', '').replace('.', '')
                    sn += dec.replace(':', '').replace('.', '')
                oldname = sn
                events, name = Events.add_event(tasks, args, events, sn, log)

            reference = cols[12].findAll('a')[0].contents[0].strip()
            refurl = cols[12].findAll('a')[0]['href'].strip()
            source = events[name].add_source(srcname=reference, url=refurl)
            sec_source = events[name].add_source(srcname=sec_ref, url=sec_refurl, secondary=True)
            sources = uniq_cdl(list(filter(None, [source, sec_source])))
            events[name].add_quantity('alias', oldname, sources)
            events[name].add_quantity('alias', sn, sources)

            if cols[14].contents:
                if aka == 'SNR G1.9+0.3':
                    aka = 'G001.9+00.3'
                if aka[:4] == 'PS1 ':
                    aka = 'PS1-' + aka[4:]
                if aka[:8] == 'MASTER J':
                    aka = aka.replace('MASTER J', 'MASTER OT J').replace('SNHunt', 'SNhunt')
                if 'POSSIBLE' in aka.upper() and ra and dec:
                    aka = 'PSN J' + ra.replace(':', '').replace('.', '')
                    aka += dec.replace(':', '').replace('.', '')
                events[name].add_quantity('alias', aka, sources)

            if str(cols[1].contents[0]).strip() != 'unk':
                type = str(cols[1].contents[0]).strip(' :,')
                events[name].add_quantity('claimedtype', type, sources)
            if str(cols[2].contents[0]).strip() != 'anonymous':
                events[name].add_quantity('host', str(cols[2].contents[0]).strip(), sources)
            events[name].add_quantity('ra', ra, sources)
            events[name].add_quantity('dec', dec, sources)
            if str(cols[6].contents[0]).strip() not in ['2440587', '2440587.292']:
                astrot = astrotime(float(str(cols[6].contents[0]).strip()), format='jd').datetime
                ddate = make_date_string(astrot.year, astrot.month, astrot.day)
                events[name].add_quantity('discoverdate', ddate, sources)
            if str(cols[7].contents[0]).strip() not in ['2440587', '2440587.292']:
                astrot = astrotime(float(str(cols[7].contents[0]).strip()), format='jd')
                if ((float(str(cols[8].contents[0]).strip()) <= 90.0 and
                     not any('GRB' in xx for xx in events[name].get_aliases()))):
                    mag = str(cols[8].contents[0]).strip()
                    add_photometry(
                        events, name, time=str(astrot.mjd), magnitude=mag, source=sources)
            if cols[11].contents[0] != 'n/a':
                events[name].add_quantity('redshift', str(cols[11].contents[0]).strip(), sources)
            events[name].add_quantity('discoverer', str(cols[13].contents[0]).strip(), sources)
            if args.update:
                events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    if not args.update:
        vsnetfiles = ['latestsne.dat']
        for vsnetfile in vsnetfiles:
            file_name = os.path.join(PATH.REPO_EXTERNAL, "" + vsnetfile)
            with open(file_name, 'r', encoding='latin1') as csv_file:
                tsvin = csv.reader(csv_file, delimiter=' ', skipinitialspace=True)
            for rr, row in enumerate(tsvin):
                if not row or row[0][:4] in ['http', 'www.'] or len(row) < 3:
                    continue
                name = row[0].strip()
                if name[:4].isdigit():
                    name = 'SN' + name
                if name.startswith('PSNJ'):
                    name = 'PSN J' + name[4:]
                if name.startswith('MASTEROTJ'):
                    name = name.replace('MASTEROTJ', 'MASTER OT J')
                events, name = Events.add_event(tasks, args, events, name, log)
                sec_source = events[name].add_source(
                    refname=sec_ref, url=sec_refurl, secondary=True)
                events[name].add_quantity('alias', name, sec_source)

                if not is_number(row[1]):
                    continue
                year = row[1][:4]
                month = row[1][4:6]
                day = row[1][6:]
                if '.' not in day:
                    day = day[:2] + '.' + day[2:]
                mjd = astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd
                mjd += float(day) - floor(float(day))
                magnitude = row[2].rstrip(ascii_letters)
                if not is_number(magnitude):
                    continue
                if magnitude.isdigit():
                    if int(magnitude) > 100:
                        magnitude = magnitude[:2] + '.' + magnitude[2:]

                if float(str(cols[8].contents[0]).strip()) >= 90.0:
                    continue

                if len(row) >= 4:
                    if is_number(row[3]):
                        e_magnitude = row[3]
                        refind = 4
                    else:
                        e_magnitude = ''
                        refind = 3

                    if refind >= len(row):
                        sources = sec_source
                    else:
                        reference = ' '.join(row[refind:])
                        source = events[name].add_source(srcname=reference)
                        events[name].add_quantity('alias', name, sec_source)
                        sources = uniq_cdl([source, sec_source])
                else:
                    sources = sec_source

                band = row[2].lstrip('1234567890.')

                add_photometry(
                    events, name, time=mjd, band=band, magnitude=magnitude,
                    e_magnitude=e_magnitude, source=sources)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events
