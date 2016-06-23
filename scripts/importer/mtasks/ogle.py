"""Import data from OGLE.
"""
from bs4 import BeautifulSoup, Tag, NavigableString
import os
import re
import urllib

from scripts import PATH
from ... utils import is_number, pbar, Decimal
from .. funcs import add_event, add_photometry, add_source, add_quantity, \
    jd_to_mjd, journal_events, load_cached_url, uniq_cdl


def do_ogle(events, stubs, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    basenames = ['transients', 'transients/2014b', 'transients/2014',
                 'transients/2013', 'transients/2012']
    oglenames = []
    ogleupdate = [True, False, False, False, False]
    for b, bn in enumerate(pbar(basenames, current_task)):
        if args.update and not ogleupdate[b]:
            continue

        filepath = os.path.join(PATH.REPO_EXTERNAL, 'OGLE-')
        filepath += bn.replace('/', '-') + '-transients.html'
        htmltxt = load_cached_url(args, current_task, 'http://ogle.astrouw.edu.pl/ogle4/' + bn +
                                  '/transients.html', filepath)
        if not htmltxt:
            continue

        soup = BeautifulSoup(htmltxt, 'html5lib')
        links = soup.findAll('a')
        breaks = soup.findAll('br')
        datalinks = []
        datafnames = []
        for a in links:
            if a.has_attr('href'):
                if '.dat' in a['href']:
                    datalinks.append('http://ogle.astrouw.edu.pl/ogle4/' + bn + '/' + a['href'])
                    datafnames.append(bn.replace('/', '-') + '-' + a['href'].replace('/', '-'))

        ec = -1
        reference = 'OGLE-IV Transient Detection System'
        refurl = 'http://ogle.astrouw.edu.pl/ogle4/transients/transients.html'
        for br in pbar(breaks, current_task):
            sibling = br.nextSibling
            if 'Ra,Dec=' in sibling:
                line = sibling.replace('\n', '').split('Ra,Dec=')
                name = line[0].strip()
                ec += 1

                if 'NOVA' in name or 'dupl' in name:
                    continue

                if name in oglenames:
                    continue
                oglenames.append(name)

                events, name = add_event(tasks, args, events, name, log)

                mySibling = sibling.nextSibling
                atelref = ''
                claimedtype = ''
                while 'Ra,Dec=' not in mySibling:
                    if isinstance(mySibling, NavigableString):
                        if 'Phot.class=' in str(mySibling):
                            claimedtype = re.sub(r'\([^)]*\)', '', str(mySibling).split('=')[-1])
                            claimedtype = claimedtype.replace('SN', '').strip()
                    if isinstance(mySibling, Tag):
                        atela = mySibling
                        if atela and atela.has_attr('href') and 'astronomerstelegram' in atela['href']:
                            atelref = atela.contents[0].strip()
                            atelurl = atela['href']
                    mySibling = mySibling.nextSibling
                    if mySibling is None:
                        break

                # nextSibling = sibling.nextSibling
                # if ((isinstance(nextSibling, Tag) and nextSibling.has_attr('alt') and
                #      nextSibling.contents[0].strip() != 'NED')):
                #     radec = nextSibling.contents[0].strip().split()
                # else:
                #     radec = line[-1].split()
                # ra = radec[0]
                # dec = radec[1]

                fname = os.path.join(PATH.REPO_EXTERNAL, 'OGLE/') + datafnames[ec]
                if task_obj.load_archive(args) and os.path.isfile(fname):
                    with open(fname, 'r') as f:
                        csvtxt = f.read()
                else:
                    response = urllib.request.urlopen(datalinks[ec])
                    with open(fname, 'w') as f:
                        csvtxt = response.read().decode('utf-8')
                        f.write(csvtxt)

                lcdat = csvtxt.splitlines()
                sources = [events[name].add_source(srcname=reference, url=refurl)]
                events[name].add_quantity('alias', name, sources[0])
                if atelref and atelref != 'ATel#----':
                    sources.append(events[name].add_source(srcname=atelref, url=atelurl))
                sources = uniq_cdl(sources)

                if name.startswith('OGLE'):
                    if name[4] == '-':
                        if is_number(name[5:9]):
                            events[name].add_quantity('discoverdate', name[5:9], sources)
                    else:
                        if is_number(name[4:6]):
                            events[name].add_quantity('discoverdate', '20' + name[4:6], sources)

                # RA and Dec from OGLE pages currently not reliable
                # events[name].add_quantity('ra', ra, sources)
                # events[name].add_quantity('dec', dec, sources)
                if claimedtype and claimedtype != '-':
                    events[name].add_quantity('claimedtype', claimedtype, sources)
                elif 'SN' not in name and 'claimedtype' not in events[name]:
                    events[name].add_quantity('claimedtype', 'Candidate', sources)
                for row in lcdat:
                    row = row.split()
                    mjd = str(jd_to_mjd(Decimal(row[0])))
                    magnitude = row[1]
                    if float(magnitude) > 90.0:
                        continue
                    e_mag = row[2]
                    upperlimit = False
                    if e_mag == '-1' or float(e_mag) > 10.0:
                        e_mag = ''
                        upperlimit = True
                    add_photometry(
                        events, name, time=mjd, band='I', magnitude=magnitude, e_magnitude=e_mag,
                        system='Vega', source=sources, upperlimit=upperlimit)
                if args.update:
                    events, stubs = journal_events(tasks, args, events, stubs, log)

        events, stubs = journal_events(tasks, args, events, stubs, log)
    return events
