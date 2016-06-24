"""Import data from Pan-STARRS.
"""
from astropy.time import Time as astrotime
from bs4 import BeautifulSoup
from glob import glob
import json
import os
import urllib
import warnings

from scripts import PATH
from .. import Events
from .. funcs import add_photometry, load_cached_url, make_date_string, uniq_cdl
from ... utils import is_number, pbar


def do_ps_mds(events, stubs, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    with open(os.path.join(PATH.REPO_EXTERNAL, 'MDS/apj506838t1_mrt.txt')) as f:
        for ri, row in enumerate(pbar(f.read().splitlines(), current_task)):
            if ri < 35:
                continue
            cols = [x.strip() for x in row.split(',')]
            events, name = Events.add_event(tasks, args, events, cols[0], log)
            source = events[name].add_source(bibcode='2015ApJ...799..208S')
            events[name].add_quantity('alias', name, source)
            events[name].add_quantity('ra', cols[2], source)
            events[name].add_quantity('dec', cols[3], source)
            astrot = astrotime(float(cols[4]), format='mjd').datetime
            ddate = make_date_string(astrot.year, astrot.month, astrot.day)
            events[name].add_quantity('discoverdate', ddate, source)
            events[name].add_quantity('redshift', cols[5], source, kind='spectroscopic')
            events[name].add_quantity('claimedtype', 'II P', source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_ps_threepi(events, stubs, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    teles = 'Pan-STARRS1'
    fname = os.path.join(PATH.REPO_EXTERNAL, '3pi/page00.html')
    ps_url = 'http://psweb.mp.qub.ac.uk/ps1threepi/psdb/public/?page=1&sort=followup_flag_date'
    html = load_cached_url(args, current_task, ps_url, fname, write=False)
    if not html:
        return events

    bs = BeautifulSoup(html, 'html5lib')
    div = bs.find('div', {'class': 'pagination'})
    offline = False
    if not div:
        offline = True
    else:
        links = div.findAll('a')
        if not links:
            offline = True

    if offline:
        if args.update:
            return events
        warnings.warn('Pan-STARRS 3pi offline, using local files only.')
        with open(fname, 'r') as f:
            html = f.read()
        bs = BeautifulSoup(html, 'html5lib')
        div = bs.find('div', {'class': 'pagination'})
        links = div.findAll('a')
    else:
        with open(fname, 'w') as f:
            f.write(html)

    numpages = int(links[-2].contents[0])
    oldnumpages = len(glob(os.path.join(PATH.REPO_EXTERNAL, '3pi/page*')))
    for page in pbar(range(1, numpages), current_task):
        fname = os.path.join(PATH.REPO_EXTERNAL, '3pi/page') + str(page).zfill(2) + '.html'
        if ((task_obj.load_archive(args) and
             os.path.isfile(fname) and page < oldnumpages)):
            with open(fname, 'r') as f:
                html = f.read()
        elif not offline:
            use_url = ('http://psweb.mp.qub.ac.uk/ps1threepi/psdb/public/?page=' + str(page) +
                       '&sort=followup_flag_date')
            response = urllib.request.urlopen(use_url)
            with open(fname, 'w') as f:
                html = response.read().decode('utf-8')
                f.write(html)
        else:
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
                    psname = td.contents[0]
                    pslink = psname['href']
                    psname = psname.text
                elif tdi == 1:
                    ra = td.contents[0]
                elif tdi == 2:
                    dec = td.contents[0]
                elif tdi == 3:
                    ttype = td.contents[0]
                    if ttype != 'sn' and ttype != 'orphan':
                        break
                elif tdi == 5:
                    if not td.contents:
                        continue
                    ctype = td.contents[0]
                    if ctype == 'Observed':
                        ctype = ''
                elif tdi == 16:
                    if td.contents:
                        crossrefs = td.findAll('a')
                        for cref in crossrefs:
                            if 'atel' in cref.contents[0].lower():
                                refs.append([cref.contents[0], cref['href']])
                            elif is_number(cref.contents[0][:4]):
                                continue
                            else:
                                aliases.append(cref.contents[0])

            if ttype != 'sn' and ttype != 'orphan':
                continue

            name = ''
            for alias in aliases:
                if alias[:2] == 'SN':
                    name = alias
            if not name:
                name = psname
            events, name = Events.add_event(tasks, args, events, name, log)
            sources = [events[name].add_source(srcname='Pan-STARRS 3Pi',
                                               url='http://psweb.mp.qub.ac.uk/ps1threepi/psdb/')]
            events[name].add_quantity('alias', name, sources[0])
            for ref in refs:
                sources.append(events[name].add_source(srcname=ref[0], url=ref[1]))
            source = uniq_cdl(sources)
            for alias in aliases:
                newalias = alias
                if alias[:3] in ['CSS', 'SSS', 'MLS']:
                    newalias = alias.replace('-', ':', 1)
                newalias = newalias.replace('PSNJ', 'PSN J')
                events[name].add_quantity('alias', newalias, source)
            events[name].add_quantity('ra', ra, source)
            events[name].add_quantity('dec', dec, source)
            events[name].add_quantity('claimedtype', ctype, source)

            fname2 = os.path.join(PATH.REPO_EXTERNAL, '3pi/candidate-')
            fname2 += pslink.rstrip('/').split('/')[-1] + '.html'
            if task_obj.load_archive(args) and os.path.isfile(fname2):
                with open(fname2, 'r') as f:
                    html2 = f.read()
            elif not offline:
                pslink = 'http://psweb.mp.qub.ac.uk/ps1threepi/psdb/public/' + pslink
                with open(fname2, 'w') as f:
                    response2 = urllib.request.urlopen(pslink)
                    html2 = response2.read().decode('utf-8')
                    f.write(html2)
            else:
                continue

            bs2 = BeautifulSoup(html2, 'html5lib')
            scripts = bs2.findAll('script')
            nslines = []
            nslabels = []
            for script in scripts:
                if 'jslcdata.push' not in script.text:
                    continue
                slines = script.text.splitlines()
                for line in slines:
                    if 'jslcdata.push' in line:
                        json_fname = line.strip().replace('jslcdata.push(', '').replace(');', '')
                        nslines.append(json.loads(json_fname))
                    if 'jslabels.push' in line and 'blanks' not in line and 'non det' not in line:
                        json_fname = line.strip().replace('jslabels.push(', '').replace(');', '')
                        nslabels.append(json.loads(json_fname)['label'])
            for li, line in enumerate(nslines[:len(nslabels)]):
                if not line:
                    continue
                for obs in line:
                    add_photometry(
                        events, name, time=str(obs[0]), band=nslabels[li], magnitude=str(obs[1]),
                        e_magnitude=str(obs[2]), source=source, telescope=teles)
            for li, line in enumerate(nslines[2*len(nslabels):]):
                if not line:
                    continue
                for obs in line:
                    add_photometry(
                        events, name, time=str(obs[0]), band=nslabels[li], magnitude=str(obs[1]),
                        upperlimit=True, source=source, telescope=teles)
            assoctab = bs2.find('table', {'class': 'generictable'})
            hostname = ''
            redshift = ''
            if assoctab:
                trs = assoctab.findAll('tr')
                headertds = [x.contents[0] for x in trs[1].findAll('td')]
                tds = trs[1].findAll('td')
                for tdi, td in enumerate(tds):
                    if tdi == 1:
                        hostname = td.contents[0].strip()
                    elif tdi == 4:
                        if 'z' in headertds:
                            redshift = td.contents[0].strip()
            # Skip galaxies with just SDSS id
            if is_number(hostname):
                continue
            events[name].add_quantity('host', hostname, source)
            if redshift:
                events[name].add_quantity('redshift', redshift, source, kind='host')
            if args.update:
                events, stubs = Events.journal_events(tasks, args, events, stubs, log)
        events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    return events
