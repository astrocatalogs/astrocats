"""General data import tasks.
"""
from bs4 import BeautifulSoup
import os

from .. scripts import PATH
from ... utils import pbar
from .. funcs import add_event, add_source, add_quantity, \
    journal_events, load_cached_url


def do_asassn(events, args, tasks):
    current_task = 'asassn'
    asn_url = 'http://www.astronomy.ohio-state.edu/~assassin/sn_list.html'
    html = load_cached_url(args, asn_url, os.path.join(PATH.REPO_EXTERNAL, 'ASASSN/sn_list.html'))
    if not html:
        return events
    bs = BeautifulSoup(html, 'html5lib')
    trs = bs.find('table').findAll('tr')
    for tri, tr in enumerate(pbar(trs, current_task)):
        name = ''
        ra = ''
        dec = ''
        redshift = ''
        hostoff = ''
        claimedtype = ''
        host = ''
        atellink = ''
        typelink = ''
        if tri == 0:
            continue
        tds = tr.findAll('td')
        for tdi, td in enumerate(tds):
            if tdi == 1:
                name = add_event(tasks, args, events, td.text.strip())
                atellink = td.find('a')
                if atellink:
                    atellink = atellink['href']
                else:
                    atellink = ''
            if tdi == 2:
                discdate = td.text.replace('-', '/')
            if tdi == 3:
                ra = td.text
            if tdi == 4:
                dec = td.text
            if tdi == 5:
                redshift = td.text
            if tdi == 8:
                hostoff = td.text
            if tdi == 9:
                claimedtype = td.text
                typelink = td.find('a')
                if typelink:
                    typelink = typelink['href']
                else:
                    typelink = ''
            if tdi == 12:
                host = td.text

        sources = [add_source(events, name, url=asn_url, refname='ASAS-SN Supernovae')]
        typesources = sources[:]
        if atellink:
            sources.append(
                add_source(events, name, refname='ATel ' + atellink.split('=')[-1], url=atellink))
        if typelink:
            typesources.append(
                add_source(events, name, refname='ATel ' + typelink.split('=')[-1], url=typelink))
        sources = ','.join(sources)
        typesources = ','.join(typesources)
        add_quantity(events, name, 'alias', name, sources)
        add_quantity(events, name, 'discoverdate', discdate, sources)
        add_quantity(events, name, 'ra', ra, sources, unit='floatdegrees')
        add_quantity(events, name, 'dec', dec, sources, unit='floatdegrees')
        add_quantity(events, name, 'redshift', redshift, sources)
        add_quantity(events, name, 'hostoffset', hostoff, sources, unit='arcseconds')
        for ct in claimedtype.split('/'):
            if ct != 'Unk':
                add_quantity(events, name, 'claimedtype', ct, typesources)
        if host != 'Uncatalogued':
            add_quantity(events, name, 'host', host, sources)
    events = journal_events(tasks, args, events)
    return events
