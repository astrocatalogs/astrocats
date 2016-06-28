"""General data import tasks.
"""
import os

from bs4 import BeautifulSoup

from scripts import PATH

from .. import Events
from ...utils import pbar
from ..funcs import load_cached_url


def do_asassn(catalog):
    current_task = catalog.current_task
    asn_url = 'http://www.astronomy.ohio-state.edu/~assassin/sn_list.html'
    html = load_cached_url(catalog.args, current_task, asn_url, os.path.join(
        PATH.REPO_EXTERNAL, 'ASASSN/sn_list.html'))
    if not html:
        return
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
                name = Events.add_event(td.text.strip())
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

        sources = [catalog.events[name].add_source(
            url=asn_url, srcname='ASAS-SN Supernovae')]
        typesources = sources[:]
        if atellink:
            sources.append(
                (catalog.events[name]
                 .add_source(srcname='ATel ' +
                             atellink.split('=')[-1], url=atellink)))
        if typelink:
            typesources.append(
                (catalog.events[name]
                 .add_source(srcname='ATel ' +
                             typelink.split('=')[-1], url=typelink)))
        sources = ','.join(sources)
        typesources = ','.join(typesources)
        catalog.events[name].add_quantity('alias', name, sources)
        catalog.events[name].add_quantity('discoverdate', discdate, sources)
        catalog.events[name].add_quantity('ra', ra, sources,
                                          unit='floatdegrees')
        catalog.events[name].add_quantity('dec', dec, sources,
                                          unit='floatdegrees')
        catalog.events[name].add_quantity('redshift', redshift, sources)
        catalog.events[name].add_quantity(
            'hostoffsetang', hostoff, sources, unit='arcseconds')
        for ct in claimedtype.split('/'):
            if ct != 'Unk':
                catalog.events[name].add_quantity('claimedtype', ct,
                                                  typesources)
        if host != 'Uncatalogued':
            catalog.events[name].add_quantity('host', host, sources)
    catalog.journal_events()
    return
