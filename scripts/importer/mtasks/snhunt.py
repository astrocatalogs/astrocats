"""Import tasks for the Supernova Hunt.
"""
import os
import re

from bs4 import BeautifulSoup
from scripts import PATH
from scripts.utils import pbar

from ..funcs import load_cached_url


def do_snhunt(catalog):
    current_task = catalog.get_current_task_str()
    snh_url = 'http://nesssi.cacr.caltech.edu/catalina/current.html'
    html = load_cached_url(catalog.args, current_task, snh_url, os.path.join(
        PATH.REPO_EXTERNAL, 'SNhunt/current.html'))
    if not html:
        return
    text = html.splitlines()
    findtable = False
    for ri, row in enumerate(text):
        if 'Supernova Discoveries' in row:
            findtable = True
        if findtable and '<table' in row:
            tstart = ri + 1
        if findtable and '</table>' in row:
            tend = ri - 1
    tablestr = '<html><body><table>'
    for row in text[tstart:tend]:
        if row[:3] == 'tr>':
            tablestr = tablestr + '<tr>' + row[3:]
        else:
            tablestr = tablestr + row
    tablestr = tablestr + '</table></body></html>'
    bs = BeautifulSoup(tablestr, 'html5lib')
    trs = bs.find('table').findAll('tr')
    for tr in pbar(trs, current_task):
        cols = [str(xx.text) for xx in tr.findAll('td')]
        if not cols:
            continue
        name = re.sub('<[^<]+?>', '', cols[4]
                      ).strip().replace(' ', '').replace('SNHunt', 'SNhunt')
        name = catalog.add_entry(name)
        source = catalog.events[name].add_source(srcname='Supernova Hunt', url=snh_url)
        catalog.events[name].add_quantity('alias', name, source)
        host = re.sub('<[^<]+?>', '', cols[1]).strip().replace('_', ' ')
        catalog.events[name].add_quantity('host', host, source)
        catalog.events[name].add_quantity('ra', cols[2], source, unit='floatdegrees')
        catalog.events[name].add_quantity('dec', cols[3], source, unit='floatdegrees')
        dd = cols[0]
        discoverdate = dd[:4] + '/' + dd[4:6] + '/' + dd[6:8]
        catalog.events[name].add_quantity('discoverdate', discoverdate, source)
        discoverers = cols[5].split('/')
        for discoverer in discoverers:
            catalog.events[name].add_quantity('discoverer', 'CRTS', source)
            catalog.events[name].add_quantity('discoverer', discoverer, source)
        if catalog.args.update:
            catalog.journal_events()

    catalog.journal_events()
    return
