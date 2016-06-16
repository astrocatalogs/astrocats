"""Import data from the Palomar Transient Factory (PTF).
"""
from bs4 import BeautifulSoup
import os
import requests

from scripts import PATH
from .. funcs import add_event, add_source, add_quantity, archived_task, \
    journal_events
from ... utils import is_number


def do_ptf(events, args, tasks):
    # response = urllib.request.urlopen('http://wiserep.weizmann.ac.il/objects/list')
    # bs = BeautifulSoup(response, 'html5lib')
    # select = bs.find('select', {'name': 'objid'})
    # options = select.findAll('option')
    # for option in options:
    #    print(option.text)
    #    name = option.text
    #    if ((name.startswith('PTF') and is_number(name[3:5])) or
    #        name.startswith('PTFS') or name.startswith('iPTF')):
    #        name = add_event(tasks, args, events, name)

    if archived_task(tasks, args, 'ptf'):
        with open(os.path.join(PATH.REPO_EXTERNAL, 'PTF/update.html'), 'r') as f:
            html = f.read()
    else:
        session = requests.Session()
        response = session.get('http://wiserep.weizmann.ac.il/spectra/update')
        html = response.text
        with open(os.path.join(PATH.REPO_EXTERNAL, 'PTF/update.html'), 'w') as f:
            f.write(html)

    bs = BeautifulSoup(html, 'html5lib')
    select = bs.find('select', {'name': 'objid'})
    options = select.findAll('option')
    for option in options:
        name = option.text
        if (((name.startswith('PTF') and is_number(name[3:5])) or
             name.startswith('PTFS') or name.startswith('iPTF'))):
            if '(' in name:
                alias = name.split('(')[0].strip(' ')
                name = name.split('(')[-1].strip(') ').replace('sn', 'SN')
                name = add_event(tasks, args, events, name)
                source = add_source(events, name, bibcode='2012PASP..124..668Y')
                add_quantity(events, name, 'alias', alias, source)
            else:
                name = add_event(tasks, args, events, name)

    with open(os.path.join(PATH.REPO_EXTERNAL, 'PTF/old-ptf-events.csv')) as f:
        for suffix in f.read().splitlines():
            name = add_event(tasks, args, events, 'PTF' + suffix)
    with open(os.path.join(PATH.REPO_EXTERNAL, 'PTF/perly-2016.csv')) as f:
        for row in f.read().splitlines():
            cols = [x.strip() for x in row.split(',')]
            alias = ''
            if cols[8]:
                name = cols[8]
                alias = 'PTF' + cols[0]
            else:
                name = 'PTF' + cols[0]
            name = add_event(tasks, args, events, name)
            source = add_source(events, name, bibcode='2016arXiv160408207P')
            add_quantity(events, name, 'alias', name, source)
            if alias:
                add_quantity(events, name, 'alias', alias, source)
            add_quantity(events, name, 'ra', cols[1], source)
            add_quantity(events, name, 'dec', cols[2], source)
            add_quantity(events, name, 'claimedtype', 'SLSN-' + cols[3], source)
            add_quantity(events, name, 'redshift', cols[4], source, kind='spectroscopic')
            maxdate = cols[6].replace('-', '/')
            upl = maxdate.startswith('<')
            add_quantity(events, name, 'maxdate', maxdate.lstrip('<'), source, upperlimit=upl)
            add_quantity(events, name, 'ebv', cols[7], source, kind='spectroscopic')
            name = add_event(tasks, args, events, 'PTF' + suffix)

    events = journal_events(tasks, args, events)
    return events
