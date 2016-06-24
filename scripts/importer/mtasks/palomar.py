"""Import data from the Palomar Transient Factory (PTF).
"""
from bs4 import BeautifulSoup
import os
import requests

from scripts import PATH
from .. import Events
from ... utils import is_number


def do_ptf(events, stubs, args, tasks, task_obj, log):
    # response = urllib.request.urlopen('http://wiserep.weizmann.ac.il/objects/list')
    # bs = BeautifulSoup(response, 'html5lib')
    # select = bs.find('select', {'name': 'objid'})
    # options = select.findAll('option')
    # for option in options:
    #    print(option.text)
    #    name = option.text
    #    if ((name.startswith('PTF') and is_number(name[3:5])) or
    #        name.startswith('PTFS') or name.startswith('iPTF')):
    # events, name = Events.add_event(tasks, args, events, name, log)

    if task_obj.load_archive(args):
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
                events, name = Events.add_event(tasks, args, events, name, log)
                source = events[name].add_source(bibcode='2012PASP..124..668Y')
                events[name].add_quantity('alias', alias, source)
            else:
                # events, name = Events.add_event(tasks, args, events, name, log)
                events, name, source = Events.new_event(bibcode='2012PASP..124..668Y')

    with open(os.path.join(PATH.REPO_EXTERNAL, 'PTF/old-ptf-events.csv')) as f:
        for suffix in f.read().splitlines():
            events, name = Events.add_event(tasks, args, events, 'PTF' + suffix, log)
    with open(os.path.join(PATH.REPO_EXTERNAL, 'PTF/perly-2016.csv')) as f:
        for row in f.read().splitlines():
            cols = [x.strip() for x in row.split(',')]
            alias = ''
            if cols[8]:
                name = cols[8]
                alias = 'PTF' + cols[0]
            else:
                name = 'PTF' + cols[0]
            events, name = Events.add_event(tasks, args, events, name, log)
            source = events[name].add_source(bibcode='2016arXiv160408207P')
            events[name].add_quantity('alias', name, source)
            if alias:
                events[name].add_quantity('alias', alias, source)
            events[name].add_quantity('ra', cols[1], source)
            events[name].add_quantity('dec', cols[2], source)
            events[name].add_quantity('claimedtype', 'SLSN-' + cols[3], source)
            events[name].add_quantity('redshift', cols[4], source, kind='spectroscopic')
            maxdate = cols[6].replace('-', '/')
            upl = maxdate.startswith('<')
            events[name].add_quantity('maxdate', maxdate.lstrip('<'), source, upperlimit=upl)
            events[name].add_quantity('ebv', cols[7], source, kind='spectroscopic')
            events, name = Events.add_event(tasks, args, events, 'PTF' + suffix, log)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events
