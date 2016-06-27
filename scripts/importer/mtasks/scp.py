"""General data import tasks.
"""
import csv
import os

from scripts import PATH
from scripts.utils import pbar

from .. import Events


def do_scp(events, stubs, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    tsvin = list(csv.reader(open(os.path.join(PATH.REPO_EXTERNAL, 'SCP09.csv'),
                                 'r'), delimiter=','))
    for ri, row in enumerate(pbar(tsvin, current_task)):
        if ri == 0:
            continue
        name = row[0].replace('SCP', 'SCP-')
        events, name = Events.add_event(tasks, args, events, name, log)
        source = (events[name]
                  .add_source(srcname='Supernova Cosmology Project',
                              url=('http://supernova.lbl.gov/'
                                   '2009ClusterSurvey/')))
        events[name].add_quantity('alias', name, source)
        if row[1]:
            events[name].add_quantity('alias', row[1], source)
        if row[2]:
            kind = 'spectroscopic' if row[3] == 'sn' else 'host'
            events[name].add_quantity('redshift', row[2], source, kind=kind)
        if row[4]:
            events[name].add_quantity(
                'redshift', row[2], source, kind='cluster')
        if row[6]:
            claimedtype = row[6].replace('SN ', '')
            kind = ('spectroscopic/light curve' if 'a' in row[7] and 'c' in
                    row[7] else
                    'spectroscopic' if 'a' in row[7] else
                    'light curve' if 'c' in row[7]
                    else '')
            if claimedtype != '?':
                events[name].add_quantity(
                    'claimedtype', claimedtype, source, kind=kind)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events
