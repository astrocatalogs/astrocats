"""General data import tasks.
"""
import csv
import os

from scripts import PATH
from scripts.utils import pbar

from .. import Events


def do_fermi(catalog):
    current_task = task_obj.current_task(args)
    with open(os.path.join(PATH.REPO_EXTERNAL,
                           '1SC_catalog_v01.asc'), 'r') as ff:
        tsvin = list(csv.reader(ff, delimiter=','))
        for ri, row in enumerate(pbar(tsvin, current_task)):
            if row[0].startswith('#'):
                if len(row) > 1 and 'UPPER_LIMITS' in row[1]:
                    break
                continue
            if 'Classified' not in row[1]:
                continue
            name = row[0].replace('SNR', 'G')
            name = catalog.add_event(name)
            source = events[name].add_source(bibcode='2016ApJS..224....8A')
            events[name].add_quantity('alias', name, source)
            events[name].add_quantity(
                'alias', row[0].replace('SNR', 'MWSNR'), source)
            events[name].add_quantity(
                'ra', row[2], source, unit='floatdegrees')
            events[name].add_quantity(
                'dec', row[3], source, unit='floatdegrees')
    events = Events.journal_events(tasks, args, events, log)
    return events
