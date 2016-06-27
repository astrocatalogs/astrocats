"""General data import tasks.
"""
import csv
import os

from scripts import PATH
from scripts.utils import pbar

from .. import Events
from ..funcs import load_cached_url


def do_grb(catalog):
    current_task = 'GRB'
    file_path = os.path.join(PATH.REPO_EXTERNAL, 'GRB-catalog/catalog.csv')
    csvtxt = load_cached_url(args,
                             current_task,
                             ('http://grb.pa.msu.edu/grbcatalog/'
                              'download_data?cut_0_min=10&cut_0=BAT%20T90'
                              '&cut_0_max=100000&num_cuts=1&no_date_cut=True'),
                             file_path)
    if not csvtxt:
        return events
    data = list(csv.reader(csvtxt.splitlines(), delimiter=',',
                           quotechar='"', skipinitialspace=True))
    for r, row in enumerate(pbar(data, current_task)):
        if r == 0:
            continue
        (events,
         name,
         source) = Events.new_event(tasks, args, events, 'GRB ' +
                                    row[0], log,
                                    srcname='Gamma-ray Bursts Catalog',
                                    url='http://grbcatalog.org')
        events[name].add_quantity('ra', row[2], source, unit='floatdegrees')
        events[name].add_quantity('dec', row[3], source, unit='floatdegrees')
        events[name].add_quantity('redshift', row[8], source)

    events = Events.journal_events(tasks, args, events, log)
    return events
