"""General data import tasks.
"""
import os
from collections import OrderedDict
from glob import glob

from .. scripts import PATH
from ... utils import pbar_strings
from .. funcs import add_event, add_photometry, add_source, add_quantity, \
    load_event_from_file, journal_events


def do_internal(events, args, tasks):
    """Load events from files in the 'internal' repository, and save them.
    """
    current_task = 'Internal'
    path_pattern = os.path.join(PATH.REPO_INTERNAL, '*.json')
    files = glob(path_pattern)
    for datafile in pbar_strings(files, desc=current_task):
        if args.update:
            if not load_event_from_file(events, args, tasks, path=datafile,
                                        clean=True, delete=False, append=True):
                raise IOError('Failed to find specified file.')
        else:
            if not load_event_from_file(events, args, tasks, path=datafile,
                                        clean=True, delete=False):
                raise IOError('Failed to find specified file.')

    events = journal_events(tasks, args, events)
    return events


def do_external_radio(events, args, tasks):
    current_task = 'External Radio'
    path_pattern = os.path.join(PATH.REPO_EXTERNAL_RADIO, '*.txt')
    for datafile in pbar_strings(glob(path_pattern), desc=current_task):
        name = add_event(tasks, args, events, os.path.basename(datafile).split('.')[0])
        radiosourcedict = OrderedDict()
        with open(datafile, 'r') as f:
            for li, line in enumerate([x.strip() for x in f.read().splitlines()]):
                if line.startswith('(') and li <= len(radiosourcedict):
                    radiosourcedict[line.split()[0]] = add_source(events, name, bibcode=line.split()[-1])
                elif li in [x + len(radiosourcedict) for x in range(3)]:
                    continue
                else:
                    cols = list(filter(None, line.split()))
                    source = radiosourcedict[cols[6]]
                    add_photometry(
                        events, name, time=cols[0], frequency=cols[2], u_frequency='GHz', fluxdensity=cols[3],
                        e_fluxdensity=cols[4], u_fluxdensity='µJy', instrument=cols[5], source=source)
                    add_quantity(events, name, 'alias', name, source)

    events = journal_events(tasks, args, events)
    return events


def do_external_xray(events, args, tasks):
    current_task = 'External X-ray'
    path_pattern = os.path.join(PATH.REPO_EXTERNAL_XRAY, '*.txt')
    for datafile in pbar_strings(glob(path_pattern), desc=current_task):
        name = add_event(tasks, args, events, os.path.basename(datafile).split('.')[0])
        with open(datafile, 'r') as f:
            for li, line in enumerate(f.read().splitlines()):
                if li == 0:
                    source = add_source(events, name, bibcode=line.split()[-1])
                elif li in [1, 2, 3]:
                    continue
                else:
                    cols = list(filter(None, line.split()))
                    add_photometry(
                        events, name, time=cols[:2],
                        energy=cols[2:4], u_energy='keV', counts=cols[4], flux=cols[6],
                        unabsorbedflux=cols[8], u_flux='ergs/s/cm^2',
                        photonindex=cols[15], instrument=cols[17], nhmw=cols[11],
                        upperlimit=(float(cols[5]) < 0), source=source)
                    add_quantity(events, name, 'alias', name, source)

    events = journal_events(tasks, args, events)
    return events


'''
def do_simbad(events, args, tasks):
    Simbad.list_votable_fields()
    customSimbad = Simbad()
    customSimbad.add_votable_fields('otype', 'id(opt)')
    result = customSimbad.query_object('SN 20[0-9][0-9]*', wildcard=True)
    for r, row in enumerate(result):
        if row['OTYPE'].decode() != 'SN':
            continue
        name = row['MAIN_ID'].decode()
        aliases = Simbad.query_objectids(name)
        print(aliases)
        if name[:3] == 'SN ':
            name = 'SN' + name[3:]
        if name[:2] == 'SN' and is_number(name[2:]):
            name = name + 'A'
        name = add_event(tasks, args, events, name)
    events = journal_events(tasks, args, events)
    return events
'''
