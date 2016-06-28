"""General data import tasks.
"""
import os
from collections import OrderedDict
from glob import glob

from scripts import PATH
from scripts.utils import pbar_strings

from ..Events import KEYS


def do_external_radio(catalog):
    current_task = catalog.get_current_task_str()
    path_pattern = os.path.join(PATH.REPO_EXTERNAL_RADIO, '*.txt')
    for datafile in pbar_strings(glob(path_pattern), desc=current_task):
        oldname = os.path.basename(datafile).split('.')[0]
        name = catalog.add_event(oldname)
        radiosourcedict = OrderedDict()
        with open(datafile, 'r') as ff:
            for li, line in enumerate([xx.strip() for xx in
                                       ff.read().splitlines()]):
                if line.startswith('(') and li <= len(radiosourcedict):
                    key = line.split()[0]
                    bibc = line.split()[-1]
                    radiosourcedict[key] = catalog.events[
                        name].add_source(bibcode=bibc)
                elif li in [xx + len(radiosourcedict) for xx in range(3)]:
                    continue
                else:
                    cols = list(filter(None, line.split()))
                    source = radiosourcedict[cols[6]]
                    catalog.events[name].add_photometry(
                        time=cols[0], frequency=cols[2], u_frequency='GHz',
                        fluxdensity=cols[3], e_fluxdensity=cols[4],
                        u_fluxdensity='µJy',
                        instrument=cols[5], source=source)
                    catalog.events[name].add_quantity('alias', oldname, source)

    catalog.journal_events()
    return


def do_external_xray(catalog):
    current_task = catalog.get_current_task_str()
    path_pattern = os.path.join(PATH.REPO_EXTERNAL_XRAY, '*.txt')
    for datafile in pbar_strings(glob(path_pattern), desc=current_task):
        oldname = os.path.basename(datafile).split('.')[0]
        name = catalog.add_event(oldname)
        with open(datafile, 'r') as ff:
            for li, line in enumerate(ff.read().splitlines()):
                if li == 0:
                    source = catalog.events[name].add_source(
                        bibcode=line.split()[-1])
                elif li in [1, 2, 3]:
                    continue
                else:
                    cols = list(filter(None, line.split()))
                    catalog.events[name].add_photometry(
                        time=cols[:2],
                        energy=cols[2:4], u_energy='keV', counts=cols[4],
                        flux=cols[6],
                        unabsorbedflux=cols[8], u_flux='ergs/ss/cm^2',
                        photonindex=cols[15], instrument=cols[
                            17], nhmw=cols[11],
                        upperlimit=(float(cols[5]) < 0), source=source)
                    catalog.events[name].add_quantity('alias', oldname, source)

    catalog.journal_events()
    return


def do_internal(catalog):
    """Load events from files in the 'internal' repository, and save them.
    """
    from ..Events import EVENT
    current_task = catalog.get_current_task_str()
    path_pattern = os.path.join(PATH.REPO_INTERNAL, '*.json')
    files = glob(path_pattern)
    catalog.log.debug("found {} files matching '{}'".format(
        len(files), path_pattern))
    for datafile in pbar_strings(files, desc=current_task):
        new_event = EVENT.init_from_file(path=datafile, clean=True)
        catalog.events.update({new_event[KEYS.NAME]: new_event})

    return
