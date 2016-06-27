"""General data import tasks.
"""
import os
import re
from glob import glob

from astropy.time import Time as astrotime

from cdecimal import Decimal
from scripts import PATH
from scripts.utils import pbar

from .. import Events
from ..funcs import add_spectrum, event_exists


def do_superfit_spectra(events, args, tasks, task_obj, log):
    from .. funcs import get_max_light, get_preferred_name
    superfit_url = 'http://www.dahowell.com/superfit.html'
    current_task = task_obj.current_task(args)
    sfdirs = list(glob(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'superfit/*')))
    for sfdir in pbar(sfdirs, desc=current_task):
        sffiles = sorted(glob(sfdir + '/*.dat'))
        lastname = ''
        oldname = ''
        for sffile in pbar(sffiles, desc=current_task):
            basename = os.path.basename(sffile)
            name = basename.split('.')[0]
            if name.startswith('sn'):
                name = 'SN' + name[2:]
                if len(name) == 7:
                    name = name[:6] + name[6].upper()
            elif name.startswith('ptf'):
                name = 'PTF' + name[3:]

            if 'theory' in name:
                continue
            if event_exists(events, name):
                prefname = get_preferred_name(events, name)
                if 'spectra' in events[prefname] and lastname != prefname:
                    continue
            if oldname and name != oldname:
                events = Events.journal_events(
                    tasks, args, events, log)
            oldname = name
            events, name = Events.add_event(tasks, args, events, name, log)
            epoch = basename.split('.')[1]
            (mldt, mlmag, mlband, mlsource) = get_max_light(events, name)
            if mldt:
                if epoch == 'max':
                    epoff = Decimal(0.0)
                elif epoch[0] == 'p':
                    epoff = Decimal(epoch[1:])
                else:
                    epoff = -Decimal(epoch[1:])
            else:
                epoff = ''

            source = events[name].add_source(
                srcname='Superfit', url=superfit_url, secondary=True)
            events[name].add_quantity('alias', oldname, source)

            with open(sffile) as ff:
                rows = ff.read().splitlines()
            specdata = []
            for row in rows:
                if row.strip():
                    specdata.append(
                        list(filter(None, re.split('\t+|\s+', row,
                                                   maxsplit=0))))
            specdata = [[xx.replace('D', 'E') for xx in list(ii)]
                        for ii in zip(*specdata)]
            wavelengths = specdata[0]
            fluxes = specdata[1]

            if epoff != '':
                mlmjd = astrotime(
                    '-'.join([str(mldt.year), str(mldt.month),
                              str(mldt.day)])).mjd
                mlmjd = str(Decimal(mlmjd) + epoff)
            else:
                mlmjd = ''
            add_spectrum(events, name, 'Angstrom', 'Uncalibrated', u_time='MJD'
                         if mlmjd else '', time=mlmjd,
                         wavelengths=wavelengths, fluxes=fluxes, source=source)

            lastname = name

        events = Events.journal_events(tasks, args, events, log)
    return events
