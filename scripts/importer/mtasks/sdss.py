"""General data import tasks.
"""
import csv
import os
import re
from glob import glob

from scripts import PATH
from scripts.utils import pbar_strings

from .. import Events
from ..funcs import add_photometry


def do_sdss(catalog):
    current_task = catalog.current_task
    with open(os.path.join(PATH.REPO_EXTERNAL,
                           'SDSS/2010ApJ...708..661D.txt'), 'r') as sdss_file:
        bibcodes2010 = sdss_file.read().split('\n')
    sdssbands = ['u', 'g', 'r', 'i', 'z']
    file_names = list(glob(os.path.join(PATH.REPO_EXTERNAL, 'SDSS/*.sum')))
    for fname in pbar_strings(file_names, desc=current_task):
        tsvin = csv.reader(open(fname, 'r'), delimiter=' ',
                           skipinitialspace=True)
        basename = os.path.basename(fname)
        if basename in bibcodes2010:
            bibcode = '2010ApJ...708..661D'
        else:
            bibcode = '2008AJ....136.2306H'

        for rr, row in enumerate(tsvin):
            if rr == 0:
                if row[5] == 'RA:':
                    name = 'SDSS-II SN ' + row[3]
                else:
                    name = 'SN' + row[5]
                name = catalog.add_event(name)
                source = catalog.events[name].add_source(bibcode=bibcode)
                catalog.events[name].add_quantity('alias', name, source)
                catalog.events[name].add_quantity(
                    'alias', 'SDSS-II SN ' + row[3], source)

                if row[5] != 'RA:':
                    year = re.findall(r'\d+', name)[0]
                    catalog.events[name].add_quantity('discoverdate', year, source)

                catalog.events[name].add_quantity(
                    'ra', row[-4], source, unit='floatdegrees')
                catalog.events[name].add_quantity(
                    'dec', row[-2], source, unit='floatdegrees')
            if rr == 1:
                error = row[4] if float(row[4]) >= 0.0 else ''
                catalog.events[name].add_quantity('redshift', row[2], source,
                                          error=error,
                                          kind='heliocentric')
            if rr >= 19:
                # Skip bad measurements
                if int(row[0]) > 1024:
                    continue

                mjd = row[1]
                band = sdssbands[int(row[2])]
                magnitude = row[3]
                e_mag = row[4]
                telescope = 'SDSS'
                add_photometry(catalog.events, name, time=mjd, telescope=telescope,
                               band=band, magnitude=magnitude,
                               e_magnitude=e_mag, source=source, system='SDSS')

    catalog.journal_events()
    return
