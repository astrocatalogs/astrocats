"""General data import tasks.
"""
import csv
import os
import re
import urllib

from cdecimal import Decimal
from scripts import PATH
from scripts.utils import pbar

from .. import Events
from ..funcs import add_photometry, jd_to_mjd, load_cached_url


def do_gaia(catalog):
    current_task = catalog.current_task
    fname = os.path.join(PATH.REPO_EXTERNAL, 'GAIA/alerts.csv')
    csvtxt = load_cached_url(catalog.args, current_task,
                             'http://gsaweb.ast.cam.ac.uk/alerts/alerts.csv',
                             fname)
    if not csvtxt:
        return
    tsvin = list(csv.reader(csvtxt.splitlines(),
                            delimiter=',', skipinitialspace=True))
    reference = 'Gaia Photometric Science Alerts'
    refurl = 'http://gsaweb.ast.cam.ac.uk/alerts/alertsindex'
    for ri, row in enumerate(pbar(tsvin, current_task)):
        if ri == 0 or not row:
            continue
        name = catalog.add_event(row[0])
        source = catalog.events[name].add_source(srcname=reference, url=refurl)
        catalog.events[name].add_quantity('alias', name, source)
        year = '20' + re.findall(r'\d+', row[0])[0]
        catalog.events[name].add_quantity('discoverdate', year, source)
        catalog.events[name].add_quantity('ra', row[2], source, unit='floatdegrees')
        catalog.events[name].add_quantity('dec', row[3], source, unit='floatdegrees')
        if row[7] and row[7] != 'unknown':
            type = row[7].replace('SNe', '').replace('SN', '').strip()
            catalog.events[name].add_quantity('claimedtype', type, source)
        elif any([xx in row[9].upper() for xx in
                  ['SN CANDIATE', 'CANDIDATE SN', 'HOSTLESS SN']]):
            catalog.events[name].add_quantity('claimedtype', 'Candidate', source)

        if ('aka' in row[9].replace('gakaxy', 'galaxy').lower() and
                'AKARI' not in row[9]):
            commentsplit = (row[9]
                            .replace('_', ' ')
                            .replace('MLS ', 'MLS')
                            .replace('CSS ', 'CSS')
                            .replace('SN iPTF', 'iPTF')
                            .replace('SN ', 'SN')
                            .replace('AT ', 'AT'))
            commentsplit = commentsplit.split()
            for csi, cs in enumerate(commentsplit):
                if 'aka' in cs.lower() and csi < len(commentsplit) - 1:
                    alias = commentsplit[
                        csi + 1].strip('(),:.').replace('PSNJ', 'PSN J')
                    if alias[:6] == 'ASASSN' and alias[6] != '-':
                        alias = 'ASASSN-' + alias[6:]
                    catalog.events[name].add_quantity('alias', alias, source)
                    break

        fname = os.path.join(PATH.REPO_EXTERNAL, 'GAIA/') + row[0] + '.csv'
        if task_obj.load_archive(args) and os.path.isfile(fname):
            with open(fname, 'r') as ff:
                csvtxt = ff.read()
        else:
            response = urllib.request.urlopen('http://gsaweb.ast.cam.ac.uk/'
                                              'alerts/alert/' +
                                              row[0] + '/lightcurve.csv')
            with open(fname, 'w') as ff:
                csvtxt = response.read().decode('utf-8')
                ff.write(csvtxt)

        tsvin2 = csv.reader(csvtxt.splitlines())
        for ri2, row2 in enumerate(tsvin2):
            if ri2 <= 1 or not row2:
                continue
            mjd = str(jd_to_mjd(Decimal(row2[1].strip())))
            magnitude = row2[2].strip()
            if magnitude == 'null':
                continue
            e_mag = 0.
            telescope = 'GAIA'
            band = 'G'
            add_photometry(catalog.events, name, time=mjd, telescope=telescope,
                           band=band, magnitude=magnitude,
                           e_magnitude=e_mag, source=source)
        if args.update:
            events = Events.journal_events(
                tasks, args, events, log)
    catalog.journal_events()
    return
