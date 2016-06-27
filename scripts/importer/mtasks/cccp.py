"""General data import tasks.
"""
import csv
import json
import os
from collections import OrderedDict
from glob import glob

import requests
from bs4 import BeautifulSoup

from cdecimal import Decimal
from scripts import PATH
from scripts.utils import is_number, pbar, pbar_strings, round_sig

from .. import Events
from ..funcs import add_photometry, event_exists, load_cached_url, uniq_cdl


def do_cccp(catalog):
    current_task = task_obj.current_task(args)
    cccpbands = ['B', 'V', 'R', 'I']
    file_names = list(
        glob(os.path.join(PATH.REPO_EXTERNAL, 'CCCP/apj407397*.txt')))
    for datafile in pbar_strings(file_names, current_task + ': apj407397...'):
        with open(datafile, 'r') as ff:
            tsvin = csv.reader(ff, delimiter='\t', skipinitialspace=True)
            for rr, row in enumerate(tsvin):
                if rr == 0:
                    continue
                elif rr == 1:
                    name = 'SN' + row[0].split('SN ')[-1]
                    events, name = Events.add_event(
                        tasks, args, events, name, log)
                    source = events[name].add_source(
                        bibcode='2012ApJ...744...10K')
                    events[name].add_quantity('alias', name, source)
                elif rr >= 5:
                    mjd = str(Decimal(row[0]) + 53000)
                    for bb, band in enumerate(cccpbands):
                        if row[2 * bb + 1]:
                            mag = row[2 * bb + 1].strip('>')
                            upl = (not row[2 * bb + 2])
                            add_photometry(events, name, time=mjd, band=band,
                                           magnitude=mag,
                                           e_magnitude=row[2 * bb + 2],
                                           upperlimit=upl, source=source)

    if task_obj.load_archive(args):
        with open(os.path.join(PATH.REPO_EXTERNAL,
                               'CCCP/sc_cccp.html'), 'r') as ff:
            html = ff.read()
    else:
        session = requests.Session()
        response = session.get(
            'https://webhome.weizmann.ac.il/home/iair/sc_cccp.html')
        html = response.text
        with open(os.path.join(PATH.REPO_EXTERNAL,
                               'CCCP/sc_cccp.html'), 'w') as ff:
            ff.write(html)

    soup = BeautifulSoup(html, 'html5lib')
    links = soup.body.findAll("a")
    for link in pbar(links, current_task + ': links'):
        if 'sc_sn' in link['href']:
            events, name = Events.add_event(
                tasks, args, events, link.text.replace(' ', ''), log)
            source = (events[name]
                      .add_source(srcname='CCCP',
                                  url=('https://webhome.weizmann.ac.il'
                                       '/home/iair/sc_cccp.html')))
            events[name].add_quantity('alias', name, source)

            if task_obj.load_archive(args):
                fname = os.path.join(PATH.REPO_EXTERNAL,
                                     'CCCP/') + link['href'].split('/')[-1]
                with open(fname, 'r') as ff:
                    html2 = ff.read()
            else:
                response2 = session.get(
                    'https://webhome.weizmann.ac.il/home/iair/' + link['href'])
                html2 = response2.text
                fname = os.path.join(PATH.REPO_EXTERNAL,
                                     'CCCP/') + link['href'].split('/')[-1]
                with open(fname, 'w') as ff:
                    ff.write(html2)

            soup2 = BeautifulSoup(html2, 'html5lib')
            links2 = soup2.body.findAll("a")
            for link2 in links2:
                if '.txt' in link2['href'] and '_' in link2['href']:
                    band = link2['href'].split('_')[1].split('.')[0].upper()
                    if task_obj.load_archive(args):
                        fname = os.path.join(PATH.REPO_EXTERNAL, 'CCCP/')
                        fname += link2['href'].split('/')[-1]
                        if not os.path.isfile(fname):
                            continue
                        with open(fname, 'r') as ff:
                            html3 = ff.read()
                    else:
                        response3 = (session
                                     .get('https://webhome.weizmann.ac.il'
                                          '/home/iair/cccp/' +
                                          link2['href']))
                        if response3.status_code == 404:
                            continue
                        html3 = response3.text
                        fname = os.path.join(PATH.REPO_EXTERNAL, 'CCCP/')
                        fname += link2['href'].split('/')[-1]
                        with open(fname, 'w') as ff:
                            ff.write(html3)
                    table = [[str(Decimal(yy.strip())).rstrip('0') for yy in
                              xx.split(',')]
                             for xx in list(filter(None, html3.split('\n')))]
                    for row in table:
                        add_photometry(events, name,
                                       time=str(Decimal(row[0]) + 53000),
                                       band=band, magnitude=row[1],
                                       e_magnitude=row[2], source=source)

    events = Events.journal_events(tasks, args, events, log)
    return events


def do_cpcs(catalog):
    current_task = task_obj.current_task(args)
    cpcs_url = ('http://gsaweb.ast.cam.ac.uk/'
                'followup/list_of_alerts?format=json&num=100000&'
                'published=1&observed_only=1&'
                'hashtag=JG_530ad9462a0b8785bfb385614bf178c6')
    jsontxt = load_cached_url(args, current_task, cpcs_url, os.path.join(
        PATH.REPO_EXTERNAL, 'CPCS/index.json'))
    if not jsontxt:
        return events
    alertindex = json.loads(jsontxt, object_pairs_hook=OrderedDict)
    ids = [xx['id'] for xx in alertindex]
    for ii, ai in enumerate(pbar(ids, current_task)):
        name = alertindex[ii]['ivorn'].split('/')[-1].strip()
        # Skip aa few weird entries
        if name == 'ASASSNli':
            continue
        # Just use aa whitelist for now since naming seems inconsistent
        white_list = ['GAIA', 'OGLE', 'ASASSN', 'MASTER', 'OTJ', 'PS1', 'IPTF']
        if True in [xx in name.upper() for xx in white_list]:
            name = name.replace('Verif', '').replace('_', ' ')
            if 'ASASSN' in name and name[6] != '-':
                name = 'ASASSN-' + name[6:]
            if 'MASTEROTJ' in name:
                name = name.replace('MASTEROTJ', 'MASTER OT J')
            if 'OTJ' in name:
                name = name.replace('OTJ', 'MASTER OT J')
            if name.upper().startswith('IPTF'):
                name = 'iPTF' + name[4:]
            # Only add events that are classified as SN.
            if event_exists(events, name):
                continue
            oldname = name
            events, name = Events.add_event(tasks, args, events, name, log)
        else:
            continue

        sec_source = events[name].add_source(
            srcname='Cambridge Photometric Calibration Server',
            url='http://gsaweb.ast.cam.ac.uk/followup/', secondary=True)
        events[name].add_quantity('alias', oldname, sec_source)
        unit_deg = 'floatdegrees'
        events[name].add_quantity(
            'ra', str(alertindex[ii]['ra']), sec_source, unit=unit_deg)
        events[name].add_quantity('dec', str(
            alertindex[ii]['dec']), sec_source, unit=unit_deg)

        alerturl = ('http://gsaweb.ast.cam.ac.uk/'
                    'followup/get_alert_lc_data?alert_id=' +
                    str(ai))
        source = events[name].add_source(
            srcname='CPCS Alert ' + str(ai), url=alerturl)
        fname = os.path.join(PATH.REPO_EXTERNAL,
                             'CPCS/alert-') + str(ai).zfill(2) + '.json'
        if task_obj.load_archive(args) and os.path.isfile(fname):
            with open(fname, 'r') as ff:
                jsonstr = ff.read()
        else:
            session = requests.Session()
            response = session.get(
                alerturl + '&hashtag=JG_530ad9462a0b8785bfb385614bf178c6')
            with open(fname, 'w') as ff:
                jsonstr = response.text
                ff.write(jsonstr)

        try:
            cpcsalert = json.loads(jsonstr)
        except:
            continue

        mjds = [round_sig(xx, sig=9) for xx in cpcsalert['mjd']]
        mags = [round_sig(xx, sig=6) for xx in cpcsalert['mag']]
        errs = [round_sig(xx, sig=6) if (is_number(xx) and float(xx) > 0.0)
                else '' for xx in cpcsalert['magerr']]
        bnds = cpcsalert['filter']
        obs = cpcsalert['observatory']
        for mi, mjd in enumerate(mjds):
            add_photometry(events, name, time=mjd, magnitude=mags[mi],
                           e_magnitude=errs[mi],
                           band=bnds[mi], observatory=obs[mi],
                           source=uniq_cdl([source, sec_source]))
        if args.update:
            events = Events.journal_events(
                tasks, args, events, log)

    events = Events.journal_events(tasks, args, events, log)
    return events
