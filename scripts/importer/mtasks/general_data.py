"""General data import tasks.
"""
import csv
import json
import os
import re
import urllib
from collections import OrderedDict
from datetime import datetime
from glob import glob
from html import unescape
from math import ceil, log10

import requests
from astropy.time import Time as astrotime
from bs4 import BeautifulSoup

from cdecimal import Decimal
from scripts import PATH
from scripts.utils import (is_number, pbar, pbar_strings, pretty_num,
                           round_sig, single_spaces)

from .. import Events
from ..constants import TRAVIS_QUERY_LIMIT
from ..Events import load_event_from_file
from ..funcs import (add_photometry, add_spectrum, event_exists, jd_to_mjd,
                     load_cached_url, make_date_string, uniq_cdl)


def do_cccp(events, stubs, args, tasks, task_obj, log):
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
                            add_photometry(
                                events, name, time=mjd, band=band, magnitude=mag,
                                e_magnitude=row[2 * bb + 2], upperlimit=upl, source=source)

    if task_obj.load_archive(args):
        with open(os.path.join(PATH.REPO_EXTERNAL, 'CCCP/sc_cccp.html'), 'r') as ff:
            html = ff.read()
    else:
        session = requests.Session()
        response = session.get(
            'https://webhome.weizmann.ac.il/home/iair/sc_cccp.html')
        html = response.text
        with open(os.path.join(PATH.REPO_EXTERNAL, 'CCCP/sc_cccp.html'), 'w') as ff:
            ff.write(html)

    soup = BeautifulSoup(html, 'html5lib')
    links = soup.body.findAll("a")
    for link in pbar(links, current_task + ': links'):
        if 'sc_sn' in link['href']:
            events, name = Events.add_event(
                tasks, args, events, link.text.replace(' ', ''), log)
            source = events[name].add_source(
                srcname='CCCP', url='https://webhome.weizmann.ac.il/home/iair/sc_cccp.html')
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
                        response3 = session.get('https://webhome.weizmann.ac.il/home/iair/cccp/' +
                                                link2['href'])
                        if response3.status_code == 404:
                            continue
                        html3 = response3.text
                        fname = os.path.join(PATH.REPO_EXTERNAL, 'CCCP/')
                        fname += link2['href'].split('/')[-1]
                        with open(fname, 'w') as ff:
                            ff.write(html3)
                    table = [[str(Decimal(yy.strip())).rstrip('0') for yy in xx.split(',')]
                             for xx in list(filter(None, html3.split('\n')))]
                    for row in table:
                        add_photometry(
                            events, name, time=str(Decimal(row[0]) + 53000), band=band,
                            magnitude=row[1], e_magnitude=row[2], source=source)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_cpcs(events, stubs, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    cpcs_url = ('http://gsaweb.ast.cam.ac.uk/followup/list_of_alerts?format=json&num=100000&'
                'published=1&observed_only=1&hashtag=JG_530ad9462a0b8785bfb385614bf178c6')
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

        alerturl = 'http://gsaweb.ast.cam.ac.uk/followup/get_alert_lc_data?alert_id=' + \
            str(ai)
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
            add_photometry(
                events, name, time=mjd, magnitude=mags[
                    mi], e_magnitude=errs[mi],
                band=bnds[mi], observatory=obs[mi], source=uniq_cdl([source, sec_source]))
        if args.update:
            events, stubs = Events.journal_events(
                tasks, args, events, stubs, log)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_crts(events, stubs, args, tasks, task_obj, log):
    crtsnameerrors = ['2011ax']
    current_task = task_obj.current_task(args)
    folders = ['catalina', 'MLS', 'SSS']
    for fold in pbar(folders, current_task):
        html = load_cached_url(args, current_task, 'http://nesssi.cacr.caltech.edu/' + fold + '/AllSN.html',
                               os.path.join(PATH.REPO_EXTERNAL, 'CRTS', fold + '.html'))
        if not html:
            continue
        bs = BeautifulSoup(html, 'html5lib')
        trs = bs.findAll('tr')
        for tri, tr in enumerate(pbar(trs, current_task)):
            tds = tr.findAll('td')
            if not tds:
                continue
            # refs = []
            aliases = []
            crtsname = ''
            ra = ''
            dec = ''
            lclink = ''
            # ttype = ''
            # ctype = ''
            for tdi, td in enumerate(tds):
                if tdi == 0:
                    crtsname = td.contents[0].text.strip()
                elif tdi == 1:
                    ra = td.contents[0]
                elif tdi == 2:
                    dec = td.contents[0]
                elif tdi == 11:
                    lclink = td.find('a')['onclick']
                    lclink = lclink.split("'")[1]
                elif tdi == 13:
                    aliases = re.sub('[()]', '', re.sub(
                        '<[^<]+?>', '', td.contents[-1].strip()))
                    aliases = [xx.strip('; ') for xx in list(
                        filter(None, aliases.split(' ')))]

            name = ''
            hostmag = ''
            hostupper = False
            validaliases = []
            for ai, alias in enumerate(aliases):
                if alias in ['SN', 'SDSS']:
                    continue
                if alias in crtsnameerrors:
                    continue
                if alias == 'mag':
                    if ai < len(aliases) - 1:
                        ind = ai + 1
                        if aliases[ai + 1] in ['SDSS']:
                            ind = ai + 2
                        elif aliases[ai + 1] in ['gal', 'obj', 'object', 'source']:
                            ind = ai - 1
                        if '>' in aliases[ind]:
                            hostupper = True
                        hostmag = aliases[ind].strip('>~').replace(
                            ',', '.').replace('m', '.')
                    continue
                if is_number(alias[:4]) and alias[:2] == '20' and len(alias) > 4:
                    name = 'SN' + alias
                if ((('asassn' in alias and len(alias) > 6) or
                     ('ptf' in alias and len(alias) > 3) or
                     ('ps1' in alias and len(alias) > 3) or 'snhunt' in alias or
                     ('mls' in alias and len(alias) > 3) or 'gaia' in alias or
                     ('lsq' in alias and len(alias) > 3))):
                    alias = alias.replace('SNHunt', 'SNhunt')
                    validaliases.append(alias)

            if not name:
                name = crtsname
            events, name = Events.add_event(tasks, args, events, name, log)
            source = events[name].add_source(
                srcname='Catalina Sky Survey', bibcode='2009ApJ...696..870D',
                url='http://nesssi.cacr.caltech.edu/catalina/AllSN.html')
            events[name].add_quantity('alias', name, source)
            for alias in validaliases:
                events[name].add_quantity('alias', alias, source)
            events[name].add_quantity('ra', ra, source, unit='floatdegrees')
            events[name].add_quantity('dec', dec, source, unit='floatdegrees')

            if hostmag:
                # 1.0 magnitude error based on Drake 2009 assertion that SN are only considered
                #    real if they are 2 mags brighter than host.
                add_photometry(
                    events, name, band='C', magnitude=hostmag, e_magnitude=1.0, source=source,
                    host=True, telescope='Catalina Schmidt', upperlimit=hostupper)

            fname2 = (PATH.REPO_EXTERNAL + '/' + fold + '/' +
                      lclink.split('.')[-2].rstrip('p').split('/')[-1] + '.html')
            if task_obj.load_archive(args) and os.path.isfile(fname2):
                with open(fname2, 'r') as ff:
                    html2 = ff.read()
            else:
                with open(fname2, 'w') as ff:
                    response2 = urllib.request.urlopen(lclink)
                    html2 = response2.read().decode('utf-8')
                    ff.write(html2)

            lines = html2.splitlines()
            teles = 'Catalina Schmidt'
            for line in lines:
                if 'javascript:showx' in line:
                    mjdstr = re.search("showx\('(.*?)'\)",
                                       line).group(1).split('(')[0].strip()
                    if not is_number(mjdstr):
                        continue
                    mjd = str(Decimal(mjdstr) + Decimal(53249.0))
                else:
                    continue
                if 'javascript:showy' in line:
                    mag = re.search("showy\('(.*?)'\)", line).group(1)
                if 'javascript:showz' in line:
                    err = re.search("showz\('(.*?)'\)", line).group(1)
                e_mag = err if float(err) > 0.0 else ''
                upl = (float(err) == 0.0)
                add_photometry(
                    events, name, time=mjd, band='C', magnitude=mag, source=source,
                    includeshost=True, telescope=teles, e_magnitude=e_mag, upperlimit=upl)
            if args.update:
                events, stubs = Events.journal_events(
                    tasks, args, events, stubs, log)

        if args.travis and tri > TRAVIS_QUERY_LIMIT:
            break

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_des(events, stubs, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    des_url = 'https://portal.nersc.gov/des-sn/'
    des_trans_url = des_url + 'transients/'
    ackn_url = 'http://www.noao.edu/noao/library/NOAO_Publications_Acknowledgments.html#DESdatause'
    # Make sure there is aa trailing slash
    des_path = os.path.join(PATH.REPO_EXTERNAL, 'DES', '')
    html = load_cached_url(
        args, current_task, des_trans_url, des_path + 'transients.html')
    if not html:
        return events
    bs = BeautifulSoup(html, 'html5lib')
    trs = bs.find('tbody').findAll('tr')
    for tri, tr in enumerate(pbar(trs, current_task)):
        name = ''
        # source = ''
        if tri == 0:
            continue
        tds = tr.findAll('td')
        for tdi, td in enumerate(tds):
            if tdi == 0:
                events, name = Events.add_event(
                    tasks, args, events, td.text.strip(), log)
            if tdi == 1:
                (ra, dec) = [xx.strip() for xx in td.text.split('\xa0')]
            if tdi == 6:
                atellink = td.find('a')
                if atellink:
                    atellink = atellink['href']
                else:
                    atellink = ''

        sources = [events[name].add_source(url=des_url, srcname='DES Bright Transients',
                                           acknowledgment=ackn_url)]
        if atellink:
            sources.append(
                events[name].add_source(srcname='ATel ' + atellink.split('=')[-1], url=atellink))
        sources += [events[name].add_source(bibcode='2012ApJ...753..152B'),
                    events[name].add_source(bibcode='2015AJ....150..150F'),
                    events[name].add_source(bibcode='2015AJ....150...82G'),
                    events[name].add_source(bibcode='2015AJ....150..172K')]
        sources = ','.join(sources)
        events[name].add_quantity('alias', name, sources)
        events[name].add_quantity('ra', ra, sources)
        events[name].add_quantity('dec', dec, sources)

        html2 = load_cached_url(
            args, current_task, des_trans_url + name, des_path + name + '.html')
        if not html2:
            continue
        lines = html2.splitlines()
        for line in lines:
            if 'var data = ' in line:
                jsontxt = json.loads(line.split('=')[-1].rstrip(';'))
                for ii, band in enumerate(jsontxt['band']):
                    upl = True if float(jsontxt['snr'][ii]) <= 3.0 else ''
                    add_photometry(
                        events, name, time=jsontxt['mjd'][
                            ii], magnitude=jsontxt['mag'][ii],
                        e_magnitude=jsontxt['mag_error'][ii],
                        band=band, observatory='CTIO', telescope='Blanco 4m', instrument='DECam',
                        upperlimit=upl, source=sources)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_external_radio(events, stubs, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    path_pattern = os.path.join(PATH.REPO_EXTERNAL_RADIO, '*.txt')
    for datafile in pbar_strings(glob(path_pattern), desc=current_task):
        oldname = os.path.basename(datafile).split('.')[0]
        events, name = Events.add_event(tasks, args, events, oldname, log)
        radiosourcedict = OrderedDict()
        with open(datafile, 'r') as ff:
            for li, line in enumerate([xx.strip() for xx in ff.read().splitlines()]):
                if line.startswith('(') and li <= len(radiosourcedict):
                    key = line.split()[0]
                    bibc = line.split()[-1]
                    radiosourcedict[key] = events[
                        name].add_source(bibcode=bibc)
                elif li in [xx + len(radiosourcedict) for xx in range(3)]:
                    continue
                else:
                    cols = list(filter(None, line.split()))
                    source = radiosourcedict[cols[6]]
                    add_photometry(
                        events, name, time=cols[0], frequency=cols[
                            2], u_frequency='GHz',
                        fluxdensity=cols[3], e_fluxdensity=cols[
                            4], u_fluxdensity='µJy',
                        instrument=cols[5], source=source)
                    events[name].add_quantity('alias', oldname, source)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_external_xray(events, stubs, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    path_pattern = os.path.join(PATH.REPO_EXTERNAL_XRAY, '*.txt')
    for datafile in pbar_strings(glob(path_pattern), desc=current_task):
        oldname = os.path.basename(datafile).split('.')[0]
        events, name = Events.add_event(tasks, args, events, oldname, log)
        with open(datafile, 'r') as ff:
            for li, line in enumerate(ff.read().splitlines()):
                if li == 0:
                    source = events[name].add_source(bibcode=line.split()[-1])
                elif li in [1, 2, 3]:
                    continue
                else:
                    cols = list(filter(None, line.split()))
                    add_photometry(
                        events, name, time=cols[:2],
                        energy=cols[2:4], u_energy='keV', counts=cols[4], flux=cols[6],
                        unabsorbedflux=cols[8], u_flux='ergs/ss/cm^2',
                        photonindex=cols[15], instrument=cols[
                            17], nhmw=cols[11],
                        upperlimit=(float(cols[5]) < 0), source=source)
                    events[name].add_quantity('alias', oldname, source)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_fermi(events, stubs, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    with open(os.path.join(PATH.REPO_EXTERNAL, '1SC_catalog_v01.asc'), 'r') as ff:
        tsvin = list(csv.reader(ff, delimiter=','))
        for ri, row in enumerate(pbar(tsvin, current_task)):
            if row[0].startswith('#'):
                if len(row) > 1 and 'UPPER_LIMITS' in row[1]:
                    break
                continue
            if 'Classified' not in row[1]:
                continue
            name = row[0].replace('SNR', 'G')
            events, name = Events.add_event(tasks, args, events, name, log)
            source = events[name].add_source(bibcode='2016ApJS..224....8A')
            events[name].add_quantity('alias', name, source)
            events[name].add_quantity(
                'alias', row[0].replace('SNR', 'MWSNR'), source)
            events[name].add_quantity(
                'ra', row[2], source, unit='floatdegrees')
            events[name].add_quantity(
                'dec', row[3], source, unit='floatdegrees')
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_gaia(events, stubs, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    fname = os.path.join(PATH.REPO_EXTERNAL, 'GAIA/alerts.csv')
    csvtxt = load_cached_url(
        args, current_task, 'http://gsaweb.ast.cam.ac.uk/alerts/alerts.csv', fname)
    if not csvtxt:
        return events
    tsvin = list(csv.reader(csvtxt.splitlines(),
                            delimiter=',', skipinitialspace=True))
    reference = 'Gaia Photometric Science Alerts'
    refurl = 'http://gsaweb.ast.cam.ac.uk/alerts/alertsindex'
    for ri, row in enumerate(pbar(tsvin, current_task)):
        if ri == 0 or not row:
            continue
        events, name = Events.add_event(tasks, args, events, row[0], log)
        source = events[name].add_source(srcname=reference, url=refurl)
        events[name].add_quantity('alias', name, source)
        year = '20' + re.findall(r'\d+', row[0])[0]
        events[name].add_quantity('discoverdate', year, source)
        events[name].add_quantity('ra', row[2], source, unit='floatdegrees')
        events[name].add_quantity('dec', row[3], source, unit='floatdegrees')
        if row[7] and row[7] != 'unknown':
            type = row[7].replace('SNe', '').replace('SN', '').strip()
            events[name].add_quantity('claimedtype', type, source)
        elif any([xx in row[9].upper() for xx in ['SN CANDIATE', 'CANDIDATE SN', 'HOSTLESS SN']]):
            events[name].add_quantity('claimedtype', 'Candidate', source)

        if 'aka' in row[9].replace('gakaxy', 'galaxy').lower() and 'AKARI' not in row[9]:
            commentsplit = (row[9].replace('_', ' ').replace('MLS ', 'MLS').replace('CSS ', 'CSS').
                            replace('SN iPTF', 'iPTF').replace('SN ', 'SN').replace('AT ', 'AT'))
            commentsplit = commentsplit.split()
            for csi, cs in enumerate(commentsplit):
                if 'aka' in cs.lower() and csi < len(commentsplit) - 1:
                    alias = commentsplit[
                        csi + 1].strip('(),:.').replace('PSNJ', 'PSN J')
                    if alias[:6] == 'ASASSN' and alias[6] != '-':
                        alias = 'ASASSN-' + alias[6:]
                    events[name].add_quantity('alias', alias, source)
                    break

        fname = os.path.join(PATH.REPO_EXTERNAL, 'GAIA/') + row[0] + '.csv'
        if task_obj.load_archive(args) and os.path.isfile(fname):
            with open(fname, 'r') as ff:
                csvtxt = ff.read()
        else:
            response = urllib.request.urlopen('http://gsaweb.ast.cam.ac.uk/alerts/alert/' +
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
            add_photometry(
                events, name, time=mjd, telescope=telescope, band=band, magnitude=magnitude,
                e_magnitude=e_mag, source=source)
        if args.update:
            events, stubs = Events.journal_events(
                tasks, args, events, stubs, log)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_internal(events, stubs, args, tasks, task_obj, log):
    """Load events from files in the 'internal' repository, and save them.
    """
    current_task = task_obj.current_task(args)
    path_pattern = os.path.join(PATH.REPO_INTERNAL, '*.json')
    files = glob(path_pattern)
    log.debug("found {} files matching '{}'".format(len(files), path_pattern))
    log.error("`do_internal` 'update' section is disabled")
    for datafile in pbar_strings(files, desc=current_task):
        # FIX: do we still need this difference?
        '''
        if args.update:
            if not load_event_from_file(events, args, tasks, path=datafile,
                                        clean=True, delete=False, append=True):
                raise IOError('Failed to find specified file.')
        else:
            if not load_event_from_file(events, args, tasks, path=datafile,
                                        clean=True, delete=False):
                raise IOError('Failed to find specified file.')
        '''
        new_event = load_event_from_file(events, args, tasks, log, path=datafile,
                                         clean=True, delete=False)
        events.update({new_event.name: new_event})

    return events


def do_itep(events, stubs, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    itepbadsources = ['2004ApJ...602..571B']
    needsbib = []
    with open(os.path.join(PATH.REPO_EXTERNAL, 'itep-refs.txt'), 'r') as refs_file:
        refrep = refs_file.read().splitlines()
    refrepf = dict(list(zip(refrep[1::2], refrep[::2])))
    fname = os.path.join(PATH.REPO_EXTERNAL, 'itep-lc-cat-28dec2015.txt')
    tsvin = list(csv.reader(open(fname, 'r'),
                            delimiter='|', skipinitialspace=True))
    curname = ''
    for rr, row in enumerate(pbar(tsvin, current_task)):
        if rr <= 1 or len(row) < 7:
            continue
        oldname = 'SN' + row[0].strip()
        mjd = str(jd_to_mjd(Decimal(row[1].strip())))
        band = row[2].strip()
        magnitude = row[3].strip()
        e_magnitude = row[4].strip()
        reference = row[6].strip().strip(',')

        if curname != oldname:
            curname = oldname
            events, name = Events.add_event(tasks, args, events, oldname, log)

            sec_reference = 'Sternberg Astronomical Institute Supernova Light Curve Catalogue'
            sec_refurl = 'http://dau.itep.ru/sn/node/72'
            sec_source = events[name].add_source(
                srcname=sec_reference, url=sec_refurl, secondary=True)
            events[name].add_quantity('alias', oldname, sec_source)

            year = re.findall(r'\d+', name)[0]
            events[name].add_quantity('discoverdate', year, sec_source)
        if reference in refrepf:
            bibcode = unescape(refrepf[reference])
            source = events[name].add_source(bibcode=bibcode)
        else:
            needsbib.append(reference)
            source = events[name].add_source(
                srcname=reference) if reference else ''

        if bibcode not in itepbadsources:
            add_photometry(events, name, time=mjd, band=band, magnitude=magnitude,
                           e_magnitude=e_magnitude, source=sec_source + ',' + source)

    # Write out references that could use aa bibcode
    needsbib = list(OrderedDict.fromkeys(needsbib))
    with open('../itep-needsbib.txt', 'w') as bib_file:
        bib_file.writelines(['%ss\n' % ii for ii in needsbib])
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_pessto(events, stubs, args, tasks, task_obj, log):
    pessto_path = os.path.join(PATH.REPO_EXTERNAL, 'PESSTO_MPHOT.csv')
    tsvin = list(csv.reader(open(pessto_path, 'r'), delimiter=','))
    for ri, row in enumerate(tsvin):
        if ri == 0:
            bands = [xx.split('_')[0] for xx in row[3::2]]
            systems = [xx.split('_')[1].capitalize().replace(
                'Ab', 'AB') for xx in row[3::2]]
            continue
        name = row[1]
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2015A&A...579A..40S')
        events[name].add_quantity('alias', name, source)
        for hi, ci in enumerate(range(3, len(row) - 1, 2)):
            if not row[ci]:
                continue
            teles = 'Swift' if systems[hi] == 'Swift' else ''
            add_photometry(
                events, name, time=row[2], magnitude=row[
                    ci], e_magnitude=row[ci + 1],
                band=bands[hi], system=systems[hi], telescope=teles, source=source)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_scp(events, stubs, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    tsvin = list(csv.reader(
        open(os.path.join(PATH.REPO_EXTERNAL, 'SCP09.csv'), 'r'), delimiter=','))
    for ri, row in enumerate(pbar(tsvin, current_task)):
        if ri == 0:
            continue
        name = row[0].replace('SCP', 'SCP-')
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(srcname='Supernova Cosmology Project',
                                         url='http://supernova.lbl.gov/2009ClusterSurvey/')
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
            kind = ('spectroscopic/light curve' if 'a' in row[7] and 'c' in row[7] else
                    'spectroscopic' if 'a' in row[7] else
                    'light curve' if 'c' in row[7]
                    else '')
            if claimedtype != '?':
                events[name].add_quantity(
                    'claimedtype', claimedtype, source, kind=kind)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_sdss(events, stubs, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    with open(os.path.join(PATH.REPO_EXTERNAL, 'SDSS/2010ApJ...708..661D.txt'), 'r') as sdss_file:
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
                events, name = Events.add_event(tasks, args, events, name, log)
                source = events[name].add_source(bibcode=bibcode)
                events[name].add_quantity('alias', name, source)
                events[name].add_quantity(
                    'alias', 'SDSS-II SN ' + row[3], source)

                if row[5] != 'RA:':
                    year = re.findall(r'\d+', name)[0]
                    events[name].add_quantity('discoverdate', year, source)

                events[name].add_quantity(
                    'ra', row[-4], source, unit='floatdegrees')
                events[name].add_quantity(
                    'dec', row[-2], source, unit='floatdegrees')
            if rr == 1:
                error = row[4] if float(row[4]) >= 0.0 else ''
                events[name].add_quantity('redshift', row[2], source, error=error,
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
                add_photometry(
                    events, name, time=mjd, telescope=telescope, band=band, magnitude=magnitude,
                    e_magnitude=e_mag, source=source, system='SDSS')

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_snhunt(events, stubs, args, tasks, task_obj, log):
    current_task = task_obj.current_task(args)
    snh_url = 'http://nesssi.cacr.caltech.edu/catalina/current.html'
    html = load_cached_url(args, current_task, snh_url, os.path.join(
        PATH.REPO_EXTERNAL, 'SNhunt/current.html'))
    if not html:
        return events
    text = html.splitlines()
    findtable = False
    for ri, row in enumerate(text):
        if 'Supernova Discoveries' in row:
            findtable = True
        if findtable and '<table' in row:
            tstart = ri + 1
        if findtable and '</table>' in row:
            tend = ri - 1
    tablestr = '<html><body><table>'
    for row in text[tstart:tend]:
        if row[:3] == 'tr>':
            tablestr = tablestr + '<tr>' + row[3:]
        else:
            tablestr = tablestr + row
    tablestr = tablestr + '</table></body></html>'
    bs = BeautifulSoup(tablestr, 'html5lib')
    trs = bs.find('table').findAll('tr')
    for tr in pbar(trs, current_task):
        cols = [str(xx.text) for xx in tr.findAll('td')]
        if not cols:
            continue
        name = re.sub('<[^<]+?>', '', cols[4]
                      ).strip().replace(' ', '').replace('SNHunt', 'SNhunt')
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(srcname='Supernova Hunt', url=snh_url)
        events[name].add_quantity('alias', name, source)
        host = re.sub('<[^<]+?>', '', cols[1]).strip().replace('_', ' ')
        events[name].add_quantity('host', host, source)
        events[name].add_quantity('ra', cols[2], source, unit='floatdegrees')
        events[name].add_quantity('dec', cols[3], source, unit='floatdegrees')
        dd = cols[0]
        discoverdate = dd[:4] + '/' + dd[4:6] + '/' + dd[6:8]
        events[name].add_quantity('discoverdate', discoverdate, source)
        discoverers = cols[5].split('/')
        for discoverer in discoverers:
            events[name].add_quantity('discoverer', 'CRTS', source)
            events[name].add_quantity('discoverer', discoverer, source)
        if args.update:
            events, stubs = Events.journal_events(
                tasks, args, events, stubs, log)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_superfit_spectra(events, stubs, args, tasks, task_obj, log):
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
                events, stubs = Events.journal_events(
                    tasks, args, events, stubs, log)
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
                        list(filter(None, re.split('\tt+|\ss+', row, maxsplit=0))))
            specdata = [[xx.replace('D', 'E') for xx in list(ii)]
                        for ii in zip(*specdata)]
            wavelengths = specdata[0]
            fluxes = specdata[1]

            if epoff != '':
                mlmjd = astrotime(
                    '-'.join([str(mldt.year), str(mldt.month), str(mldt.day)])).mjd
                mlmjd = str(Decimal(mlmjd) + epoff)
            else:
                mlmjd = ''
            add_spectrum(
                events, name, 'Angstrom', 'Uncalibrated', u_time='MJD' if mlmjd else '', time=mlmjd,
                wavelengths=wavelengths, fluxes=fluxes, source=source)

            lastname = name

        events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_tns(events, stubs, args, tasks, task_obj, log):
    from datetime import timedelta
    session = requests.Session()
    current_task = task_obj.current_task(args)
    tns_url = 'https://wis-tns.weizmann.ac.il/'
    search_url = tns_url + \
        'search?&num_page=1&format=html&sort=desc&order=id&format=csv&page=0'
    csvtxt = load_cached_url(args, current_task, search_url, os.path.join(
        PATH.REPO_EXTERNAL, 'TNS/index.csv'))
    if not csvtxt:
        return events
    maxid = csvtxt.splitlines()[1].split(',')[0].strip('"')
    maxpages = ceil(int(maxid) / 1000.)

    for page in pbar(range(maxpages), current_task):
        fname = os.path.join(PATH.REPO_EXTERNAL, 'TNS/page-') + \
            str(page).zfill(2) + '.csv'
        if task_obj.load_archive(args) and os.path.isfile(fname) and page < 7:
            with open(fname, 'r') as tns_file:
                csvtxt = tns_file.read()
        else:
            with open(fname, 'w') as tns_file:
                session = requests.Session()
                ses_url = (tns_url + 'search?&num_page=1000&format=html&edit'
                           '[type]=&edit[objname]=&edit[id]=&sort=asc&order=id&display[redshift]=1'
                           '&display[hostname]=1&display[host_redshift]=1'
                           '&display[source_group_name]=1&display[programs_name]=1'
                           '&display[internal_name]=1&display[isTNS_AT]=1&display[public]=1'
                           '&display[end_pop_period]=0&display[spectra_count]=1'
                           '&display[discoverymag]=1&display[discmagfilter]=1'
                           '&display[discoverydate]=1&display[discoverer]=1&display[sources]=1'
                           '&display[bibcode]=1&format=csv&page=' + str(page))
                response = session.get(ses_url)
                csvtxt = response.text
                tns_file.write(csvtxt)

        tsvin = list(csv.reader(csvtxt.splitlines(), delimiter=','))
        for ri, row in enumerate(pbar(tsvin, current_task, leave=False)):
            if ri == 0:
                continue
            if row[4] and 'SN' not in row[4]:
                continue
            name = row[1].replace(' ', '')
            events, name = Events.add_event(tasks, args, events, name, log)
            source = events[name].add_source(
                srcname='Transient Name Server', url=tns_url)
            events[name].add_quantity('alias', name, source)
            if row[2] and row[2] != '00:00:00.00':
                events[name].add_quantity('ra', row[2], source)
            if row[3] and row[3] != '+00:00:00.00':
                events[name].add_quantity('dec', row[3], source)
            if row[4]:
                events[name].add_quantity(
                    'claimedtype', row[4].replace('SN', '').strip(), source)
            if row[5]:
                events[name].add_quantity(
                    'redshift', row[5], source, kind='spectroscopic')
            if row[6]:
                events[name].add_quantity('host', row[6], source)
            if row[7]:
                events[name].add_quantity(
                    'redshift', row[7], source, kind='host')
            if row[8]:
                events[name].add_quantity('discoverer', row[8], source)
            # Currently, all events listing all possible observers. TNS bug?
            # if row[9]:
            #    observers = row[9].split(',')
            #    for observer in observers:
            #        events[name].add_quantity('observer', observer.strip(), source)
            if row[10]:
                events[name].add_quantity('alias', row[10], source)
            if row[8] and row[14] and row[15] and row[16]:
                survey = row[8]
                magnitude = row[14]
                band = row[15].split('-')[0]
                mjd = astrotime(row[16]).mjd
                add_photometry(events, name, time=mjd, magnitude=magnitude, band=band,
                               survey=survey, source=source)
            if row[16]:
                date = row[16].split()[0].replace('-', '/')
                if date != '0000/00/00':
                    date = date.replace('/00', '')
                    time = row[16].split()[1]
                    if time != '00:00:00':
                        ts = time.split(':')
                        dt = timedelta(hours=int(ts[0]), minutes=int(
                            ts[1]), seconds=int(ts[2]))
                        date += pretty_num(dt.total_seconds() /
                                           (24 * 60 * 60), sig=6).lstrip('0')
                    events[name].add_quantity('discoverdate', date, source)
            if args.update:
                events, stubs = Events.journal_events(
                    tasks, args, events, stubs, log)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_simbad(events, stubs, args, tasks, task_obj, log):
    # Simbad.list_votable_fields()
    # Some coordinates that SIMBAD claims belong to the SNe actually belong to
    # the host.
    simbadmirrors = ['http://simbad.harvard.edu/simbad/sim-script',
                     'http://simbad.u-strasbg.fr/simbad/sim-script']
    simbadbadcoordbib = ['2013ApJ...770..107C']
    simbadbadnamebib = ['2004AJ....127.2809W', '2005MNRAS.364.1419Z',
                        '2015A&A...574A.112D', '2011MNRAS.417..916G', '2002ApJ...566..880G']
    simbadbannedcats = ['[TBV2008]', 'OGLE-MBR']
    customSimbad = Simbad()
    customSimbad.ROW_LIMIT = -1
    customSimbad.TIMEOUT = 120
    customSimbad.add_votable_fields('otype', 'sptype', 'sp_bibcode', 'id')
    for mirror in simbadmirrors:
        customSimbad.SIMBAD_URL = mirror
        try:
            table = customSimbad.query_criteria('maintype=SN | maintype="SN?"')
        except:
            continue
        else:
            break

    # 2000A&AS..143....9W
    for brow in pbar(table, current_task=current_task):
        row = {x: re.sub(r'b\'(.*)\'', r'\1',
                         str(brow[x])) for x in brow.colnames}
        # Skip items with no bibliographic info aside from SIMBAD, too
        # error-prone
        if row['OTYPE'] == 'Candidate_SN*' and not row['SP_TYPE']:
            continue
        if not row['COO_BIBCODE'] and not row['SP_BIBCODE'] and not row['SP_BIBCODE_2']:
            continue
        if any([x in row['MAIN_ID'] for x in simbadbannedcats]):
            continue
        if row['COO_BIBCODE'] and row['COO_BIBCODE'] in simbadbadnamebib:
            continue
        name = single_spaces(re.sub(r'\[[^)]*\]', '', row['MAIN_ID']).strip())
        if name == 'SN':
            continue
        if is_number(name):
            continue
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(srcname='SIMBAD astronomical database', bibcode="2000A&AS..143....9W",
                                         url="http://simbad.u-strasbg.fr/", secondary=True)
        aliases = row['ID'].split(',')
        for alias in aliases:
            if any([x in alias for x in simbadbannedcats]):
                continue
            ali = single_spaces(re.sub(r'\[[^)]*\]', '', alias).strip())
            if is_number(ali):
                continue
            ali = name_clean(ali)
            events[name].add_quantity('alias', ali, source)
        if row['COO_BIBCODE'] and row['COO_BIBCODE'] not in simbadbadcoordbib:
            csources = ','.join(
                [source, events[name].add_source(bibcode=row['COO_BIBCODE'])])
            events[name].add_quantity('ra', row['RA'], csources)
            events[name].add_quantity('dec', row['DEC'], csources)
        if row['SP_BIBCODE']:
            ssources = uniq_cdl([source, events[name].add_source(bibcode=row['SP_BIBCODE'])] +
                                ([events[name].add_source(bibcode=row['SP_BIBCODE_2'])] if row['SP_BIBCODE_2'] else []))
            events[name].add_quantity('claimedtype', row['SP_TYPE'].replace(
                'SN.', '').replace('SN', '').replace('(~)', '').strip(': '), ssources)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events
