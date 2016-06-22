#!/usr/local/bin/python3.5

import json
import re
import os
import math
import codecs
import urllib
import requests
import ads
import gzip
import statistics
import time
from html import unescape
from glob import glob
from tqdm import tqdm
from collections import OrderedDict
from astropy.coordinates import SkyCoord as coord
from astropy import units as un
from astropy.time import Time as astrotime
from copy import deepcopy
from math import sqrt
from digits import *

hosts = OrderedDict()

def get_event_filename(name):
    return(name.replace('/', '_'))

def touch(fname, times=None):
    with open(fname, 'a'):
        os.utime(fname, times)

with open('rep-folders.txt', 'r') as f:
    repfolders = f.read().splitlines()

files = []
for rep in repfolders:
    files += glob('../' + rep + "/*.json") + glob('../' + rep + "/*.json.gz")

for fcnt, eventfile in enumerate(tqdm(sorted(files, key=lambda s: s.lower()))):
    #if fcnt > 1000:
    #    break
    fileeventname = os.path.splitext(os.path.basename(eventfile))[0].replace('.json','')

    if not os.path.isfile(eventfile):
        continue

    if eventfile.split('.')[-1] == 'gz':
        with gzip.open(eventfile, 'rt') as f:
            filetext = f.read()
    else:
        with open(eventfile, 'r') as f:
            filetext = f.read()

    item = json.loads(filetext, object_pairs_hook=OrderedDict)
    item = item[list(item.keys())[0]]

    if 'host' in item:
        hngs = [x['value'] for x in item['host'] if ((x['kind'] != 'cluster') if 'kind' in x else True)]
        hncs = [x['value'] for x in item['host'] if ((x['kind'] == 'cluster') if 'kind' in x else False)]
        hng = ''
        hnc = ''
        for ho in hosts:
            hog = [x for x in hosts[ho]['host'] if hosts[ho]['kind'] != 'cluster']
            hoc = [x for x in hosts[ho]['host']]
            if len(list(set(hngs).intersection(hog))):
                hng = ho
                hosts[ho]['host'] = list(set(hosts[ho]['host'] + hngs))
            if len(list(set(hncs).intersection(hoc))):
                hnc = ho
                hosts[ho]['host'] = list(set(hosts[ho]['host'] + hncs + hngs))
            if hng and hnc:
                break

        if not hng and hngs:
            hng = hngs[0]
            hosts[hng] = OrderedDict([('host', hngs), ('kind', 'galaxy'), ('events', []), ('eventdates', []),
                ('types', []), ('photocount', 0), ('spectracount', 0), ('lumdist', ''),
                ('redshift', ''), ('hostra', ''), ('hostdec', '')])

        if not hnc and hncs:
            hnc = hncs[0]
            hosts[hnc] = OrderedDict([('host', hncs + hngs), ('kind', 'cluster'), ('events', []), ('eventdates', []),
                ('types', []), ('photocount', 0), ('spectracount', 0), ('lumdist', ''),
                ('redshift', ''), ('hostra', ''), ('hostdec', '')])

        for hi, hn in enumerate([hng, hnc]):
            if not hn:
                continue
            hosts[hn]['events'].append({'name':item['name'],'img':('ra' in item and 'dec' in item)})

            if (not hosts[hn]['lumdist'] or '*' in hosts[hn]['lumdist']) and 'lumdist' in item:
                ldkinds = [x['kind'] if 'kind' in x else '' for x in item['lumdist']]
                try:
                    ind = ldkinds.index('host')
                except ValueError:
                    hosts[hn]['lumdist'] = item['lumdist'][0]['value'] + '*'
                else:
                    hosts[hn]['lumdist'] = item['lumdist'][ind]['value']

            if (not hosts[hn]['redshift'] or '*' in hosts[hn]['redshift']) and 'redshift' in item:
                zkinds = [x['kind'] if 'kind' in x else '' for x in item['redshift']]
                try:
                    ind = zkinds.index('host')
                except ValueError:
                    hosts[hn]['redshift'] = item['redshift'][0]['value'] + '*'
                else:
                    hosts[hn]['redshift'] = item['redshift'][ind]['value']

            if not hosts[hn]['hostra'] and 'hostra' in item:
                hosts[hn]['hostra'] = item['hostra'][0]['value']
            if not hosts[hn]['hostdec'] and 'hostdec' in item:
                hosts[hn]['hostdec'] = item['hostdec'][0]['value']

            if 'discoverdate' in item and item['discoverdate']:
                datestr = item['discoverdate'][0]['value'].replace('/', '-')
                if datestr.count('-') == 1:
                    datestr += '-01'
                elif datestr.count('-') == 0:
                    datestr += '-01-01'
                try:
                    hosts[hn]['eventdates'].append(astrotime(datestr, format = 'isot').unix)
                except:
                    hosts[hn]['eventdates'].append(float("inf"))
            else:
                hosts[hn]['eventdates'].append(float("inf"))

            if 'claimedtype' in item:
                cts = []
                for ct in item['claimedtype']:
                    sct = ct['value'].strip('?')
                    if sct:
                        cts.append(sct)
                hosts[hn]['types'] = list(set(hosts[hn]['types']).union(cts))

            if 'photometry' in item:
                hosts[hn]['photocount'] += len(item['photometry'])

            if 'spectra' in item:
                hosts[hn]['spectracount'] += len(item['spectra'])

curtime = time.time()
centrate = 100.0*365.25*24.0*60.0*60.0

for hn in hosts:
    finitedates = sorted([x for x in hosts[hn]['eventdates'] if x != float("inf")])
    if len(finitedates) >= 2:
        datediff = curtime - finitedates[0]
        lamb = float(len(finitedates))/(curtime - finitedates[0])*centrate
        hosts[hn]['rate'] = (pretty_num(lamb, sig = 3) + ',' +
            pretty_num(lamb/sqrt(float(len(finitedates))), sig = 3))
    else:
        hosts[hn]['rate'] = ''
    hosts[hn]['events'] = [x for (y,x) in sorted(zip(hosts[hn]['eventdates'], hosts[hn]['events']), key = lambda ev: ev[0])]
    del(hosts[hn]['eventdates'])

# Convert to array since that's what datatables expects
hosts = list(hosts.values())

jsonstring = json.dumps(hosts, indent='\t', separators=(',', ':'), ensure_ascii=False)
with open('../hosts.json', 'w') as f:
    f.write(jsonstring)

minjsonstring = json.dumps(hosts, separators=(',', ':'), ensure_ascii=False)
with gzip.open("../hosts.min.json.gz", 'wt') as fff:
    touch("../hosts.min.json")
    fff.write(minjsonstring)
