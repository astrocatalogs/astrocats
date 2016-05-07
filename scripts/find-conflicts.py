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
from digits import *
from html import unescape
from glob import glob
from tqdm import tqdm
from collections import OrderedDict
from astropy.coordinates import SkyCoord as coord
from astropy import units as un
from astropy.time import Time as astrotime
from copy import deepcopy

conflicts = []

def get_event_filename(name):
    return(name.replace('/', '_'))

with open('rep-folders.txt', 'r') as f:
    repfolders = f.read().splitlines()

files = []
for rep in repfolders:
    files += glob('../' + rep + "/*.json") + glob('../' + rep + "/*.json.gz")

for fcnt, eventfile in enumerate(tqdm(sorted(files, key=lambda s: s.lower()))):
    #if fcnt > 100:
    #    break
    fileeventname = os.path.splitext(os.path.basename(eventfile))[0].replace('.json','')

    if eventfile.split('.')[-1] == 'gz':
        with gzip.open(eventfile, 'rt') as f:
            filetext = f.read()
    else:
        with open(eventfile, 'r') as f:
            filetext = f.read()

    item = json.loads(filetext, object_pairs_hook=OrderedDict)
    item = item[list(item.keys())[0]]

    ras = []
    decs = []
    zs = []
    rasources = []
    decsources = []
    zsources = []
    for key in list(item.keys()):
        lc = 0
        if key in ['name', 'sources', 'photometry', 'spectra']:
            continue
        if len(item[key]) == 1:
            continue
        for quantum in item[key]:
            if key == 'ra':
                newsources = []
                for alias in quantum['source'].split(','):
                    for source in item['sources']:
                        if source['alias'] == alias:
                            newsources.append({'idtype':'bibcode' if 'bibcode' in source else 'name',
                                'id':source['bibcode'] if 'bibcode' in source else source['name']})
                if newsources:
                    ras.append(quantum['value'])
                    rasources.append({'idtype':','.join([x['idtype'] for x in newsources]), 'id':','.join([x['id'] for x in newsources])})
            elif key == 'dec':
                newsources = []
                for alias in quantum['source'].split(','):
                    for source in item['sources']:
                        if source['alias'] == alias:
                            newsources.append({'idtype':'bibcode' if 'bibcode' in source else 'name',
                                'id':source['bibcode'] if 'bibcode' in source else source['name']})
                if newsources:
                    decs.append(quantum['value'])
                    decsources.append({'idtype':','.join([x['idtype'] for x in newsources]), 'id':','.join([x['id'] for x in newsources])})
            elif key == 'redshift':
                newsources = []
                for alias in quantum['source'].split(','):
                    for source in item['sources']:
                        if source['alias'] == alias:
                            newsources.append({'idtype':'bibcode' if 'bibcode' in source else 'name',
                                'id':source['bibcode'] if 'bibcode' in source else source['name']})
                if newsources:
                    zs.append(float(quantum['value']))
                    zsources.append({'idtype':','.join([x['idtype'] for x in newsources]), 'id':','.join([x['id'] for x in newsources])})
                
    edit = True if os.path.isfile('../sne-internal/' + get_event_filename(item['name']) + '.json') else False

    if ras and decs and item['name'] not in ['SN1996D', 'SN1998ew', 'SN2003an', 'SN2011in', 'SN2012ac', 'SN2012ht', 'SN2013bz']:
        oralen = len(ras)
        odeclen = len(decs)
        if len(ras) > len(decs):
            decs = decs + [decs[0] for x in range(len(ras) - len(decs))]
        elif len(ras) < len(decs):
            ras = ras + [ras[0] for x in range(len(decs) - len(ras))]

        coo = coord(ras, decs, unit = (un.hourangle, un.deg))
        ras = ras[:oralen]
        decs = decs[:odeclen]
        radegs = coo.ra.deg[:oralen]
        decdegs = coo.dec.deg[:odeclen]

        if len(radegs) > 1:
            maxradiff = max([abs((radegs[i+1]-radegs[i])/radegs[i+1]) for i in range(len(radegs)-1)])
            if maxradiff > 0.01:
                tqdm.write('R.A. difference greater than a % for ' + item['name'])
                conflicts.append(OrderedDict([('name', item['name']), ('alias', item['alias']), ('edit',edit),
                    ('quantity', 'ra'), ('difference', str(round_sig(maxradiff))), ('values', ras), ('sources', rasources)]))
        if len(decdegs) > 1:
            maxdecdiff = max([abs((decdegs[i+1]-decdegs[i])/decdegs[i+1]) for i in range(len(decdegs)-1)])
            if maxdecdiff > 0.01:
                tqdm.write('Dec. difference greater than a % for ' + item['name'])
                conflicts.append(OrderedDict([('name', item['name']), ('alias', item['alias']), ('edit',edit),
                    ('quantity', 'dec'), ('difference', str(round_sig(maxdecdiff))), ('values', decs), ('sources', decsources)]))

    if zs:
        maxzdiff = max([abs((zs[i+1]-zs[i])/zs[i+1]) for i in range(len(zs)-1)])
        if maxzdiff > 0.01:
            tqdm.write('Redshift difference greater than a % for ' + item['name'])
            conflicts.append(OrderedDict([('name', item['name']), ('alias', item['alias']), ('edit',edit),
                ('quantity', 'redshift'), ('difference', str(round_sig(maxzdiff))), ('values', zs), ('sources', zsources)]))

# Convert to array since that's what datatables expects
jsonstring = json.dumps(conflicts, indent='\t', separators=(',', ':'), ensure_ascii=False)
with open('../conflicts.json', 'w') as f:
    f.write(jsonstring)
