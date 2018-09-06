#!/usr/local/bin/python3.5

import argparse
import codecs
import json
import urllib
from collections import OrderedDict

from tqdm import tqdm

from astrocats.structures.catalog import Catalog

parser = argparse.ArgumentParser(
    description='Generate a sky location map AstroCats data.'
)
parser.add_argument(
    '--catalog',
    '-c',
    dest='catalog',
    help='Select which catalog to generate',
    default='sne',
    type=str)
args = parser.parse_args()

if args.catalog == 'tde':
    moduledir = 'tidaldisruptions'
    modulename = 'tde'
    moduleurl = 'tde.space'
    moduletitle = 'TDE'
elif args.catalog == 'sne':
    moduledir = 'supernovae'
    modulename = 'sne'
    moduleurl = 'sne.space'
    moduletitle = 'Supernova'
elif args.catalog == 'kne':
    moduledir = 'kilonovae'
    modulename = 'kne'
    moduleurl = 'kilonova.space'
    moduletitle = 'Kilonova'
else:
    raise ValueError('Unknown catalog!')

outdir = "astrocats/" + moduledir + "/output/html/"

atels = OrderedDict()
yearseps = range(1990, 2020, 2)
for i in tqdm(range(len(yearseps) - 1)):
    adsquery = (Catalog.ADS_BIB_URL + '*ATel.*' +
                '&data_type=Custom&nr_to_return=10000&format=%25R' +
                '&start_mon=&start_year=' +
                str(yearseps[i]) + '&end_mon=&end_year=' +
                str(yearseps[i + 1] + 1))
    response = urllib.request.urlopen(adsquery)
    newatels = list(
        filter(None, response.read().decode('utf-8').splitlines()))[2:]
    atels.update(OrderedDict(sorted(
        [(list(filter(None, x.split('ATel.')[-1].split('.')))[0],
          x.strip()) for x in newatels])))

jsonstring = json.dumps(atels, separators=(',', ':'), ensure_ascii=False)
path = outdir + '/atels.json'
with codecs.open(path, 'w', encoding='utf8') as f:
    f.write(jsonstring)

cbets = OrderedDict()
yearseps = range(2000, 2020, 2)
for i in tqdm(range(len(yearseps) - 1)):
    adsquery = (Catalog.ADS_BIB_URL + '*CBET.*' +
                '&data_type=Custom&nr_to_return=10000&format=%25R' +
                '&start_mon=&start_year=' +
                str(yearseps[i]) + '&end_mon=&end_year=' +
                str(yearseps[i + 1] + 1))
    response = urllib.request.urlopen(adsquery)
    newcbets = list(
        filter(None, response.read().decode('utf-8').splitlines()))[2:]
    cbets.update(OrderedDict(sorted(
        [(list(filter(None, x.split('CBET.')[-1].split('.')))[0],
          x.strip()) for x in newcbets])))

jsonstring = json.dumps(cbets, separators=(',', ':'), ensure_ascii=False)
path = outdir + '/cbets.json'
with codecs.open(path, 'w', encoding='utf8') as f:
    f.write(jsonstring)

iaucs = OrderedDict()
yearseps = range(1920, 2020, 2)
for i in tqdm(range(len(yearseps) - 1)):
    adsquery = (Catalog.ADS_BIB_URL + '*IAUC.*' +
                '&data_type=Custom&nr_to_return=10000&format=%25R' +
                '&start_mon=&start_year=' +
                str(yearseps[i]) + '&end_mon=&end_year=' +
                str(yearseps[i + 1] + 1))
    response = urllib.request.urlopen(adsquery)
    newiaucs = list(
        filter(None, response.read().decode('utf-8').splitlines()))[2:]
    iaucs.update(OrderedDict(sorted(
        [(list(filter(None, x.split('IAUC.')[-1].split('.')))[0],
          x.strip()) for x in newiaucs])))

jsonstring = json.dumps(iaucs, separators=(',', ':'), ensure_ascii=False)
path = outdir + 'iaucs.json'
with codecs.open(path, 'w', encoding='utf8') as f:
    f.write(jsonstring)
