#!/usr/local/bin/python3.5

import json
import re
import os
import math
import codecs
from glob import glob
from tqdm import tqdm
from collections import OrderedDict
from astropy.coordinates import SkyCoord as coord
from astropy import units as un
from astropy.time import Time as astrotime
from copy import deepcopy

biblio = OrderedDict()

def get_event_filename(name):
    return(name.replace('/', '_'))

with open('rep-folders.txt', 'r') as f:
    repfolders = f.read().splitlines()

files = []
for rep in repfolders:
    files += glob('../' + rep + "/*.json") + glob('../' + rep + "/*.json.gz")

for fcnt, eventfile in enumerate(sorted(files, key=lambda s: s.lower())):
    fileeventname = os.path.splitext(os.path.basename(eventfile))[0].replace('.json','')

    if eventfile.split('.')[-1] == 'gz':
        with gzip.open(eventfile, 'rt') as f:
            filetext = f.read()
    else:
        with open(eventfile, 'r') as f:
            filetext = f.read()

    item = json.loads(filetext, object_pairs_hook=OrderedDict)
    item = item[list(item.keys())[0]]

    if 'sources' in item:
        for source in item['sources']:
            if 'bibcode' in source:
                bc = source['bibcode']
                if bc in biblio:
                    biblio[bc]['events'].append(item['name'])
                else:
                    print(bc)
                    biblio[bc] = OrderedDict([('bibcode', bc), ('events', [item['name']])])

# Convert to array since that's what datatables expects
biblio = list(biblio.values())
jsonstring = json.dumps(biblio, indent='\t', separators=(',', ':'), ensure_ascii=False)
with open('../biblio.json', 'w') as f:
    f.write(jsonstring)
