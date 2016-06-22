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
from html import unescape
from glob import glob
from tqdm import tqdm
from collections import OrderedDict
from astropy.coordinates import SkyCoord as coord
from astropy import units as un
from astropy.time import Time as astrotime
from copy import deepcopy
from events import *
from repos import *

errata = []

files = repo_file_list(bones = False)

for fcnt, eventfile in enumerate(tqdm(sorted(files, key=lambda s: s.lower()))):
    #if fcnt > 100:
    #    break
    fileeventname = os.path.splitext(os.path.basename(eventfile))[0].replace('.json','')

    filetext = get_event_text(eventfile)

    item = json.loads(filetext, object_pairs_hook=OrderedDict)
    item = item[list(item.keys())[0]]

    if 'errors' in item:
        for error in item['errors']:
            quantity = error['extra']
            likelyvalue = ''
            if quantity in list(item.keys()) and 'value' in item[quantity][0]:
                likelyvalue = item[quantity][0]['value']
            errata.append(OrderedDict([('name', item['name']), ('alias', item['alias']), ('ident', error['value']), ('kind', error['kind']),
                ('quantity', error['extra']), ('likelyvalue', likelyvalue)]))

jsonstring = json.dumps(errata, indent='\t', separators=(',', ':'), ensure_ascii=False)
with open('../errata.json', 'w') as f:
    f.write(jsonstring)
