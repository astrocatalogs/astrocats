#!/usr/local/bin/python3.5

import json
import os
import gzip
from glob import glob
from tqdm import tqdm
from collections import OrderedDict

def get_event_filename(name):
    return(name.replace('/', '_'))

with open('rep-folders.txt', 'r') as f:
    repfolders = f.read().splitlines()

files = []
for rep in repfolders:
    files += glob('../' + rep + "/*.json") + glob('../' + rep + "/*.json.gz")

spectracount = 0
photocount = 0
eventswithspectra = 0
eventswithphoto = 0

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

    if 'spectra' in item:
        eventswithspectra += 1
        spectracount += len(item['spectra'])

    if 'photometry' in item:
        eventswithphoto += 1
        photocount += len(item['photometry'])

print('Event count: ' + str(len(files)))
print('Events with spectra: ' + str(eventswithspectra))
print('Events with photometry: ' + str(eventswithphoto))
print('Total spectra: ' + str(spectracount))
print('Total photometry: ' + str(photocount))
