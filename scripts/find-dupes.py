#!/usr/local/bin/python3.5

import json
import re
import os
import math
import codecs
import gzip
from repos import *
from glob import glob
from tqdm import tqdm
from collections import OrderedDict
from astropy.coordinates import SkyCoord as coord
from astropy import units as un
from astropy.time import Time as astrotime
from copy import deepcopy

dupes = OrderedDict()

def get_event_filename(name):
    return(name.replace('/', '_'))

def get_event_filename(name):
    return(name.replace('/', '_'))

files = repo_file_list(bones = False)

newcatalog = []

for fcnt, eventfile in enumerate(tqdm(sorted(files, key=lambda s: s.lower()))):
    #if fcnt > 1000:
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
    newitem = OrderedDict()

    if 'maxdate' in item and item['maxdate']:
        date = item['maxdate'][0]['value'].replace('/', '-')
        negdate = date.startswith('-')
        datesplit = date.lstrip('-').split('-')
        if len(datesplit) >= 1:
            if '<' in datesplit[0]:
                print(item['name'])
                datesplit[0] = datesplit[0].strip('<')
            maxyear = float(datesplit[0])
        if len(datesplit) >= 2:
            maxyear += float(datesplit[1])/12.
        if len(datesplit) >= 3:
            maxyear += float(datesplit[2])/(12.*30.)
        if negdate:
            maxyear = -maxyear
        newitem['maxyear'] = maxyear
    if 'ra' in item and 'dec' in item and item['ra'] and item['dec']:
        newitem['name'] = item['name']
        newitem['alias'] = [x['value'] for x in item['alias']]
        newitem['ra'] = item['ra'][0]['value']
        newitem['dec'] = item['dec'][0]['value']
        # Temporary fix for David's typo
        #if newitem['dec'].count('.') == 2:
        #    newitem['dec'] = newitem['dec'][:newitem['dec'].rfind('.')]
        if 'distinctfrom' in item:
            newitem['distinctfrom'] = [x['value'] for x in item['distinctfrom']]
        newcatalog.append(newitem)

coo = coord([x['ra'] for x in newcatalog], [x['dec'] for x in newcatalog], unit = (un.hourangle, un.deg))
radegs = coo.ra.deg
decdegs = coo.dec.deg

for i, item in enumerate(newcatalog):
    newcatalog[i]['radeg'] = radegs[i]
    newcatalog[i]['decdeg'] = decdegs[i]

newcatalog2 = deepcopy(newcatalog)

for item1 in tqdm(newcatalog):
    name1 = item1['name']

    maxyear1 = None
    if 'maxyear' in item1 and item1['maxyear']:
        maxyear1 = item1['maxyear']

    for item2 in newcatalog2[:]:
        name2 = item2['name']
        if name1 == name2:
            newcatalog2.remove(item2)
            continue

        aliases1 = item1['alias']
        aliases2 = item2['alias']

        distinctfrom1 = item1['distinctfrom'] if 'distinctfrom' in item1 else []
        distinctfrom2 = item2['distinctfrom'] if 'distinctfrom' in item2 else []

        if (len(set(aliases1).intersection(distinctfrom2))):
            tqdm.write('Found ' + name2 + ' in distinct from list of ' + name1 + '.')
            continue
        if (len(set(aliases2).intersection(distinctfrom1))):
            tqdm.write('Found ' + name1 + ' in distinct from list of ' + name2 + '.')
            continue

        maxyear2 = None
        if 'maxyear' in item2 and item2['maxyear']:
            maxyear2 = item2['maxyear']

        ra1 = item1['ra']
        ra2 = item2['ra']
        dec1 = item1['dec']
        dec2 = item2['dec']
        radeg1 = item1['radeg']
        radeg2 = item2['radeg']
        decdeg1 = item1['decdeg']
        decdeg2 = item2['decdeg']

        diffyear = ''
        if radeg1 == radeg2 and decdeg1 == decdeg2:
            distdeg = 0.0
            tqdm.write(name1 + ' has an exact coordinate match to ' + name2)
        else:
            distdeg = math.hypot((radeg1 - radeg2), (decdeg1 - decdeg2))
            if distdeg < 10./3600.:
                if maxyear1 and maxyear2:
                    diffyear = abs(maxyear1 - maxyear2)
                    if diffyear <= 2.0:
                        tqdm.write(name1 + ' has a close coordinate and date match to ' + name2 + " [" + str(distdeg) + ', ' +
                              str(diffyear) + ']')
                    else:
                        tqdm.write(name1 + ' has a close coordinate, but significantly different date, to ' + name2 + " [" + str(distdeg) + ', ' +
                              str(diffyear) + ']')
                        #continue
                else:
                    tqdm.write(name1 + ' has a close coordinate match to ' + name2 + " [" + str(distdeg) + "]")
                if (not name1.startswith(('SN', 'AT')) and name2.startswith(('SN', 'AT')) or
                    (maxyear1 and maxyear2 and maxyear2 < maxyear1 and not name1.startswith(('SN', 'AT')))):
                    name1,name2 = name2,name1
                    aliases1,aliases2 = aliases2,aliases1
                    ra1,ra2 = ra2,ra1
                    dec1,dec2 = dec2,dec1
            else:
                continue

        edit = True if os.path.isfile('../sne-internal/' + get_event_filename(name1) + '.json') else False

        dupes[name1] = OrderedDict([('name1',name1), ('aliases1',aliases1), ('name2',name2), ('aliases2',aliases2), ('ra1',ra1), ('dec1',dec1),
            ('ra2',ra2), ('dec2',dec2), ('distdeg',str(distdeg)), ('diffyear',str(diffyear)), ('edit',edit)])

# Convert to array since that's what datatables expects
dupes = list(dupes.values())
jsonstring = json.dumps(dupes, indent='\t', separators=(',', ':'), ensure_ascii=False)
with open('../dupes.json', 'w') as f:
    f.write(jsonstring)
