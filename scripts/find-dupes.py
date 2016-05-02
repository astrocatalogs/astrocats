#!/usr/local/bin/python3.5

import json
import re
import os
import math
import codecs
from tqdm import tqdm
from collections import OrderedDict
from astropy.coordinates import SkyCoord as coord
from astropy import units as un
from astropy.time import Time as astrotime
from copy import deepcopy

dupes = OrderedDict()

def get_event_filename(name):
    return(name.replace('/', '_'))

with open('../catalog.min.json', 'r') as f:
    catalog = json.loads(f.read())
    newcatalog = []
    skipnames = []
    for i, item in enumerate(catalog):
        if 'maxdate' in item and item['maxdate']:
            date = item['maxdate'][0]['value'].replace('/', '-')
            datesplit = date.split('-')
            if len(datesplit) >= 1:
                if '<' in datesplit[0]:
                    print(item['name'])
                    datesplit[0] = datesplit[0].strip('<')
                maxyear = float(datesplit[0])
            if len(datesplit) >= 2:
                maxyear += float(datesplit[1])/12.
            if len(datesplit) >= 3:
                maxyear += float(datesplit[2])/(12.*30.)
            catalog[i]['maxdate'] = maxyear
        if 'ra' in item and 'dec' in item and item['ra'] and item['dec']:
            ra = item['ra'][0]['value']
            dec = item['dec'][0]['value']
            catalog[i]['ra'] = ra
            catalog[i]['dec'] = dec
            newcatalog.append(item)

    coo = coord([x['ra'] for x in newcatalog], [x['dec'] for x in newcatalog], unit = (un.hourangle, un.deg))
    radegs = coo.ra.deg
    decdegs = coo.dec.deg

    for i, item in enumerate(newcatalog):
        newcatalog[i]['radeg'] = radegs[i]
        newcatalog[i]['decdeg'] = decdegs[i]

    newcatalog2 = deepcopy(newcatalog)

    for item1 in tqdm(newcatalog):
        name1 = item1['name']
        aliases = item1['aliases']

        maxyear1 = None
        if 'maxdate' in item1 and item1['maxdate']:
            maxyear1 = item1['maxdate']

        for item2 in newcatalog2[:]:
            name2 = item2['name']
            if name1 == name2:
                newcatalog2.remove(item2)
                continue

            maxyear2 = None
            if 'maxdate' in item2 and item2['maxdate']:
                maxyear2 = item2['maxdate']

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
                    if (not name1.startswith('SN') and name2.startswith('SN') or
                        (maxyear1 and maxyear2 and maxyear2 < maxyear1 and not name1.startswith('SN'))):
                        name1,name2 = name2,name1
                        ra1,ra2 = ra2,ra1
                        dec1,dec2 = dec2,dec1
                else:
                    continue

            edit = True if os.path.isfile('../sne-internal/' + get_event_filename(name1) + '.json') else False

            aliases1 = item1['aliases']
            aliases2 = item2['aliases']
            dupes[name1] = OrderedDict([('name1',name1), ('aliases1',aliases1), ('name2',name2), ('aliases2',aliases2), ('ra1',ra1), ('dec1',dec1),
                ('ra2',ra2), ('dec2',dec2), ('distdeg',str(distdeg)), ('diffyear',str(diffyear)), ('edit',edit)])

# Convert to array since that's what datatables expects
dupes = list(dupes.values())
jsonstring = json.dumps(dupes, indent='\t', separators=(',', ':'), ensure_ascii=False)
with open('../dupes.json', 'w') as f:
    f.write(jsonstring)
