#!/usr/local/bin/python3.5

import json
import re
import math
from astropy.time import Time as astrotime

with open('../catalog.min.json', 'r') as f:
    filetext = f.read()
    catalog = json.loads(filetext)
    for i, item in enumerate(catalog):
        catalog[i]['name'] = re.sub('<[^<]+?>', '', item['name']) 
        if 'discoverdate' in item and item['discoverdate']:
            date = item['discoverdate'][0]['value'].replace('/', '-')
            datesplit = date.split('-')
            if len(datesplit) >= 1:
                discoveryear = float(datesplit[0])
            if len(datesplit) >= 2:
                discoveryear += float(datesplit[1])/12.
            if len(datesplit) >= 3:
                discoveryear += float(datesplit[2])/(12.*30.)
            catalog[i]['discoverdate'] = discoveryear

    catalog1 = catalog
    catalog2 = catalog

    for item1 in catalog1:
        name1 = item1['name']
        aliases = item1['aliases']
        if 'ra' not in item1 or 'dec' not in item1 or not item1['ra'] or not item1['dec']:
            continue
        ra1 = item1['ra'][0]['value']
        dec1 = item1['dec'][0]['value']

        if 'discoverdate' in item1 and item1['discoverdate']:
            discoveryear1 = item1['discoverdate']

        radeg1 = 0.0
        rasp1 = ra1.split(':')
        if len(rasp1) >= 1:
            radeg1 = float(rasp1[0].strip())*360./24.
        if len(rasp1) >= 2:
            radeg1 += float(rasp1[1].strip())*360./24./60.
        if len(rasp1) >= 3:
            radeg1 += float(rasp1[2].strip())*360./24./60./60.

        decdeg1 = 0.0
        decsp1 = dec1.split(':')
        sign = 1.0
        if len(decsp1) >= 1:
            if decsp1[0][0] == '-':
                sign = -1.0
            decdeg1 = float(decsp1[0].strip('- '))
        if len(decsp1) >= 2:
            decdeg1 += float(decsp1[1].strip())/60.
        if len(decsp1) >= 3:
            decdeg1 += float(decsp1[2].strip())/60./60.
        decdeg1 *= sign

        for item2 in catalog2[:]:
            name2 = item2['name']
            if name2 in aliases and name1 != name2:
                print (name2 + ' in alias list of ' + name1)
            if name1 == name2 or 'ra' not in item2 or 'dec' not in item2 or not item2['ra'] or not item2['dec']:
                catalog2.remove(item2)
                continue
            discoveryear2 = None
            if 'discoverdate' in item2 and item2['discoverdate']:
                discoveryear2 = item2['discoverdate']

            ra2 = item2['ra'][0]['value']
            dec2 = item2['dec'][0]['value']
            radeg2 = 0.0
            rasp2 = ra2.split(':')
            if len(rasp2) >= 1:
                radeg2 = float(rasp2[0].strip())*360./24.
            if len(rasp2) >= 2:
                radeg2 += float(rasp2[1].strip())*360./24./60.
            if len(rasp2) >= 3:
                radeg2 += float(rasp2[2].strip())*360./24./60./60.

            decdeg2 = 0.0
            decsp2 = dec2.split(':')
            sign = 1.0
            if len(decsp2) >= 1:
                if decsp2[0][0] == '-':
                    sign = -1.0
                decdeg2 = float(decsp2[0].strip('- '))
            if len(decsp2) >= 2:
                decdeg2 += float(decsp2[1].strip())/60.
            if len(decsp2) >= 3:
                decdeg2 += float(decsp2[2].strip())/60./60.
            decdeg2 *= sign

            if ra1 == ra2 and dec1 == dec2:
                print(name1 + ' has an exact coordinate match to ' + name2)
                continue

            distdeg = math.sqrt(pow(radeg1 - radeg2, 2) + pow(decdeg1 - decdeg2, 2))
            if distdeg < 2./3600.:
                if discoveryear1 and discoveryear2:
                    if abs(discoveryear1 - discoveryear2) <= 0.5:
                        print(name1 + ' has a close coordinate and date match to ' + name2 + " [" + str(distdeg) + ', ' +
                              str(abs(discoveryear1) - abs(discoveryear2)) + ']')
                else:
                    print(name1 + ' has a close coordinate match to ' + name2 + " [" + str(distdeg) + "]")
