"""Utility functions for OSC import.
"""

import json
import os
from collections import OrderedDict
from math import floor, hypot, log10, pi, sqrt

from astropy import units as un
from astropy.coordinates import SkyCoord as coord
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value

from cdecimal import Decimal

from ..utils import (get_sig_digits, is_number, pretty_num, round_sig, tprint,
                     zpad)
from .constants import (ADS_BIB_URL, CLIGHT, KM, OSC_BIBCODE, OSC_NAME,
                        OSC_URL, PREF_KINDS)


def convert_aq_output(row):
    return OrderedDict([(x, str(row[x]) if is_number(row[x]) else row[x])
                        for x in row.colnames])


def read_json_dict(filename):
    # path = '../atels.json'
    if os.path.isfile(filename):
        with open(filename, 'r') as f:
            mydict = json.loads(f.read(), object_pairs_hook=OrderedDict)
    else:
        mydict = OrderedDict()
    return mydict


def read_json_arr(filename):
    if os.path.isfile(filename):
        with open(filename, 'r') as f:
            myarr = json.loads(f.read())
    else:
        myarr = []
    return myarr


def get_preferred_name(entries, name):
    if name not in entries:
        # matches = []
        for event in entries:
            aliases = entries[event].get_aliases()
            if len(aliases) > 1 and name in aliases:
                return event
        return name
    else:
        return name


def get_source_year(source):
    if 'bibcode' in source:
        if is_number(source['bibcode'][:4]):
            return int(source['bibcode'][:4])
        else:
            return -10000
    raise ValueError('No bibcode available for source!')


def has_task(tasks, args, task):
    return task in tasks and (not args.update or tasks[task]['update'])


def jd_to_mjd(jd):
    return jd - Decimal(2400000.5)


def load_cached_url(args, current_task, url, filepath, timeout=120, write=True,
                    failhard=False):
    import codecs
    from hashlib import md5
    filemd5 = ''
    filetxt = ''
    if not args.refresh and os.path.isfile(filepath):
        with codecs.open(filepath, 'r', encoding='utf8') as f:
            filetxt = f.read()
            if args.update:
                filemd5 = md5(filetxt.encode('utf-8')).hexdigest()

    try:
        import requests
        session = requests.Session()
        response = session.get(url, timeout=timeout)
        response.raise_for_status()
        for x in response.history:
            x.raise_for_status()
            if (x.status_code == 500 or x.status_code == 307 or
                    x.status_code == 404):
                raise
        txt = response.text
        newmd5 = md5(txt.encode('utf-8')).hexdigest()
        # tprint(filemd5 + ": " + newmd5)
        if args.update and newmd5 == filemd5:
            tprint('Skipping file in "' + current_task +
                   '," local and remote copies identical [' + newmd5 + '].')
            return False
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        if failhard:
            return ''
        return filetxt
    else:
        if write:
            with codecs.open(filepath, 'w', encoding='utf8') as f:
                f.write(txt if txt else filetxt)
    return txt


def make_date_string(year, month='', day=''):
    if not year:
        raise ValueError(
            "At least the year must be specified when constructing date "
            "string")
    datestring = str(year)
    if month:
        datestring = datestring + '/' + str(month).zfill(2)
    if day:
        datestring = datestring + '/' + str(day).zfill(2)

    return datestring


def name_clean(name):
    newname = name.strip(' ;,*')
    if newname.startswith('NAME '):
        newname = newname.replace('NAME ', '', 1)
    if newname.endswith(' SN'):
        newname = newname.replace(' SN', '')
    if newname.endswith(':SN'):
        newname = newname.replace(':SN', '')
    if newname.startswith('MASJ'):
        newname = newname.replace('MASJ', 'MASTER OT J', 1)
    if newname.startswith('MASTER') and is_number(newname[7]):
        newname = newname.replace('MASTER', 'MASTER OT J', 1)
    if newname.startswith('MASTER OT J '):
        newname = newname.replace('MASTER OT J ', 'MASTER OT J', 1)
    if newname.startswith('OGLE '):
        newname = newname.replace('OGLE ', 'OGLE-', 1)
    if newname.startswith('OGLE-') and len(newname) != 16:
        namesp = newname.split('-')
        if (len(namesp[1]) == 4 and is_number(namesp[1]) and
                is_number(namesp[3])):
            newname = 'OGLE-' + namesp[1] + '-SN-' + namesp[3].zfill(3)
    if newname.startswith('SN SDSS'):
        newname = newname.replace('SN SDSS ', 'SDSS', 1)
    if newname.startswith('SDSS '):
        newname = newname.replace('SDSS ', 'SDSS', 1)
    if newname.startswith('SDSS'):
        namesp = newname.split('-')
        if (len(namesp) == 3 and is_number(namesp[0][4:]) and
                is_number(namesp[1]) and is_number(namesp[2])):
            newname = namesp[0] + '-' + namesp[1] + '-' + namesp[2].zfill(3)
    if newname.startswith('SDSS-II SN'):
        namesp = newname.split()
        if len(namesp) == 3 and is_number(namesp[2]):
            newname = 'SDSS-II SN ' + namesp[2].lstrip('0')
    if newname.startswith('SN CL'):
        newname = newname.replace('SN CL', 'CL', 1)
    if newname.startswith('SN HiTS '):
        newname = newname.replace('SN HiTS ', 'SNHiTS', 1)
    if newname.startswith('GAIA'):
        newname = newname.replace('GAIA', 'Gaia', 1)
    if newname.startswith('Gaia '):
        newname = newname.replace('Gaia ', 'Gaia', 1)
    if newname.startswith('Gaia'):
        newname = 'Gaia' + newname[4:].lower()
    if newname.startswith('GRB'):
        newname = newname.replace('GRB', 'GRB ', 1)
    if newname.startswith('GRB ') and is_number(newname[4:].strip()):
        newname = 'GRB ' + newname[4:].strip() + 'A'
    if newname.startswith('LSQ '):
        newname = newname.replace('LSQ ', 'LSQ', 1)
    if newname.startswith('KSN '):
        newname = newname.replace('KSN ', 'KSN-', 1)
    if newname.startswith('SNSDF '):
        newname = newname.replace(' ', '')
    if newname.startswith('SNSDF'):
        namesp = newname.split('.')
        if len(namesp[0]) == 9:
            newname = namesp[0] + '-' + namesp[1].zfill(2)
    if newname.startswith('HFF '):
        newname = newname.replace(' ', '')
    if newname.startswith('SN HST'):
        newname = newname.replace('SN HST', 'HST', 1)
    if newname.startswith('HST ') and newname[4] != 'J':
        newname = newname.replace('HST ', 'HST J', 1)
    if newname.startswith('SNLS') and newname[4] != '-':
        newname = newname.replace('SNLS', 'SNLS-', 1)
    if newname.startswith('SNLS- '):
        newname = newname.replace('SNLS- ', 'SNLS-', 1)
    if newname.startswith('CRTS CSS'):
        newname = newname.replace('CRTS CSS', 'CSS', 1)
    if newname.startswith('CRTS MLS'):
        newname = newname.replace('CRTS MLS', 'MLS', 1)
    if newname.startswith('CRTS SSS'):
        newname = newname.replace('CRTS SSS', 'SSS', 1)
    if newname.startswith(('CSS', 'MLS', 'SSS')):
        newname = newname.replace(' ', ':').replace('J', '')
    if newname.startswith('SN HFF'):
        newname = newname.replace('SN HFF', 'HFF', 1)
    if newname.startswith('SN GND'):
        newname = newname.replace('SN GND', 'GND', 1)
    if newname.startswith('SN SCP'):
        newname = newname.replace('SN SCP', 'SCP', 1)
    if newname.startswith('SN UDS'):
        newname = newname.replace('SN UDS', 'UDS', 1)
    if newname.startswith('SCP') and newname[3] != '-':
        newname = newname.replace('SCP', 'SCP-', 1)
    if newname.startswith('SCP- '):
        newname = newname.replace('SCP- ', 'SCP-', 1)
    if newname.startswith('PS 1'):
        newname = newname.replace('PS 1', 'PS1', 1)
    if newname.startswith('PS1 SN PS'):
        newname = newname.replace('PS1 SN PS', 'PS', 1)
    if newname.startswith('PS1 SN'):
        newname = newname.replace('PS1 SN', 'PS1', 1)
    if newname.startswith('PSN K'):
        newname = newname.replace('PSN K', 'K', 1)
    if newname.startswith('K') and is_number(newname[1:5]):
        namesp = newname.split('-')
        if len(namesp[0]) == 5:
            newname = namesp[0] + '-' + namesp[1].zfill(3)
    if newname.startswith('Psn'):
        newname = newname.replace('Psn', 'PSN', 1)
    if newname.startswith('PSNJ'):
        newname = newname.replace('PSNJ', 'PSN J', 1)
    if newname.startswith('TCPJ'):
        newname = newname.replace('TCPJ', 'TCP J', 1)
    if newname.startswith('SMTJ'):
        newname = newname.replace('SMTJ', 'SMT J', 1)
    if newname.startswith('PSN20J'):
        newname = newname.replace('PSN20J', 'PSN J', 1)
    if newname.startswith('SN ASASSN'):
        newname = newname.replace('SN ASASSN', 'ASASSN', 1)
    if newname.startswith('ASASSN '):
        newname = newname.replace('ASASSN ', 'ASASSN-', 1).replace('--', '-')
    if newname.startswith('ASASSN') and newname[6] != '-':
        newname = newname.replace('ASASSN', 'ASASSN-', 1)
    if newname.startswith('ROTSE3J'):
        newname = newname.replace('ROTSE3J', 'ROTSE3 J', 1)
    if newname.startswith('MACSJ'):
        newname = newname.replace('MACSJ', 'MACS J', 1)
    if newname.startswith('MWSNR'):
        newname = newname.replace('MWSNR', 'MWSNR ', 1)
    if newname.startswith('SN HUNT'):
        newname = newname.replace('SN HUNT', 'SNhunt', 1)
    if newname.startswith('SN Hunt'):
        newname = newname.replace(' ', '')
    if newname.startswith('SNHunt'):
        newname = newname.replace('SNHunt', 'SNhunt', 1)
    if newname.startswith('SNhunt '):
        newname = newname.replace('SNhunt ', 'SNhunt', 1)
    if newname.startswith('ptf'):
        newname = newname.replace('ptf', 'PTF', 1)
    if newname.startswith('SN PTF'):
        newname = newname.replace('SN PTF', 'PTF', 1)
    if newname.startswith('PTF '):
        newname = newname.replace('PTF ', 'PTF', 1)
    if newname.startswith('iPTF '):
        newname = newname.replace('iPTF ', 'iPTF', 1)
    if newname.startswith('PESSTOESO'):
        newname = newname.replace('PESSTOESO', 'PESSTO ESO ', 1)
    if newname.startswith('snf'):
        newname = newname.replace('snf', 'SNF', 1)
    if newname.startswith('SNF '):
        newname = newname.replace('SNF ', 'SNF', 1)
    if (newname.startswith('SNF') and
            is_number(newname[3:]) and len(newname) >= 12):
        newname = 'SNF' + newname[3:11] + '-' + newname[11:]
    if newname.startswith(('MASTER OT J', 'ROTSE3 J')):
        prefix = newname.split('J')[0]
        coords = newname.split('J')[-1].strip()
        decsign = '+' if '+' in coords else '-'
        coordsplit = coords.replace('+', '-').split('-')
        if ('.' not in coordsplit[0] and
                len(coordsplit[0]) > 6 and '.' not in coordsplit[1] and
                len(coordsplit[1]) > 6):
            newname = (prefix + 'J' + coordsplit[0][:6] + '.' +
                       coordsplit[0][6:] + decsign + coordsplit[1][:6] +
                       '.' + coordsplit[1][6:])
    if (newname.startswith('Gaia ') and
            is_number(newname[3:4]) and len(newname) > 5):
        newname = newname.replace('Gaia ', 'Gaia', 1)
    if len(newname) <= 4 and is_number(newname):
        newname = 'SN' + newname + 'A'
    if (len(newname) > 4 and is_number(newname[:4]) and not
            is_number(newname[4:])):
        newname = 'SN' + newname
    if (newname.startswith('Sn ') and
            is_number(newname[3:7]) and len(newname) > 7):
        newname = newname.replace('Sn ', 'SN', 1)
    if (newname.startswith('sn') and
            is_number(newname[2:6]) and len(newname) > 6):
        newname = newname.replace('sn', 'SN', 1)
    if (newname.startswith('SN ') and
            is_number(newname[3:7]) and len(newname) > 7):
        newname = newname.replace('SN ', 'SN', 1)
    if (newname.startswith('SN') and
            is_number(newname[2:6]) and len(newname) == 7 and
            newname[6].islower()):
        newname = 'SN' + newname[2:6] + newname[6].upper()
    elif (newname.startswith('SN') and is_number(newname[2:6]) and
          (len(newname) == 8 or len(newname) == 9) and newname[6:].isupper()):
        newname = 'SN' + newname[2:6] + newname[6:].lower()

    newname = (' '.join(newname.split())).strip()
    return newname


def radec_clean(svalue, quantity, unit=''):
    if unit == 'floatdegrees':
        if not is_number(svalue):
            return (svalue, unit)
        deg = float('%g' % Decimal(svalue))
        sig = get_sig_digits(svalue)
        if 'ra' in quantity:
            flhours = deg / 360.0 * 24.0
            hours = floor(flhours)
            minutes = floor((flhours - hours) * 60.0)
            seconds = (flhours * 60.0 - (hours * 60.0 + minutes)) * 60.0
            hours = 0 if hours < 1.e-6 else hours
            minutes = 0 if minutes < 1.e-6 else minutes
            seconds = 0.0 if seconds < 1.e-6 else seconds
            if seconds > 60.0:
                raise(ValueError('Invalid seconds value for ' + quantity))
            svalue = str(hours).zfill(2) + ':' + str(minutes).zfill(2) + \
                ':' + zpad(pretty_num(seconds, sig=sig - 1))
        elif 'dec' in quantity:
            fldeg = abs(deg)
            degree = floor(fldeg)
            minutes = floor((fldeg - degree) * 60.0)
            seconds = (fldeg * 60.0 - (degree * 60.0 + minutes)) * 60.0
            if seconds > 60.0:
                raise(ValueError('Invalid seconds value for ' + quantity))
            svalue = (('+' if deg >= 0.0 else '-') +
                      str(degree).strip('+-').zfill(2) + ':' +
                      str(minutes).zfill(2) + ':' +
                      zpad(pretty_num(seconds, sig=sig - 1)))
    elif unit == 'nospace' and 'ra' in quantity:
        svalue = svalue[:2] + ':' + svalue[2:4] + \
            ((':' + zpad(svalue[4:])) if len(svalue) > 4 else '')
    elif unit == 'nospace' and 'dec' in quantity:
        if svalue.startswith(('+', '-')):
            svalue = svalue[:3] + ':' + svalue[3:5] + \
                ((':' + zpad(svalue[5:])) if len(svalue) > 5 else '')
        else:
            svalue = '+' + svalue[:2] + ':' + svalue[2:4] + \
                ((':' + zpad(svalue[4:])) if len(svalue) > 4 else '')
    else:
        svalue = svalue.replace(' ', ':')
        if 'dec' in quantity:
            valuesplit = svalue.split(':')
            svalue = (('-' if valuesplit[0].startswith('-') else '+') +
                      valuesplit[0].strip('+-').zfill(2) +
                      (':' + valuesplit[1].zfill(2) if
                       len(valuesplit) > 1 else '') +
                      (':' + zpad(valuesplit[2]) if
                       len(valuesplit) > 2 else ''))

    if 'ra' in quantity:
        sunit = 'hours'
    elif 'dec' in quantity:
        sunit = 'degrees'

    # Correct case of arcseconds = 60.0.
    valuesplit = svalue.split(':')
    if len(valuesplit) == 3 and valuesplit[-1] in ["60.0", "60.", "60"]:
        svalue = valuesplit[0] + ':' + str(Decimal(valuesplit[1]) +
                                           Decimal(1.0)) + ':' + "00.0"

    # Strip trailing dots.
    svalue = svalue.rstrip('.')

    return (svalue, sunit)


def host_clean(name):
    newname = name.strip(' ;,*')

    # Handle some special cases
    hostcases = {'M051a': 'M51A', 'M051b': 'M51B'}
    for k in hostcases:
        if newname == k:
            newname = hostcases[k]

    # Some general cases
    newname = newname.strip("()").replace('  ', ' ', 1)
    newname = newname.replace("ABELL", "Abell", 1)
    newname = newname.replace("Abell", "Abell ", 1)
    newname = newname.replace("APMUKS(BJ)", "APMUKS(BJ) ", 1)
    newname = newname.replace("ARP", "ARP ", 1)
    newname = newname.replace("CGCG", "CGCG ", 1)
    newname = newname.replace("HOLM", "HOLM ", 1)
    newname = newname.replace("IC", "IC ", 1)
    newname = newname.replace("Intergal.", "Intergalactic", 1)
    newname = newname.replace("MCG+", "MCG +", 1)
    newname = newname.replace("MCG-", "MCG -", 1)
    newname = newname.replace("M+", "MCG +", 1)
    newname = newname.replace("M-", "MCG -", 1)
    newname = newname.replace("MGC ", "MCG ", 1)
    newname = newname.replace("Mrk", "MRK", 1)
    newname = newname.replace("MRK", "MRK ", 1)
    newname = newname.replace("NGC", "NGC ", 1)
    newname = newname.replace("PGC", "PGC ", 1)
    newname = newname.replace("SDSS", "SDSS ", 1)
    newname = newname.replace("UGC", "UGC ", 1)
    if newname.startswith('MESSIER '):
        newname = newname.replace('MESSIER ', 'M', 1)
    if newname.startswith('M ') and is_number(newname[2:]):
        newname = newname.replace('M ', 'M', 1)
    if newname.startswith('M') and is_number(newname[1:]):
        newname = 'M' + newname[1:].lstrip(" 0")
    if len(newname) > 4 and newname.startswith("PGC "):
        newname = newname[:4] + newname[4:].lstrip(" 0")
    if len(newname) > 4 and newname.startswith("UGC "):
        newname = newname[:4] + newname[4:].lstrip(" 0")
    if len(newname) > 5 and newname.startswith(("MCG +", "MCG -")):
        newname = newname[:5] + '-'.join([x.zfill(2)
                                          for x in
                                          newname[5:].strip().split("-")])
    if len(newname) > 5 and newname.startswith("CGCG "):
        newname = newname[:5] + '-'.join([x.zfill(3)
                                          for x in
                                          newname[5:].strip().split("-")])
    if ((len(newname) > 1 and newname.startswith("E")) or
            (len(newname) > 3 and newname.startswith('ESO'))):
        if newname[0] == "E":
            esplit = newname[1:].split("-")
        else:
            esplit = newname[3:].split("-")
        if len(esplit) == 2 and is_number(esplit[0].strip()):
            if esplit[1].strip()[0] == 'G':
                parttwo = esplit[1][1:].strip()
            else:
                parttwo = esplit[1].strip()
            if is_number(parttwo.strip()):
                newname = 'ESO ' + \
                    esplit[0].lstrip('0') + '-G' + parttwo.lstrip('0')
    newname = ' '.join(newname.split())
    return newname


def same_tag_num(photo, val, tag, canbelist=False):
    issame = (
        (tag not in photo and not val) or
        (tag in photo and not val) or
        (tag in photo and
         ((not canbelist and Decimal(photo[tag]) == Decimal(val)) or
          (canbelist and
           ((isinstance(photo[tag], str) and isinstance(val, str) and
             Decimal(photo[tag]) == Decimal(val)) or
            (isinstance(photo[tag], list) and isinstance(val, list) and
             photo[tag] == val))))))
    return issame


def same_tag_str(photo, val, tag):
    issame = ((tag not in photo and not val) or (
        tag in photo and not val) or (tag in photo and photo[tag] == val))
    return issame


def clean_snname(string):
    newstring = string.replace(' ', '').upper()
    if (newstring[:2] == "SN"):
        head = newstring[:6]
        tail = newstring[6:]
        if len(tail) >= 2 and tail[1] != '?':
            tail = tail.lower()
        newstring = head + tail

    return newstring


def get_event_text(eventfile):
    import gzip
    if eventfile.split('.')[-1] == 'gz':
        with gzip.open(eventfile, 'rt') as f:
            filetext = f.read()
    else:
        with open(eventfile, 'r') as f:
            filetext = f.read()
    return filetext
