#!/usr/local/bin/python3.5

import csv
import glob
import os
import re
import urllib.request
import calendar
import sys
import json
import codecs
import resource
import argparse
from cdecimal import Decimal
from astroquery.vizier import Vizier

from astropy.time import Time as astrotime
from astropy.cosmology import Planck15 as cosmo
from collections import OrderedDict
from math import log10, floor, sqrt, isnan
from bs4 import BeautifulSoup, Tag, NavigableString
from string import ascii_letters

parser = argparse.ArgumentParser(description='Generate a catalog JSON file and plot HTML files from SNE data.')
parser.add_argument('--update', '-u', dest='update', help='Only update catalog using live sources.',    default=False, action='store_true')
args = parser.parse_args()

clight = 29979245800.

eventnames = []

tasks = {
    "internal":       {"update": False},
    "vizier":         {"update": False},
    "suspect":        {"update": False},
    "cfa":            {"update": False},
    "ucb":            {"update": False},
    "sdss":           {"update": False},
    "gaia":           {"update": False},
    "csp":            {"update": False},
    "itep":           {"update": False},
    "asiago":         {"update": False},
    "rochester":      {"update": True },
    "lennarz":        {"update": False},
    "ogle":           {"update": True },
    "snls":           {"update": False},
    "nedd":           {"update": False},
    "cfaiaspectra":   {"update": False},
    "cfaibcspectra":  {"update": False},
    "snlsspectra":    {"update": False},
    "cspspectra":     {"update": False},
    "ucbspectra":     {"update": False},
    "suspectspectra": {"update": False},
    "snfspectra":     {"update": False},
    "writeevents":    {"update": True }
}

events = OrderedDict()

with open('rep-folders.txt', 'r') as f:
    repfolders = f.read().splitlines()

repyears = [int(repfolders[x][-4:]) for x in range(len(repfolders))]
repyears[0] -= 1

typereps = {
    'I P':    ['I pec', 'I-pec', 'I Pec', 'I-Pec'],
    'Ia P':   ['Ia pec', 'Ia-pec', 'Iapec'],
    'Ib P':   ['Ib pec', 'Ib-pec'],
    'Ic P':   ['Ic pec', 'Ic-pec'],
    'Ib/c':   ['Ibc'],
    'Ib/c P': ['Ib/c-pec'],
    'II P':   ['II pec', 'IIpec', 'II Pec', 'IIPec', 'IIP', 'IIp', 'II p', 'II-pec', 'II P pec', 'II-P'],
    'II L':   ['IIL'],
    'IIn P':  ['IIn pec', 'IIn-pec'],
    'IIb P':  ['IIb-pec', 'IIb: pec']
}

repbetterquanta = {
    'redshift',
    'ebv',
    'hvel',
    'lumdist'
}

def event_attr_priority(attr):
    if attr == 'photometry' or attr == 'spectra':
        return 'zzzzzzzz'
    if attr == 'name':
        return 'aaaaaaaa'
    if attr == 'aliases':
        return 'aaaaaaab'
    if attr == 'sources':
        return 'aaaaaaac'
    return attr

def add_event(name, load = True, delete = True):
    if name not in events or 'stub' in events[name]:
        if load:
            newname = load_event_from_file(name = name, delete = delete)
            if newname:
                if 'stub' in events[newname]:
                    raise(ValueError('Failed to find event file for stubbed event'))
                print('Loaded event ' + newname)
                return newname

        matches = []
        for event in events:
            if len(events[event]['aliases']) > 1 and name in events[event]['aliases']:
                matches.append(event)
        if len(matches) == 1:
            return matches[0]
        elif len(matches) > 1:
            raise(ValueError('Error, multiple matches to event, need event merging'))
        events[name] = OrderedDict()
        events[name]['name'] = name
        add_alias(name, name)
        print('Added new event ' + name)
        return name
    else:
        return name

def event_filename(name):
    return(name.replace('/', '_'))

def add_alias(name, alias):
    if 'aliases' in events[name]:
        if alias not in events[name]['aliases']:
            events[name].setdefault('aliases',[]).append(alias)
    else:
        events[name]['aliases'] = [alias]

def snname(string):
    newstring = string.replace(' ', '').upper()
    if (newstring[:2] == "SN"):
        head = newstring[:6]
        tail = newstring[6:]
        if len(tail) >= 2 and tail[1] != '?':
            tail = tail.lower()
        newstring = head + tail

    return newstring

def get_sig_digits(x):
    return len((''.join(x.split('.'))).strip('0'))

def round_sig(x, sig=4):
    if x == 0.0:
        return 0.0
    return round(x, sig-int(floor(log10(abs(x))))-1)

def pretty_num(x, sig=4):
    return str('%g'%(round_sig(x, sig)))

def get_source(name, reference = '', url = '', bibcode = '', secondary = ''):
    nsources = len(events[name]['sources']) if 'sources' in events[name] else 0
    if not reference:
        if not bibcode:
            raise(ValueError('Bibcode must be specified if name is not.'))
        else:
            if url:
                print('Warning: Reference URL ignored if bibcode specified')
        reference = bibcode
        url = "http://adsabs.harvard.edu/abs/" + bibcode
    if 'sources' not in events[name] or reference not in [events[name]['sources'][x]['name'] for x in range(nsources)]:
        source = str(nsources + 1)
        newsource = OrderedDict()
        newsource['name'] = reference
        if url:
            newsource['url'] = url
        if bibcode:
            newsource['bibcode'] = bibcode
        newsource['alias'] =  source
        if secondary:
            newsource['secondary'] = True
        events[name].setdefault('sources',[]).append(newsource)
    else:
        sourcexs = range(nsources)
        source = [events[name]['sources'][x]['alias'] for x in sourcexs][
            [events[name]['sources'][x]['name'] for x in sourcexs].index(reference)]
    return source

def add_photometry(name, timeunit = "MJD", time = "", instrument = "", band = "", magnitude = "", e_magnitude = "", source = "", upperlimit = False, system = ""):
    if not time or not magnitude:
        print('Warning: Time or AB mag not specified when adding photometry.\n')
        print('Name : "' + name + '", Time: "' + time + '", Band: "' + band + '", AB mag: "' + magnitude + '"')
        return

    if not is_number(time) or not is_number(magnitude):
        print('Warning: Time or AB mag not numerical.\n')
        print('Name : "' + name + '", Time: "' + time + '", Band: "' + band + '", AB mag: "' + magnitude + '"')
        return

    if e_magnitude and not is_number(e_magnitude):
        print('Warning: AB error not numerical.\n')
        print('Name : "' + name + '", Time: "' + time + '", Band: "' + band + '", AB err: "' + e_magnitude + '"')
        return

    # Look for duplicate data and don't add if duplicate
    if 'photometry' in events[name]:
        for photo in events[name]['photometry']:
            if (photo['timeunit'] == timeunit and
                Decimal(photo['time']) == Decimal(time) and
                Decimal(photo['magnitude']) == Decimal(magnitude) and
                (('band' not in photo and not band) or
                 ('band' in photo and photo['band'] == band) or
                 ('band' in photo and not band)) and
                (('error' not in photo and not e_magnitude) or
                 ('error' in photo and e_magnitude and Decimal(photo['error']) == Decimal(e_magnitude)) or
                 ('error' in photo and not e_magnitude)) and
                (('system' not in photo and not system) or
                 ('system' in photo and photo['system'] == system) or
                 ('system' in photo and not system))):
                return

    photoentry = OrderedDict()
    photoentry['timeunit'] = timeunit
    photoentry['time'] = str(time)
    if band:
        photoentry['band'] = band
    if system:
        photoentry['system'] = system
    photoentry['magnitude'] = str(magnitude)
    if instrument:
        photoentry['instrument'] = instrument
    if e_magnitude:
        photoentry['error'] = str(e_magnitude)
    if source:
        photoentry['source'] = source
    if upperlimit:
        photoentry['upperlimit'] = upperlimit
    events[name].setdefault('photometry',[]).append(photoentry)

def add_spectrum(name, waveunit, fluxunit, wavelengths, fluxes, timeunit = "", time = "", instrument = "",
    deredshifted = "", dereddened = "", errorunit = "", errors = "", source = "", snr = "",
    observer = "", reducer = ""):
    if not waveunit:
        'Warning: No error unit specified, not adding spectrum.'
        return
    if not fluxunit:
        'Warning: No flux unit specified, not adding spectrum.'
        return
    spectrumentry = OrderedDict()
    if deredshifted != '':
        spectrumentry['deredshifted'] = deredshifted
    if dereddened != '':
        spectrumentry['dereddened'] = dereddened
    if instrument:
        spectrumentry['instrument'] = instrument
    if timeunit:
        spectrumentry['timeunit'] = timeunit
    if time:
        spectrumentry['time'] = time
    if snr:
        spectrumentry['snr'] = snr
    if observer:
        spectrumentry['observer'] = observer
    if reducer:
        spectrumentry['reducer'] = reducer

    spectrumentry['waveunit'] = waveunit
    spectrumentry['fluxunit'] = fluxunit
    if errors and max([float(x) for x in errors]) > 0.:
        if not errorunit:
            'Warning: No error unit specified, not adding spectrum.'
            return
        spectrumentry['errorunit'] = errorunit
        data = [wavelengths, fluxes, errors]
    else:
        data = [wavelengths, fluxes]
    spectrumentry['data'] = [list(i) for i in zip(*data)]
    if source:
        spectrumentry['source'] = source
    events[name].setdefault('spectra',[]).append(spectrumentry)

def add_quanta(name, quanta, value, sources, forcereplacebetter = False, error = '', unit = '', kind = ''):
    if not quanta:
        raise(ValueError('Quanta must be specified for add_quanta.'))
    svalue = value.strip()
    serror = error.strip()
    if not svalue or svalue == '--' or svalue == '-':
        return
    if serror and (not is_number(serror) or float(serror) < 0.):
        raise(ValueError('Quanta error value must be a number and positive.'))

    #Handle certain quanta
    if quanta == 'hvel' or quanta == 'redshift':
        if not is_number(value):
            return
    if quanta == 'host':
        svalue = svalue.strip("()")
        svalue = svalue.replace("NGC", "NGC ")
        svalue = svalue.replace("UGC", "UGC ")
        svalue = svalue.replace("IC", "IC ")
        svalue = svalue.replace("Mrk", "MRK")
        svalue = svalue.replace("MRK", "MRK ")
        svalue = svalue.replace("PGC", "PGC ")
        svalue = ' '.join(svalue.split())
    elif quanta == 'claimedtype':
        for rep in typereps:
            if svalue in typereps[rep]:
                svalue = rep
                break
    elif quanta == 'snra' or quanta == 'sndec' or quanta == 'galra' or quanta == 'galdec':
        if unit == 'decdeg' or unit == 'radeg':
            deg = float('%g' % Decimal(svalue))
            sig = get_sig_digits(svalue)
            if unit == 'radeg':
                flhours = deg / 360.0 * 24.0
                hours = floor(flhours)
                minutes = floor((flhours - hours) * 60.0)
                seconds = (flhours * 60.0 - (hours * 60.0 + minutes)) * 60.0
                svalue = str(hours).zfill(2) + ':' + str(minutes).zfill(2) + ':' + pretty_num(seconds, sig = sig - 3).zfill(2)
            if unit == 'decdeg':
                fldeg = abs(deg)
                degree = floor(fldeg)
                minutes = floor((fldeg - degree) * 60.0)
                seconds = (fldeg * 60.0 - (degree * 60.0 + minutes)) * 60.0
                svalue = ('+' if deg >= 0.0 else '-') + str(degree).strip('+-').zfill(2) + ':' + str(minutes).zfill(2) + ':' + pretty_num(seconds, sig = sig - 3).zfill(2)
        elif unit == 'decdms':
            svalue = svalue.replace(' ', ':')
            valuesplit = svalue.split(':')
            svalue = ('+' if float(valuesplit[0]) > 0.0 else '-') + valuesplit[0].strip('+-').zfill(2) + ':' + ':'.join(valuesplit[1:]) if len(valuesplit) > 1 else ''
        elif unit == 'ranospace':
            svalue = svalue[:2] + ':' + svalue[2:4] + ((':' + svalue[4:]) if len(svalue) > 4 else '')
        elif unit == 'decnospace':
            svalue = svalue[:3] + ':' + svalue[3:5] + ((':' + svalue[5:]) if len(svalue) > 5 else '')
        else:
            svalue = svalue.replace(' ', ':')

    if is_number(svalue):
        svalue = '%g' % Decimal(svalue)
    if serror:
        serror = '%g' % Decimal(serror)

    if quanta in events[name]:
        for i, ct in enumerate(events[name][quanta]):
            if ct['value'] == svalue and sources:
                for source in sources.split(','):
                    if source not in events[name][quanta][i]['source'].split(','):
                        events[name][quanta][i]['source'] += ',' + source
                        if serror and 'error' not in events[name][quanta][i]:
                            events[name][quanta][i]['error'] = serror
                return

    quantaentry = OrderedDict()
    quantaentry['value'] = svalue
    if serror:
        quantaentry['error'] = serror
    if sources:
        quantaentry['source'] = sources
    if kind:
        quantaentry['kind'] = kind
    if (forcereplacebetter or quanta in repbetterquanta) and quanta in events[name]:
        newquantas = []
        isworse = True
        newsig = get_sig_digits(svalue)
        for ct in events[name][quanta]:
            if 'error' in ct:
                if serror:
                    if float(serror) < float(ct['error']):
                        isworse = False
                        continue
                newquantas.append(ct)
            else:
                if serror:
                    isworse = False
                    continue
                oldsig = get_sig_digits(ct['value'])
                if oldsig >= newsig:
                    newquantas.append(ct)
                if newsig >= oldsig:
                    isworse = False
        if not isworse:
            newquantas.append(quantaentry)
        events[name][quanta] = newquantas
    else:
        events[name].setdefault(quanta,[]).append(quantaentry)

def get_max_light(name):
    if 'photometry' not in events[name]:
        return (None, None)

    eventphoto = [Decimal(events[name]['photometry'][x]['magnitude']) for x in range(len(events[name]['photometry']))]
    mlmag = min(eventphoto)
    mlindex = eventphoto.index(mlmag)
    mlmjd = float(events[name]['photometry'][mlindex]['time'])
    return (astrotime(mlmjd, format='mjd').datetime, mlmag)

def get_first_light(name):
    if 'photometry' not in events[name]:
        return None

    eventtime = [Decimal(events[name]['photometry'][x]['time']) for x in range(len(events[name]['photometry']))]
    flindex = eventtime.index(min(eventtime))
    flmjd = float(events[name]['photometry'][flindex]['time'])
    return astrotime(flmjd, format='mjd').datetime

def set_first_max_light(name):
    if 'maxappmag' not in events[name]:
        (mldt, mlmag) = get_max_light(name)
        if mldt:
            add_quanta(name, 'maxyear', pretty_num(mldt.year), 'D')
            add_quanta(name, 'maxmonth', pretty_num(mldt.month), 'D')
            add_quanta(name, 'maxday', pretty_num(mldt.day), 'D')
            add_quanta(name, 'maxappmag', pretty_num(mlmag), 'D')

    if 'discovermonth' not in events[name] or 'discoverday' not in events[name]:
        fldt = get_first_light(name)
        if fldt:
            add_quanta(name, 'discoveryear', pretty_num(fldt.year), 'D')
            add_quanta(name, 'discovermonth', pretty_num(fldt.month), 'D')
            add_quanta(name, 'discoverday', pretty_num(fldt.day), 'D')

def jd_to_mjd(jd):
    return jd - Decimal(2400000.5)

def utf8(x):
    return str(x, 'utf-8')

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def convert_aq_output(row):
    return OrderedDict([(x, str(row[x]) if is_number(row[x]) else row[x]) for x in row.colnames])

def derive_and_sanitize():
    # Calculate some columns based on imported data, sanitize some fields
    for name in events:
        set_first_max_light(name)
        if 'claimedtype' in events[name]:
            events[name]['claimedtype'][:] = [ct for ct in events[name]['claimedtype'] if (ct['value'] != '?' and ct['value'] != '-')]
        if 'redshift' in events[name] and 'hvel' not in events[name]:
            # Find the "best" redshift to use for this
            bestsig = 0
            for z in events[name]['redshift']:
                sig = get_sig_digits(z['value'])
                if sig > bestsig:
                    bestz = z['value']
                    bestsig = sig
            if bestsig > 0:
                bestz = float(bestz)
                add_quanta(name, 'hvel', pretty_num(clight/1.e5*((bestz + 1.)**2. - 1.)/
                    ((bestz + 1.)**2. + 1.), sig = bestsig), 'D')
        elif 'hvel' in events[name] and 'redshift' not in events[name]:
            # Find the "best" hvel to use for this
            bestsig = 0
            for hv in events[name]['hvel']:
                sig = get_sig_digits(hv['value'])
                if sig > bestsig:
                    besthv = hv['value']
                    bestsig = sig
            if bestsig > 0 and is_number(besthv):
                voc = float(besthv)*1.e5/clight
                add_quanta(name, 'redshift', pretty_num(sqrt((1. + voc)/(1. - voc)) - 1., sig = bestsig), 'D')
        if 'maxabsmag' not in events[name] and 'maxappmag' in events[name] and 'lumdist' in events[name]:
            # Find the "best" distance to use for this
            bestsig = 0
            for ld in events[name]['lumdist']:
                sig = get_sig_digits(ld['value'])
                if sig > bestsig:
                    bestld = ld['value']
                    bestsig = sig
            if bestsig > 0 and is_number(bestld) and float(bestld) > 0.:
                add_quanta(name, 'maxabsmag', pretty_num(float(events[name]['maxappmag'][0]['value']) -
                    5.0*(log10(float(bestld)*1.0e6) - 1.0), sig = bestsig), 'D')
        if 'redshift' in events[name]:
            # Find the "best" redshift to use for this
            bestsig = 0
            for z in events[name]['redshift']:
                sig = get_sig_digits(z['value'])
                if sig > bestsig:
                    bestz = z['value']
                    bestsig = sig
            if bestsig > 0 and float(bestz) > 0.:
                if 'lumdist' not in events[name]:
                    dl = cosmo.luminosity_distance(float(bestz))
                    add_quanta(name, 'lumdist', pretty_num(dl.value, sig = bestsig), 'D')
                    if 'maxabsmag' not in events[name] and 'maxappmag' in events[name]:
                        add_quanta(name, 'maxabsmag', pretty_num(float(events[name]['maxappmag'][0]['value']) -
                            5.0*(log10(dl.to('pc').value) - 1.0), sig = bestsig), 'D')
        if 'photometry' in events[name]:
            events[name]['photometry'].sort(key=lambda x: (float(x['time']),
                x['band'] if 'band' in x else '', float(x['magnitude'])))
        if 'spectra' in events[name] and list(filter(None, ['time' in x for x in events[name]['spectra']])):
            events[name]['spectra'].sort(key=lambda x: float(x['time']))
        events[name] = OrderedDict(sorted(events[name].items(), key=lambda key: event_attr_priority(key[0])))

def delete_old_event_files():
    # Delete all old event JSON files
    for folder in repfolders:
        filelist = glob.glob("../" + folder + "/*.json") + glob.glob("../" + folder + "/*.json.gz")
        for f in filelist:
            os.remove(f)

def write_all_events(empty = False):
    # Write it all out!
    for name in events:
        if 'stub' in events[name]:
            if not empty:
                continue
            else:
                del(events[name]['stub'])
        print('Writing ' + name)
        filename = event_filename(name)

        jsonstring = json.dumps({name:events[name]}, indent='\t', separators=(',', ':'), ensure_ascii=False)

        outdir = '../'
        if 'discoveryear' in events[name]:
            for r, year in enumerate(repyears):
                if int(events[name]['discoveryear'][0]['value']) <= year:
                    outdir += repfolders[r]
                    break
        else:
            outdir += str(repfolders[0])

        f = codecs.open(outdir + '/' + filename + '.json', 'w', encoding='utf8')
        f.write(jsonstring)
        f.close()

def load_event_from_file(name = '', location = '', clean = False, delete = True):
    if not name and not location:
        raise ValueError('Either event name or location must be specified to load event')

    if location:
        path = location
    else:
        indir = '../'
        path = ''
        for rep in repfolders:
            newpath = indir + rep + '/' + name + '.json'
            if os.path.isfile(newpath):
                path = newpath

    if not path:
        return False
    else:
        with open(path, 'r') as f:
            if name in events:
                del events[name]
            events.update(json.loads(f.read(), object_pairs_hook=OrderedDict))
            name = next(reversed(events))
            if clean:
                clean_event(name)
        if 'writeevents' in tasks and delete:
            os.remove(path)
        return name

def clean_event(name):
    bibcodes = []
    if 'name' not in events[name]:
        events[name]['name'] = name
    if 'aliases' not in events[name]:
        add_alias(name, name)
    if 'sources' in events[name]:
        for s, source in enumerate(events[name]['sources']):
            if 'bibcode' in source and 'name' not in source:
                bibcodes.append(source['bibcode'])
        if bibcodes:
            del events[name]['sources']
            for bibcode in bibcodes:
                get_source(name, bibcode = bibcode)
    if 'photometry' in events[name]:
        for p, photo in enumerate(events[name]['photometry']):
            if photo['timeunit'] == 'JD':
                events[name]['photometry'][p]['timeunit'] = 'MJD'
                events[name]['photometry'][p]['time'] = str(jd_to_mjd(Decimal(photo['time'])))
            if bibcodes and 'source' not in photo:
                alias = get_source(name, bibcode = bibcodes[0])
                events[name]['photometry'][p]['source'] = alias

def do_task(task):
    return task in tasks and (not args.update or tasks[task]['update'])

def journal_events(clear = True):
    if 'writeevents' in tasks:
        write_all_events()
    if clear:
        clear_events()

def clear_events():
    global events
    events = OrderedDict((k, OrderedDict((('name', events[k]['name']), ('aliases', events[k]['aliases']), ('stub', True)))) for k in events)

# Either load stubs of each event (if updating) or delete all event files (if starting fresh)
if 'writeevents' in tasks:
    if args.update:
        files = []
        for rep in repfolders:
            files += glob.glob('../' + rep + "/*.json")

        for fi in files:
            name = os.path.basename(os.path.splitext(fi)[0])
            name = add_event(name, delete = False)
            events[name] = OrderedDict((('name', events[name]['name']), ('aliases', events[name]['aliases']), ('stub', True)))
    else:
        delete_old_event_files()

# Import data provided directly to OSC
if do_task('internal'):
    for datafile in sorted(glob.glob("../sne-internal/*.json"), key=lambda s: s.lower()):
        if not load_event_from_file(location = datafile, clean = True, delete = False):
            raise IOError('Failed to find specified file.')
    journal_events()

# Import primary data sources from Vizier
if do_task('vizier'):
    Vizier.ROW_LIMIT = -1

    # 2014MNRAS.444.3258M
    result = Vizier.get_catalogs("J/MNRAS/444/3258/SNe")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    oldname = ''
    for row in table:
        name = row['SN']
        if name == oldname:
            continue
        oldname = name
        name = add_event(name)
        source = get_source(name, bibcode = '2014MNRAS.444.3258M')
        add_quanta(name, 'redshift', str(row['z']), source, kind = 'heliocentric', error = str(row['e_z']))
        add_quanta(name, 'snra', str(row['_RA']), source, unit = 'radeg')
        add_quanta(name, 'sndec', str(row['_DE']), source, unit = 'decdeg')

    # 2014MNRAS.438.1391P
    result = Vizier.get_catalogs("J/MNRAS/438/1391/table2")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in table:
        name = row['SN']
        name = add_event(name)
        source = get_source(name, bibcode = '2014MNRAS.438.1391P')
        add_quanta(name, 'redshift', str(row['zh']), source, kind = 'heliocentric')
        add_quanta(name, 'snra', row['RAJ2000'], source)
        add_quanta(name, 'sndec', row['DEJ2000'], source)

    # 2012ApJ...749...18B
    result = Vizier.get_catalogs("J/ApJ/749/18/table1")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in table:
        name = row['Name'].replace(' ','')
        name = add_event(name)
        source = get_source(name, bibcode = '2012ApJ...749...18B')
        mjd = str(astrotime(2450000.+row['JD'], format='jd').mjd)
        band = row['Filt']
        magnitude = str(row['mag'])
        e_magnitude = str(row['e_mag'])
        e_magnitude = '' if e_magnitude == '--' else e_magnitude
        upperlimit = True if row['l_mag'] == '>' else False
        add_photometry(name, time = mjd, band = band, magnitude = magnitude, e_magnitude = e_magnitude, instrument = 'Swift/UVOT',
            source = source, upperlimit = upperlimit)

    # 2010A&A...523A...7G
    result = Vizier.get_catalogs("J/A+A/523/A7/table9")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in table:
        name = 'SNLS-' + row['SNLS']
        name = add_event(name)
        source = get_source(name, bibcode = '2010A&A...523A...7G')
        astrot = astrotime(2450000.+row['Date1'], format='jd')
        add_quanta(name, 'discoverday', str(astrot.datetime.day), source)
        add_quanta(name, 'discovermonth', str(astrot.datetime.month), source)
        add_quanta(name, 'discoveryear', str(astrot.datetime.year), source)
        add_quanta(name, 'ebv', str(row['E_B-V_']), source)
        add_quanta(name, 'redshift', str(row['z']), source)
        add_quanta(name, 'claimedtype', row['Type'].replace('*', '?').replace('SN','').replace('(pec)',' P'), source)
        add_quanta(name, 'snra', row['RAJ2000'], source)
        add_quanta(name, 'sndec', row['DEJ2000'], source)

    # 2004A&A...415..863G
    result = Vizier.get_catalogs("J/A+A/415/863/table1")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in table:
        name = 'SN' + row['SN']
        name = add_event(name)
        source = get_source(name, bibcode = '2004A&A...415..863G')
        datesplit = row['Date'].split('-')
        add_quanta(name, 'discoverday', datesplit[2].lstrip('0'), source)
        add_quanta(name, 'discovermonth', datesplit[1].lstrip('0'), source)
        add_quanta(name, 'discoveryear', datesplit[0], source)
        add_quanta(name, 'host', 'Abell ' + str(row['Abell']), source)
        add_quanta(name, 'claimedtype', row['Type'], source)
        add_quanta(name, 'snra', row['RAJ2000'], source)
        add_quanta(name, 'sndec', row['DEJ2000'], source)

    # 2010ApJ...708..661D
    result = Vizier.get_catalogs("J/ApJ/708/661/sn")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in table:
        name = row['SN']
        if not name:
            name = 'SDSS-II ' + str(row['SDSS-II'])
        else:
            name = 'SN' + name
        name = add_event(name)
        source = get_source(name, bibcode = '2010ApJ...708..661D')
        add_alias(name, 'SDSS-II ' + str(row['SDSS-II']))
        add_quanta(name, 'snra', row['RAJ2000'], source)
        add_quanta(name, 'sndec', row['DEJ2000'], source)

    result = Vizier.get_catalogs("J/ApJ/708/661/table1")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in table:
        if row['f_SN'] == 'a':
            name = 'SDSS-II ' + str(row['SN'])
        else:
            name = 'SN' + row['SN']
        name = add_event(name)
        source = get_source(name, bibcode = '2010ApJ...708..661D')
        add_quanta(name, 'redshift', str(row['z']), source, error = str(row['e_z']))

    # 2014ApJ...795...44R
    result = Vizier.get_catalogs("J/ApJ/795/44/ps1_snIa")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in table:
        name = row['SN']
        name = add_event(name)
        source = get_source(name, bibcode = '2014ApJ...795...44R')
        astrot = astrotime(row['tdisc'], format='mjd')
        add_quanta(name, 'discoverday', str(astrot.datetime.day), source)
        add_quanta(name, 'discovermonth', str(astrot.datetime.month), source)
        add_quanta(name, 'discoveryear',  str(astrot.datetime.year), source)
        add_quanta(name, 'redshift', str(row['z']), source, error = str(row['e_z']))
        add_quanta(name, 'snra', row['RAJ2000'], source)
        add_quanta(name, 'sndec', row['DEJ2000'], source)

    result = Vizier.get_catalogs("J/ApJ/795/44/table6")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in table:
        name = row['SN']
        name = add_event(name)
        source = get_source(name, bibcode = '2014ApJ...795...44R')
        if row['mag'] != '--':
            add_photometry(name, time = str(row['MJD']), band = row['Filt'], magnitude = str(row['mag']),
                e_magnitude = str(row['e_mag']), source = source, system = 'AB')

    # 1990A&AS...82..145C
    result = Vizier.get_catalogs("II/189/mag")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)

    with open('../sne-external/II_189_refs.csv') as f:
        tsvin = csv.reader(f, delimiter='\t', skipinitialspace=True)
        ii189bibdict = {}
        ii189refdict = {}
        for r, row in enumerate(tsvin):
            if row[0] != '0':
                ii189bibdict[r+1] = row[1]
            else:
                ii189refdict[r+1] = row[2]

    for row in table:
        if row['band'][0] == '(':
            continue
        name = 'SN' + row['SN']
        name = add_event(name)
        source = ''
        secsource = get_source(name, bibcode = '1990A&AS...82..145C', secondary = True)
        mjd = str(jd_to_mjd(Decimal(row['JD'])))
        mag = str(row['m'])
        band = row['band'].strip("'")
        if row['r_m'] in ii189bibdict:
            source = get_source(name, bibcode = ii189bibdict[row['r_m']])
        else:
            source = get_source(name, reference = ii189refdict[row['r_m']])

        add_photometry(name, time = mjd, band = band, magnitude = mag, source = ','.join([source,secsource]))

    result = Vizier.get_catalogs("VII/272/snrs")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)

    # 2014yCat.7272....0G
    for row in table:
        row = convert_aq_output(row)
        name = ''
        if row["Names"]:
            names = row["Names"].split(',')
            for nam in names:
                if nam.strip()[:2] == 'SN':
                    name = nam.strip()
                    if is_number(name[2:]):
                        name = name + 'A'
            if not name:
                for nam in names:
                    if nam.strip('()') == nam:
                        name = nam.strip()
                        break
        if not name:
            name = row["SNR"].strip()

        name = add_event(name)
        source = get_source(name, bibcode = '2014yCat.7272....0G')

        if row["Names"]:
            names = row["Names"].split(',')
            for nam in names:
                add_alias(name, nam.strip('()'))
                if nam.strip()[:2] == 'SN':
                    add_quanta(name, 'discoveryear', nam.strip()[2:], source)

        add_quanta(name, 'host', 'Milky Way', source)
        add_quanta(name, 'snra', row['RAJ2000'], source)
        add_quanta(name, 'sndec', row['DEJ2000'], source, unit = 'decdms')

    # 2014MNRAS.442..844F
    result = Vizier.get_catalogs("J/MNRAS/442/844/table1")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in table:
        row = convert_aq_output(row)
        name = 'SN' + row['SN']
        name = add_event(name)
        source = get_source(name, bibcode = '2014MNRAS.442..844F')
        add_quanta(name, 'redshift', str(row['zhost']), source)
        add_quanta(name, 'ebv', str(row['E_B-V_']), source)

    result = Vizier.get_catalogs("J/MNRAS/442/844/table2")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in table:
        row = convert_aq_output(row)
        name = 'SN' + str(row['SN'])
        name = add_event(name)
        source = get_source(name, bibcode = "2014MNRAS.442..844F")
        if 'Bmag' in row and is_number(row['Bmag']) and not isnan(float(row['Bmag'])):
            add_photometry(name, time = row['MJD'], band = 'B', magnitude = row['Bmag'], e_magnitude = row['e_Bmag'], source = source)
        if 'Vmag' in row and is_number(row['Vmag']) and not isnan(float(row['Vmag'])):
            add_photometry(name, time = row['MJD'], band = 'V', magnitude = row['Vmag'], e_magnitude = row['e_Vmag'], source = source)
        if 'Rmag' in row and is_number(row['Rmag']) and not isnan(float(row['Rmag'])):
            add_photometry(name, time = row['MJD'], band = 'R', magnitude = row['Rmag'], e_magnitude = row['e_Rmag'], source = source)
        if 'Imag' in row and is_number(row['Imag']) and not isnan(float(row['Imag'])):
            add_photometry(name, time = row['MJD'], band = 'I', magnitude = row['Imag'], e_magnitude = row['e_Imag'], source = source)

    # 2012MNRAS.425.1789S
    result = Vizier.get_catalogs("J/MNRAS/425/1789/table1")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in table:
        row = convert_aq_output(row)
        name = ''.join(row['SimbadName'].split(' '))
        name = add_event(name)
        add_alias(name, 'SN' + row['SN'])
        source = get_source(name, bibcode = '2012MNRAS.425.1789S')
        add_quanta(name, 'host', row['Gal'], source)
        add_quanta(name, 'hvel', str(row['cz']), source)
        add_quanta(name, 'ebv', str(row['E_B-V_']), source)

    # 2015ApJS..219...13W
    result = Vizier.get_catalogs("J/ApJS/219/13/table3")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in table:
        row = convert_aq_output(row)
        name = u'LSQ' + str(row['LSQ'])
        name = add_event(name)
        source = get_source(name, bibcode = "2015ApJS..219...13W")
        add_quanta(name, 'snra', row['RAJ2000'], source)
        add_quanta(name, 'sndec', row['DEJ2000'], source)
        add_quanta(name, 'redshift', row['z'], source, error = row['e_z'])
        add_quanta(name, 'ebv', row['E_B-V_'], source)
    result = Vizier.get_catalogs("J/ApJS/219/13/table2")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in table:
        row = convert_aq_output(row)
        name = 'LSQ' + row['LSQ']
        source = get_source(name, bibcode = "2015ApJS..219...13W")
        add_photometry(name, time = str(jd_to_mjd(Decimal(row['JD']))), instrument = 'La Silla-QUEST', band = row['Filt'],
            magnitude = row['mag'], e_magnitude = row['e_mag'], system = "Swope", source = source)
    journal_events()

# Suspect catalog
if do_task('suspect'): 
    with open('../sne-external/suspectreferences.csv','r') as f:
        tsvin = csv.reader(f, delimiter='\t', skipinitialspace=True)
        suspectrefdict = {}
        for row in tsvin:
            suspectrefdict[row[0]] = row[1]

    for datafile in sorted(glob.glob("../sne-external/SUSPECT/*.html"), key=lambda s: s.lower()):
        basename = os.path.basename(datafile)
        basesplit = basename.split('-')
        name = basesplit[1]
        name = add_event(name)
        if name[:2] == 'SN' and is_number(name[2:]):
            name = name + 'A'
        band = basesplit[3].split('.')[0]
        ei = int(basesplit[2])
        bandlink = 'file://' + os.path.abspath(datafile)
        bandresp = urllib.request.urlopen(bandlink)
        bandsoup = BeautifulSoup(bandresp, "html5lib")
        bandtable = bandsoup.find('table')
        if ei == 1:
            names = bandsoup.body.findAll(text=re.compile("Name"))
            reference = ''
            for link in bandsoup.body.findAll('a'):
                if 'adsabs' in link['href']:
                    reference = str(link).replace('"', "'")

            bibcode = suspectrefdict[reference]
            source = get_source(name, bibcode = bibcode)

            year = re.findall(r'\d+', name)[0]
            add_quanta(name, 'discoveryear', year, source)
            add_quanta(name, 'host', names[1].split(':')[1].strip(), source)

            redshifts = bandsoup.body.findAll(text=re.compile("Redshift"))
            if redshifts:
                add_quanta(name, 'redshift', redshifts[0].split(':')[1].strip(), source)
            hvels = bandsoup.body.findAll(text=re.compile("Heliocentric Velocity"))
            if hvels:
                add_quanta(name, 'hvel', hvels[0].split(':')[1].strip().split(' ')[0], source)
            types = bandsoup.body.findAll(text=re.compile("Type"))

            add_quanta(name, 'claimedtype', types[0].split(':')[1].strip().split(' ')[0], source)

        secondaryreference = "SUSPECT"
        secondaryrefurl = "https://www.nhn.ou.edu/~suspect/"
        secondarysource = get_source(name, reference = secondaryreference, url = secondaryrefurl, secondary = True)

        for r, row in enumerate(bandtable.findAll('tr')):
            if r == 0:
                continue
            col = row.findAll('td')
            mjd = str(jd_to_mjd(Decimal(col[0].contents[0])))
            mag = col[3].contents[0]
            if mag.isspace():
                mag = ''
            else:
                mag = str(mag)
            e_magnitude = col[4].contents[0]
            if e_magnitude.isspace():
                e_magnitude = ''
            else:
                e_magnitude = str(e_magnitude)
            add_photometry(name, time = mjd, band = band, magnitude = mag, e_magnitude = e_magnitude, source = secondarysource + ',' + source)
    journal_events()

# CfA data
if do_task('cfa'): 
    for file in sorted(glob.glob("../sne-external/cfa-input/*.dat"), key=lambda s: s.lower()):
        f = open(file,'r')
        tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)
        csv_data = []
        for r, row in enumerate(tsvin):
            new = []
            for item in row:
                new.extend(item.split("\t"))
            csv_data.append(new)

        for r, row in enumerate(csv_data):
            for c, col in enumerate(row):
                csv_data[r][c] = col.strip()
            csv_data[r] = [_f for _f in csv_data[r] if _f]

        eventname = os.path.basename(os.path.splitext(file)[0])

        eventparts = eventname.split('_')

        name = snname(eventparts[0])
        name = add_event(name)
        secondaryname = 'CfA Supernova Archive'
        secondaryurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
        secondarysource = get_source(name, reference = secondaryname, url = secondaryurl, secondary = True)

        year = re.findall(r'\d+', name)[0]
        add_quanta(name, 'discoveryear', year, secondarysource)

        eventbands = list(eventparts[1])

        tu = 'MJD'
        jdoffset = Decimal(0.)
        for rc, row in enumerate(csv_data):
            if len(row) > 0 and row[0][0] == "#":
                if len(row[0]) > 2 and row[0][:3] == "#JD":
                    tu = 'JD'
                    rowparts = row[0].split('-')
                    jdoffset = Decimal(rowparts[1])
                elif len(row[0]) > 6 and row[0][:7] == "#Julian":
                    tu = 'JD'
                    jdoffset = Decimal(0.)
                elif len(row) > 1 and row[1].lower() == "photometry":
                    for ci, col in enumerate(row[2:]):
                        if col[0] == "(":
                            refstr = ' '.join(row[2+ci:])
                            refstr = refstr.replace('(','').replace(')','')
                            bibcode = refstr
                            source = get_source(name, bibcode = bibcode)

                elif len(row) > 1 and row[1] == "HJD":
                    tu = "HJD"

                continue
            elif len(row) > 0:
                mjd = row[0]
                for v, val in enumerate(row):
                    if v == 0:
                        if tu == 'JD':
                            mjd = str(jd_to_mjd(Decimal(val) + jdoffset))
                            tuout = 'MJD'
                        elif tu == 'HJD':
                            mjd = str(jd_to_mjd(Decimal(val)))
                            tuout = 'MJD'
                        else:
                            mjd = val
                            tuout = tu
                    elif v % 2 != 0:
                        if float(row[v]) < 90.0:
                            add_photometry(name, timeunit = tuout, time = mjd, band = eventbands[(v-1)//2], magnitude = row[v], e_magnitude = row[v+1], source = secondarysource + ',' + source)
        f.close()

    # Hicken 2012
    f = open("../sne-external/hicken-2012-standard.dat", 'r')
    tsvin = csv.reader(f, delimiter='|', skipinitialspace=True)
    for r, row in enumerate(tsvin):
        if r <= 47:
            continue

        if row[0][:2] != 'sn':
            name = 'SN' + row[0].strip()
        else:
            name = row[0].strip()

        name = add_event(name)

        source = get_source(name, bibcode = '2012ApJS..200...12H')
        add_photometry(name, timeunit = 'MJD', time = row[2].strip(), band = row[1].strip(),
            magnitude = row[6].strip(), e_magnitude = row[7].strip(), source = source)
    
    # Bianco 2014
    tsvin = open("../sne-external/bianco-2014-standard.dat", 'r')
    tsvin = csv.reader(tsvin, delimiter=' ', skipinitialspace=True)
    for row in tsvin:
        name = 'SN' + row[0]
        name = add_event(name)

        source = get_source(name, bibcode = '2014ApJS..213...19B')
        add_photometry(name, timeunit = 'MJD', time = row[2], band = row[1], magnitude = row[3],
            e_magnitude = row[4], instrument = row[5], system = "Standard", source = source)
    f.close()
    journal_events()

# Now import the UCB SNDB
if do_task('ucb'): 
    for file in sorted(glob.glob("../sne-external/SNDB/*.dat"), key=lambda s: s.lower()):
        f = open(file,'r')
        tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)

        eventname = os.path.basename(os.path.splitext(file)[0])

        eventparts = eventname.split('.')

        name = snname(eventparts[0])
        name = add_event(name)

        reference = "UCB Filippenko Group's Supernova Database (SNDB)"
        refurl = "http://heracles.astro.berkeley.edu/sndb/info"
        source = get_source(name, reference = reference, url = refurl, secondary = True)

        year = re.findall(r'\d+', name)[0]
        add_quanta(name, 'discoveryear', year, source)

        for r, row in enumerate(tsvin):
            if len(row) > 0 and row[0] == "#":
                continue
            mjd = row[0]
            magnitude = row[1]
            e_magnitude = row[2]
            band = row[4]
            instrument = row[5]
            add_photometry(name, time = mjd, instrument = instrument, band = band, magnitude = magnitude,
                e_magnitude = e_magnitude, source = source)
        f.close()
    journal_events()
    
# Import SDSS
if do_task('sdss'): 
    sdssbands = ['u', 'g', 'r', 'i', 'z']
    for file in sorted(glob.glob("../sne-external/SDSS/*.sum"), key=lambda s: s.lower()):
        f = open(file,'r')
        tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)

        for r, row in enumerate(tsvin):
            if r == 0:
                if row[5] == "RA:":
                    name = "SDSS-II " + row[3]
                else:
                    name = "SN" + row[5]
                name = add_event(name)
                add_alias(name, "SDSS-II " + row[3])

                bibcode = '2008AJ....136.2306H'
                source = get_source(name, bibcode = bibcode)

                if row[5] != "RA:":
                    year = re.findall(r'\d+', name)[0]
                    add_quanta(name, 'discoveryear', year, source)

                add_quanta(name, 'snra', row[-4], source, unit = 'radeg')
                add_quanta(name, 'sndec', row[-2], source, unit = 'decdeg')
            if r == 1:
                add_quanta(name, 'redshift', row[2], source, error = row[4])
            if r >= 19:
                # Skip bad measurements
                if int(row[0]) > 1024:
                    continue

                mjd = row[1]
                band = sdssbands[int(row[2])]
                magnitude = row[3]
                e_magnitude = row[4]
                instrument = "SDSS"
                add_photometry(name, time = mjd, instrument = instrument, band = band, magnitude = magnitude,
                    e_magnitude = e_magnitude, source = source, system = "AB")
        f.close()
    journal_events()

#Import GAIA
if do_task('gaia'): 
    #response = urllib2.urlopen('https://gaia.ac.uk/selected-gaia-science-alerts')
    path = os.path.abspath('../sne-external/selected-gaia-science-alerts')
    response = urllib.request.urlopen('file://' + path)
    html = response.read()

    soup = BeautifulSoup(html, "html5lib")
    table = soup.findAll("table")[1]
    for r, row in enumerate(table.findAll('tr')):
        if r == 0:
            continue

        col = row.findAll('td')
        classname = col[7].contents[0]

        if 'SN' not in classname:
            continue

        links = row.findAll('a')
        name = links[0].contents[0]

        if name == 'Gaia15aaaa':
            continue

        name = add_event(name)

        reference = "Gaia Photometric Science Alerts"
        refurl = "https://gaia.ac.uk/selected-gaia-science-alerts"
        source = get_source(name, reference = reference, url = refurl)

        year = '20' + re.findall(r'\d+', name)[0]
        add_quanta(name, 'discoveryear', year, source)

        add_quanta(name, 'snra', col[2].contents[0].strip(), source, unit = 'radeg')
        add_quanta(name, 'sndec', col[3].contents[0].strip(), source, unit = 'decdeg')
        add_quanta(name, 'claimedtype', classname.replace('SN', '').strip(), source)

        photfile = '../sne-external/GAIA/GAIA-' + name + '.html'
        with open(photfile, 'r') as f:
            phottxt = f.read()

        photsoup = BeautifulSoup(phottxt, "html5lib")
        photodata = str(photsoup.contents[0]).split('\n')[2:-1]
        for ph in photodata:
            photo = ph.split(',')
            mjd = str(jd_to_mjd(Decimal(photo[1].strip())))
            magnitude = photo[2].strip()
            e_magnitude = 0.
            instrument = 'GAIA'
            band = 'G'
            add_photometry(name, time = mjd, instrument = instrument, band = band, magnitude = magnitude, e_magnitude = e_magnitude, source = source)
    journal_events()

# Import CSP
if do_task('csp'): 
    cspbands = ['u', 'B', 'V', 'g', 'r', 'i', 'Y', 'J', 'H', 'K']
    for file in sorted(glob.glob("../sne-external/CSP/*.dat"), key=lambda s: s.lower()):
        f = open(file,'r')
        tsvin = csv.reader(f, delimiter='\t', skipinitialspace=True)

        eventname = os.path.basename(os.path.splitext(file)[0])

        eventparts = eventname.split('opt+')

        name = snname(eventparts[0])
        name = add_event(name)

        reference = "Carnegie Supernova Project"
        refurl = "http://csp.obs.carnegiescience.edu/data"
        source = get_source(name, reference = reference, url = refurl)

        year = re.findall(r'\d+', name)[0]
        add_quanta(name, 'discoveryear', year, source)

        for r, row in enumerate(tsvin):
            if len(row) > 0 and row[0][0] == "#":
                if r == 2:
                    add_quanta(name, 'redshift', row[0].split(' ')[-1], source)
                    add_quanta(name, 'snra', row[1].split(' ')[-1], source)
                    add_quanta(name, 'sndec', row[2].split(' ')[-1], source)
                continue
            for v, val in enumerate(row):
                if v == 0:
                    mjd = val
                elif v % 2 != 0:
                    if float(row[v]) < 90.0:
                        add_photometry(name, time = mjd, instrument = 'CSP', band = cspbands[(v-1)//2], magnitude = row[v], e_magnitude = row[v+1], source = source)
        f.close()
    journal_events()

# Import ITEP
if do_task('itep'): 
    needsbib = []
    with open("../sne-external/itep-refs.txt",'r') as f:
        refrep = f.read().splitlines()
    refrepf = dict(list(zip(refrep[1::2], refrep[::2])))
    f = open("../sne-external/itep-lc-cat-28dec2015.txt",'r')
    tsvin = csv.reader(f, delimiter='|', skipinitialspace=True)
    curname = ''
    for r, row in enumerate(tsvin):
        if r <= 1 or len(row) < 7:
            continue
        name = 'SN' + row[0].strip()
        mjd = str(jd_to_mjd(Decimal(row[1].strip())))
        band = row[2].strip()
        magnitude = row[3].strip()
        e_magnitude = row[4].strip()
        reference = row[6].strip().strip(',')

        if curname != name:
            curname = name
            name = add_event(name)

            secondaryreference = "Sternberg Astronomical Institute Supernova Light Curve Catalogue"
            secondaryrefurl = "http://dau.itep.ru/sn/node/72"
            secondarysource = get_source(name, reference = secondaryreference, url = secondaryrefurl, secondary = True)

            year = re.findall(r'\d+', name)[0]
            add_quanta(name, 'discoveryear', year, secondarysource)
        if reference in refrepf:
            bibcode = refrepf[reference]
            source = get_source(name, bibcode = bibcode)
        else:
            needsbib.append(reference)
            source = get_source(name, reference = reference) if reference else ''

        add_photometry(name, time = mjd, band = band, magnitude = magnitude, e_magnitude = e_magnitude, source = secondarysource + ',' + source)
    f.close()
    
    # Write out references that could use a bibcode
    needsbib = list(OrderedDict.fromkeys(needsbib))
    with open('../itep-needsbib.txt', 'w') as f:
        f.writelines(["%s\n" % i for i in needsbib])
    journal_events()

# Now import the Asiago catalog
if do_task('asiago'): 
    #response = urllib.request.urlopen('http://graspa.oapd.inaf.it/cgi-bin/sncat.php')
    path = os.path.abspath('../sne-external/asiago-cat.php')
    response = urllib.request.urlopen('file://' + path)
    html = response.read().decode('utf-8')
    html = html.replace("\r", "")

    soup = BeautifulSoup(html, "html5lib")
    table = soup.find("table")

    records = []
    for r, row in enumerate(table.findAll('tr')):
        if r == 0:
            continue

        col = row.findAll('td')
        records.append([utf8(x.renderContents()) for x in col])

    for record in records:
        if len(record) > 1 and record[1] != '':
            name = snname("SN" + record[1]).strip('?')
            name = add_event(name)

            reference = 'Asiago Supernova Catalogue'
            refurl = 'http://graspa.oapd.inaf.it/cgi-bin/sncat.php'
            source = get_source(name, reference = reference, url = refurl, secondary = True)

            year = re.findall(r'\d+', name)[0]
            add_quanta(name, 'discoveryear', year, source)

            hostname = record[2]
            galra = record[3]
            galdec = record[4]
            snra = record[5]
            sndec = record[6]
            redvel = record[11].strip(':')
            discoverer = record[19]

            datestr = record[18]
            if "*" in datestr:
                monthkey = 'discovermonth'
                daykey = 'discoverday'
            else:
                monthkey = 'maxmonth'
                daykey = 'maxday'

            if datestr.strip() != '':
                dayarr = re.findall(r'\d+', datestr)
                if dayarr:
                    add_quanta(name, daykey, dayarr[0], source)
                monthstr = ''.join(re.findall("[a-zA-Z]+", datestr))
                add_quanta(name, monthkey, str(list(calendar.month_abbr).index(monthstr)), source)

            hvel = ''
            redshift = ''
            if redvel != '':
                if round(float(redvel)) == float(redvel):
                    hvel = int(redvel)
                else:
                    redshift = float(redvel)
                redshift = str(redshift)
                hvel = str(hvel)

            claimedtype = record[17].strip(':')

            if (hostname != ''):
                add_quanta(name, 'host', hostname, source)
            if (claimedtype != ''):
                add_quanta(name, 'claimedtype', claimedtype, source)
            if (redshift != ''):
                add_quanta(name, 'redshift', redshift, source)
            if (hvel != ''):
                add_quanta(name, 'hvel', hvel, source)
            if (galra != ''):
                add_quanta(name, 'galra', galra, source, unit = 'ranospace')
            if (galdec != ''):
                add_quanta(name, 'galdec', galdec, source, unit = 'decnospace')
            if (snra != ''):
                add_quanta(name, 'snra', snra, source, unit = 'ranospace')
            if (sndec != ''):
                add_quanta(name, 'sndec', sndec, source, unit = 'decnospace')
            if (discoverer != ''):
                add_quanta(name, 'discoverer', discoverer, source)
    journal_events()

if do_task('rochester'): 
    rochesterpaths = ['file://'+os.path.abspath('../sne-external/snredshiftall.html'), 'http://www.rochesterastronomy.org/sn2016/snredshift.html']
    rochesterupdate = [False, True]

    for p, path in enumerate(rochesterpaths):
        if args.update and not rochesterupdate[p]:
            continue
        response = urllib.request.urlopen(path)
        html = response.read()

        soup = BeautifulSoup(html, "html5lib")
        rows = soup.findAll('tr')
        secondaryreference = "Latest Supernovae"
        secondaryrefurl = "http://www.rochesterastronomy.org/snimages/snredshiftall.html"
        for r, row in enumerate(rows):
            if r == 0:
                continue
            cols = row.findAll('td')
            if not len(cols):
                continue

            name = ''
            if cols[14].contents:
                aka = str(cols[14].contents[0]).strip()
                if is_number(aka.strip('?')):
                    aka = 'SN' + aka.strip('?') + 'A'
                    name = add_event(aka)
                elif len(aka) >= 4 and is_number(aka[:4]):
                    aka = 'SN' + aka
                    name = add_event(aka)

            sn = re.sub('<[^<]+?>', '', str(cols[0].contents[0])).strip()
            if is_number(sn.strip('?')):
                sn = 'SN' + sn.strip('?') + 'A'
            elif len(sn) >= 4 and is_number(sn[:4]):
                sn = 'SN' + sn
            if not name:
                name = add_event(sn)

            if cols[14].contents:
                if aka == 'SNR G1.9+0.3':
                    aka = 'G001.9+00.3'
                add_alias(name, aka)

            reference = cols[12].findAll('a')[0].contents[0].strip()
            refurl = cols[12].findAll('a')[0]['href'].strip()
            source = get_source(name, reference = reference, url = refurl)
            secondarysource = get_source(name, reference = secondaryreference, url = secondaryrefurl, secondary = True)
            sources = ','.join(list(filter(None, [source, secondarysource])))
            if str(cols[1].contents[0]).strip() != 'unk':
                add_quanta(name, 'claimedtype', str(cols[1].contents[0]).strip(), sources)
            if str(cols[2].contents[0]).strip() != 'anonymous':
                add_quanta(name, 'host', str(cols[2].contents[0]).strip(), sources)
            add_quanta(name, 'snra', str(cols[3].contents[0]).strip(), sources)
            add_quanta(name, 'sndec', str(cols[4].contents[0]).strip(), sources)
            if str(cols[6].contents[0]).strip() not in ['2440587', '2440587.292']:
                astrot = astrotime(float(str(cols[6].contents[0]).strip()), format='jd')
                add_quanta(name, 'discoverday', str(astrot.datetime.day), sources)
                add_quanta(name, 'discovermonth', str(astrot.datetime.month), sources)
                add_quanta(name, 'discoveryear', str(astrot.datetime.year), sources)
            if str(cols[7].contents[0]).strip() not in ['2440587', '2440587.292']:
                astrot = astrotime(float(str(cols[7].contents[0]).strip()), format='jd')
                if float(str(cols[8].contents[0]).strip()) <= 90.0:
                    add_photometry(name, time = str(astrot.mjd), magnitude = str(cols[8].contents[0]).strip(), source = sources)
            if cols[11].contents[0] != 'n/a':
                add_quanta(name, 'redshift', str(cols[11].contents[0]).strip(), sources)
            add_quanta(name, 'discoverer', str(cols[13].contents[0]).strip(), sources)

    if not args.update:
        vsnetfiles = ["latestsne.dat"]
        for vsnetfile in vsnetfiles:
            f = open("../sne-external/" + vsnetfile,'r',encoding='latin1')
            tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)
            for r, row in enumerate(tsvin):
                if not row or row[0][:4] in ['http', 'www.'] or len(row) < 3:
                    continue
                name = row[0].strip()
                if name[:4].isdigit():
                    name = 'SN' + name
                if name[:4] == 'PSNJ':
                    name = 'PSN J' + name[4:]
                name = add_event(name)
                if not is_number(row[1]):
                    continue
                year = row[1][:4]
                month = row[1][4:6]
                day = row[1][6:]
                if '.' not in day:
                    day = day[:2] + '.' + day[2:]
                mjd = astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day))
                magnitude = row[2].rstrip(ascii_letters)
                if not is_number(magnitude):
                    continue
                if magnitude.isdigit():
                    if int(magnitude) > 100:
                        magnitude = magnitude[:2] + '.' + magnitude[2:]
                secondarysource = get_source(name, reference = secondaryreference, url = secondaryrefurl, secondary = True)
                band = row[2].lstrip('1234567890.')
                if len(row) >= 4:
                    if is_number(row[3]):
                        e_magnitude = row[3]
                        refind = 4
                    else:
                        e_magnitude = ''
                        refind = 3

                    if refind >= len(row):
                        sources = secondarysource
                    else:
                        reference = ' '.join(row[refind:])
                        source = get_source(name, reference = reference)
                        sources = ','.join([source,secondarysource])
                else:
                    sources = secondarysource
                add_photometry(name, time = mjd, band = band, magnitude = magnitude, e_magnitude = e_magnitude, source = sources)
            f.close()
    journal_events()

if do_task('lennarz'): 
    Vizier.ROW_LIMIT = -1
    result = Vizier.get_catalogs("J/A+A/538/A120/usc")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)

    bibcode = "2012A&A...538A.120L"
    for row in table:
        row = convert_aq_output(row)
        name = 'SN' + row['SN']
        name = add_event(name)

        source = get_source(name, bibcode = bibcode)

        if row['Gal']:
            add_quanta(name, 'host', row['Gal'], source)
        if row['z']:
            if name != 'SN1985D':
                add_quanta(name, 'redshift', row['z'], source)
        if row['Dist']:
            add_quanta(name, 'lumdist', row['Dist'], source)

        if row['Ddate']:
            dateparts = row['Ddate'].split('-')
            if len(dateparts) == 3:
                astrot = astrotime(row['Ddate'], scale='utc')
            elif len(dateparts) == 2:
                astrot = astrotime(row['Ddate'] + '-01', scale='utc')
            else:
                astrot = astrotime(row['Ddate'] + '-01-01', scale='utc')

            if 'photometry' not in events[name]:
                if 'Dmag' in row and is_number(row['Dmag']) and not isnan(float(row['Dmag'])):
                    mjd = str(astrot.mjd)
                    add_photometry(name, time = mjd, band = row['Dband'], magnitude = row['Dmag'], source = source)
            if 'discoveryear' not in events[name] and 'discovermonth' not in events[name] and 'discoverday' not in events[name]:
                add_quanta(name, 'discoveryear', str(astrot.datetime.year), source)
                if len(dateparts) >= 2:
                    add_quanta(name, 'discovermonth', str(astrot.datetime.month), source)
                if len(dateparts) == 3:
                    add_quanta(name, 'discoverday', str(astrot.datetime.day), source)
        if row['Mdate']:
            dateparts = row['Mdate'].split('-')
            if len(dateparts) == 3:
                astrot = astrotime(row['Mdate'], scale='utc')
            elif len(dateparts) == 2:
                astrot = astrotime(row['Mdate'] + '-01', scale='utc')
            else:
                astrot = astrotime(row['Mdate'] + '-01-01', scale='utc')

            if 'photometry' not in events[name]:
                if 'MMag' in row and is_number(row['MMag']) and not isnan(float(row['MMag'])):
                    mjd = str(astrot.mjd)
                    add_photometry(name, time = mjd, band = row['Mband'], magnitude = row['Mmag'], source = source)
            if 'maxyear' not in events[name] and 'maxmonth' not in events[name] and 'maxday' not in events[name]:
                add_quanta(name, 'maxyear', str(astrot.datetime.year), source)
                if len(dateparts) >= 2:
                    add_quanta(name, 'maxmonth', str(astrot.datetime.month), source)
                if len(dateparts) == 3:
                    add_quanta(name, 'maxday', str(astrot.datetime.day), source)
    f.close()
    journal_events()

if do_task('ogle'): 
    basenames = ['transients', 'transients/2014b', 'transients/2014', 'transients/2013', 'transients/2012']
    oglenames = []
    ogleupdate = [True, False, False, False, False]
    for b, bn in enumerate(basenames):
        if args.update and not ogleupdate[b]:
            continue
        response = urllib.request.urlopen('http://ogle.astrouw.edu.pl/ogle4/' + bn + '/transients.html')
        soup = BeautifulSoup(response.read(), "html5lib")
        links = soup.findAll('a')
        breaks = soup.findAll('br')
        datalinks = []
        for a in links:
            if a.has_attr('href'):
                if '.dat' in a['href']:
                    datalinks.append('http://ogle.astrouw.edu.pl/ogle4/' + bn + '/' + a['href'])

        ec = 0
        reference = 'OGLE-IV Transient Detection System'
        refurl = 'http://ogle.astrouw.edu.pl/ogle4/transients/transients.html'
        for br in breaks:
            sibling = br.nextSibling
            if 'Ra,Dec=' in sibling:
                line = sibling.replace("\n", '').split('Ra,Dec=')
                name = line[0].strip()

                if 'NOVA' in name or 'dupl' in name:
                    continue

                if name in oglenames:
                    continue
                oglenames.append(name)

                name = add_event(name)

                mySibling = sibling.nextSibling
                atelref = ''
                claimedtype = ''
                while 'Ra,Dec=' not in mySibling:
                    if isinstance(mySibling, NavigableString):
                        if 'Phot.class=' in str(mySibling):
                            claimedtype = re.sub(r'\([^)]*\)', '', str(mySibling).split('=')[-1]).replace('SN','').strip()
                    if isinstance(mySibling, Tag):
                        atela = mySibling
                        if atela and atela.has_attr('href') and 'astronomerstelegram' in atela['href']:
                            atelref = a.contents[0].strip()
                            atelurl = a['href']
                    mySibling = mySibling.nextSibling
                    if mySibling is None:
                        break

                nextSibling = sibling.nextSibling
                if isinstance(nextSibling, Tag) and nextSibling.has_attr('alt') and nextSibling.contents[0].strip() != 'NED':
                    radec = nextSibling.contents[0].strip().split()
                else:
                    radec = line[-1].split()
                ra = radec[0]
                dec = radec[1]
                lcresponse = urllib.request.urlopen(datalinks[ec])
                lcdat = lcresponse.read().decode('utf-8').splitlines()
                sources = [get_source(name, reference = reference, url = refurl)]
                if atelref and atelref != 'ATel#----':
                    sources.append(get_source(name, reference = atelref, url = atelurl))
                sources = ','.join(sources)

                if name[:4] == 'OGLE':
                    if name[4] == '-':
                        if is_number(name[5:9]):
                            add_quanta(name, 'discoveryear', name[5:9], sources)
                    else:
                        if is_number(name[4:6]):
                            add_quanta(name, 'discoveryear', '20' + name[4:6], sources)

                add_quanta(name, 'snra', ra, sources)
                add_quanta(name, 'sndec', dec, sources)
                if claimedtype and claimedtype != '-':
                    add_quanta(name, 'claimedtype', claimedtype, sources)
                elif 'SN' not in name and 'claimedtype' not in events[name]:
                    add_quanta(name, 'claimedtype', 'Candidate', sources)
                for row in lcdat:
                    row = row.split()
                    mjd = str(jd_to_mjd(Decimal(row[0])))
                    magnitude = row[1]
                    if float(magnitude) > 90.0:
                        continue
                    e_magnitude = row[2]
                    upperlimit = False
                    if e_magnitude == '-1' or float(e_magnitude) > 10.0:
                        e_magnitude = ''
                        upperlimit = True
                    add_photometry(name, time = mjd, band = 'I', magnitude = magnitude, e_magnitude = e_magnitude, source = sources, upperlimit = upperlimit)
                ec += 1
        journal_events()

if do_task('snls'): 
    with open("../sne-external/SNLS-ugriz.dat", 'r') as f:
        data = csv.reader(f, delimiter=' ', quotechar='"', skipinitialspace = True)
        for row in data:
            flux = row[3]
            err = row[4]
            if float(flux) < float(err):
                continue
            name = 'SNLS-' + row[0]
            name = add_event(name)
            source = get_source(name, bibcode = '2010A&A...523A...7G')
            band = row[1]
            mjd = row[2]
            sig = get_sig_digits(flux.split('E')[0])
            # Conversion comes from SNLS-Readme
            magnitude = pretty_num(30.0-2.5*log10(float(flux)), sig = sig)
            e_magnitude = pretty_num(2.5*(log10(float(flux) + float(err)) - log10(float(flux))), sig = sig)
            add_photometry(name, time = mjd, band = band, magnitude = magnitude, e_magnitude = e_magnitude, source = source)

if do_task('nedd'): 
    f = open("../sne-external/NED25.12.1-D-10.4.0-20151123.csv", 'r')
    data = csv.reader(f, delimiter=',', quotechar='"')
    reference = "NED-D"
    refurl = "http://ned.ipac.caltech.edu/Library/Distances/"
    oldhostname = ''
    for r, row in enumerate(data):
        if r <= 12:
            continue
        hostname = row[3]
        #if hostname == oldhostname:
        #    continue
        distmod = row[4]
        moderr = row[5]
        dist = row[6]
        disterr = ''
        if moderr:
            sig = get_sig_digits(moderr)
            disterr = pretty_num(1.0e-6*(10.0**(0.2*(5.0 + float(distmod))) * (10.0**(0.2*float(moderr)) - 1.0)), sig = sig)
        bibcode = row[8]
        name = ''
        if hostname[:3] == 'SN ':
            if is_number(hostname[3:7]):
                name = 'SN' + hostname[3:]
            else:
                name = hostname[3:]
        if hostname[:5] == 'SNLS ':
            name = 'SNLS-' + hostname[5:].split()[0]
        if name:
            name = add_event(name)
            secondarysource = get_source(name, reference = reference, url = refurl, secondary = True)
            if bibcode:
                source = get_source(name, bibcode = bibcode)
                sources = ','.join([source, secondarysource])
            else:
                sources = secondarysource
            add_quanta(name, 'lumdist', dist, sources, error = disterr)
        #else:
        #    cleanhost = hostname.replace('MESSIER 0', 'M').replace('MESSIER ', 'M').strip()
        #    for name in events:
        #        if 'host' in events[name]:
        #            for host in events[name]['host']:
        #                if host['value'] == cleanhost:
        #                    print ('found host: ' + cleanhost)
        #                    secondarysource = get_source(name, reference = reference, url = refurl, secondary = True)
        #                    if bibcode:
        #                        source = get_source(name, bibcode = bibcode)
        #                        sources = ','.join([source, secondarysource])
        #                    else:
        #                        sources = secondarysource
        #                    add_quanta(name, 'lumdist', dist, sources, error = disterr)
        #                    break
        oldhostname = hostname
    journal_events()

if do_task('cfaiaspectra'): 
    oldname = ''
    for name in sorted(next(os.walk("../sne-external-spectra/CfA_SNIa"))[1], key=lambda s: s.lower()):
        fullpath = "../sne-external-spectra/CfA_SNIa/" + name
        if name[:2] == 'sn' and is_number(name[2:6]):
            name = 'SN' + name[2:]
        if name[:3] == 'snf' and is_number(name[3:7]):
            name = 'SNF' + name[3:]
        if name != oldname:
            clear_events()
        oldname = name
        name = add_event(name)
        reference = 'CfA Supernova Archive'
        refurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
        source = get_source(name, reference = reference, url = refurl, secondary = True)
        for file in sorted(glob.glob(fullpath + '/*'), key=lambda s: s.lower()):
            fileparts = os.path.basename(file).split('-')
            if name[:2] == "SN" and is_number(name[2:6]):
                year = fileparts[1][:4]
                month = fileparts[1][4:6]
                day = fileparts[1][6:]
                instrument = fileparts[2].split('.')[0]
            else:
                year = fileparts[2][:4]
                month = fileparts[2][4:6]
                day = fileparts[2][6:]
                instrument = fileparts[3].split('.')[0]
            time = astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day))
            f = open(file,'r')
            data = csv.reader(f, delimiter=' ', skipinitialspace=True)
            data = [list(i) for i in zip(*data)]
            wavelengths = data[0]
            fluxes = data[1]
            errors = data[2]
            add_spectrum(name = name, waveunit = 'Angstrom', fluxunit = 'erg/s/cm^2/Angstrom',
                wavelengths = wavelengths, fluxes = fluxes, timeunit = 'MJD', time = time, instrument = instrument,
                errorunit = "ergs/s/cm^2/Angstrom", errors = errors, source = source, dereddened = False, deredshifted = False)
        journal_events(clear = False)

if do_task('cfaibcspectra'): 
    oldname = ''
    for name in sorted(next(os.walk("../sne-external-spectra/CfA_SNIbc"))[1], key=lambda s: s.lower()):
        fullpath = "../sne-external-spectra/CfA_SNIbc/" + name
        if name[:2] == 'sn' and is_number(name[2:6]):
            name = 'SN' + name[2:]
        if name != oldname:
            clear_events()
        oldname = name
        name = add_event(name)
        reference = 'CfA Supernova Archive'
        refurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
        source = get_source(name, reference = reference, url = refurl, secondary = True)
        for file in sorted(glob.glob(fullpath + '/*'), key=lambda s: s.lower()):
            fileparts = os.path.basename(file).split('-')
            instrument = ''
            year = fileparts[1][:4]
            month = fileparts[1][4:6]
            day = fileparts[1][6:].split('.')[0]
            if len(fileparts) > 2:
                instrument = fileparts[-1].split('.')[0]
            time = astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day))
            f = open(file,'r')
            data = csv.reader(f, delimiter=' ', skipinitialspace=True)
            data = [list(i) for i in zip(*data)]
            wavelengths = data[0]
            fluxes = data[1]
            add_spectrum(name = name, waveunit = 'Angstrom', fluxunit = 'Uncalibrated', wavelengths = wavelengths,
                fluxes = fluxes, timeunit = 'MJD', time = time, instrument = instrument, source = source,
                dereddened = False, deredshifted = False)
        journal_events(clear = False)

if do_task('snlsspectra'): 
    result = Vizier.get_catalogs("J/A+A/507/85/table1")
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    datedict = {}
    for row in table:
        datedict['SNLS-' + row['SN']] = str(astrotime(row['Date']).mjd)

    oldname = ''
    for file in sorted(glob.glob('../sne-external-spectra/SNLS/*'), key=lambda s: s.lower()):
        fileparts = os.path.basename(file).split('_')
        name = 'SNLS-' + fileparts[1]
        if name != oldname:
            clear_events()
        oldname = name
        name = add_event(name)
        source = get_source(name, bibcode = "2009A&A...507...85B")

        add_quanta(name, 'discoveryear', '20' + fileparts[1][:2], source)

        f = open(file,'r')
        data = csv.reader(f, delimiter=' ', skipinitialspace=True)
        specdata = []
        for r, row in enumerate(data):
            if row[0] == '@TELESCOPE':
                instrument = row[1].strip()
            elif row[0] == '@REDSHIFT':
                add_quanta(name, 'redshift', row[1].strip(), source)
            if r < 14:
                continue
            specdata.append(list(filter(None, [x.strip(' \t') for x in row])))
        specdata = [list(i) for i in zip(*specdata)]
        wavelengths = specdata[1]
        
        fluxes = [pretty_num(float(x)*1.e-16, sig = get_sig_digits(x)) for x in specdata[2]]
        errors = [pretty_num(float(x)*1.e-16, sig = get_sig_digits(x)) for x in specdata[3]]

        add_spectrum(name = name, waveunit = 'Angstrom', fluxunit = 'erg/s/cm^2/Angstrom', wavelengths = wavelengths,
            fluxes = fluxes, timeunit = 'MJD', time = datedict[name], instrument = instrument, source = source)
        journal_events(clear = False)

if do_task('cspspectra'): 
    oldname = ''
    for file in sorted(glob.glob('../sne-external-spectra/CSP/*'), key=lambda s: s.lower()):
        sfile = os.path.basename(file).split('.')
        if sfile[1] == 'txt':
            continue
        sfile = sfile[0]
        fileparts = sfile.split('_')
        name = 'SN20' + fileparts[0][2:]
        if name != oldname:
            clear_events()
        oldname = name
        name = add_event(name)
        instrument = ': '.join(fileparts[-2:])
        source = get_source(name, bibcode = "2013ApJ...773...53F")

        f = open(file,'r')
        data = csv.reader(f, delimiter=' ', skipinitialspace=True)
        specdata = []
        for r, row in enumerate(data):
            if row[0] == '#JDate_of_observation:':
                jd = row[1].strip()
                time = str(jd_to_mjd(Decimal(jd)))
            elif row[0] == '#Redshift:':
                add_quanta(name, 'redshift', row[1].strip(), source)
            if r < 7:
                continue
            specdata.append(list(filter(None, [x.strip(' ') for x in row])))
        specdata = [list(i) for i in zip(*specdata)]
        wavelengths = specdata[0]
        fluxes = specdata[1]

        add_spectrum(name = name, timeunit = 'MJD', time = time, waveunit = 'Angstrom', fluxunit = 'erg/s/cm^2/Angstrom', wavelengths = wavelengths,
            fluxes = fluxes, instrument = instrument, source = source, deredshifted = True)
        journal_events(clear = False)

if do_task('ucbspectra'): 
    secondaryreference = "UCB Filippenko Group's Supernova Database (SNDB)"
    secondaryrefurl = "http://heracles.astro.berkeley.edu/sndb/info"

    path = os.path.abspath('../sne-external-spectra/UCB/sndb.html')
    response = urllib.request.urlopen('file://' + path)

    soup = BeautifulSoup(response.read(), "html5lib")
    i = 0
    oldname = ''
    for t, tr in enumerate(soup.findAll('tr')):
        if t == 0:
            continue
        for d, td in enumerate(tr.findAll('td')):
            if d == 2:
                claimedtype = td.contents[0].strip()
            elif d == 4:
                filename = td.contents[0].strip()
                name = filename.split('-')[0]
                if name[:2].upper() == 'SN':
                    name = name[:2].upper() + name[2:]
                    if len(name) == 7:
                        name = name[:6] + name[6].upper()
                if name[:3] == 'ptf':
                    name = name[:3].upper() + name[3:]
            elif d == 5:
                epoch = td.contents[0].strip()
                year = epoch[:4]
                month = epoch[4:6]
                day = epoch[6:]
                sig = get_sig_digits(day) + 5
                mjd = pretty_num(astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day)), sig = sig)
            elif d == 7:
                instrument = '' if td.contents[0].strip() == 'None' else td.contents[0].strip()
            elif d == 9:
                snr = td.contents[0].strip()
            elif d == 10:
                observerreducer = td.contents[0].strip().split('|')
                observer = '' if observerreducer[0].strip() == 'None' else observerreducer[0].strip()
                reducer = '' if observerreducer[1].strip() == 'None' else observerreducer[1].strip()
            elif d == 11:
                bibcode = td.findAll('a')[0].contents[0]

        if name != oldname:
            clear_events()
        oldname = name
        name = add_event(name)
        source = get_source(name, bibcode = bibcode)
        secondarysource = get_source(name, reference = secondaryreference, url = secondaryrefurl, secondary = True)
        sources = ','.join([source, secondarysource])
        add_quanta(name, 'claimedtype', claimedtype, sources)
        if 'discoveryear' not in events[name] and name[:2] == 'SN' and is_number(name[2:6]):
            add_quanta(name, 'discoveryear', name[2:6], sources)

        with open('../sne-external-spectra/UCB/' + filename) as f:
            specdata = list(csv.reader(f, delimiter=' ', skipinitialspace=True))
            startrow = 0
            for row in specdata:
                if row[0][0] == '#':
                    startrow += 1
                else:
                    break
            specdata = specdata[startrow:]

            haserrors = len(specdata[0]) == 3 and specdata[0][2] and specdata[0][2] != 'NaN'
            specdata = [list(i) for i in zip(*specdata)]

            wavelengths = specdata[0]
            fluxes = specdata[1]
            errors = ''
            if haserrors:
                errors = specdata[2]

            if not list(filter(None, errors)):
                errors = ''

            add_spectrum(name = name, timeunit = 'MJD', time = mjd, waveunit = 'Angstrom', fluxunit = 'Uncalibrated', wavelengths = wavelengths,
                fluxes = fluxes, errors = errors, errorunit = 'Uncalibrated', instrument = instrument, source = source, snr = snr, observer = observer, reducer = reducer,
                deredshifted = True)
        journal_events(clear = False)

if do_task('suspectspectra'): 
    with open('../sne-external-spectra/Suspect/sources.json', 'r') as f:
        sourcedict = json.loads(f.read())

    folders = next(os.walk('../sne-external-spectra/Suspect'))[1]
    for folder in folders:
        eventfolders = next(os.walk('../sne-external-spectra/Suspect/'+folder))[1]
        oldname = ''
        for eventfolder in eventfolders:
            name = eventfolder
            if is_number(name[:4]):
                name = 'SN' + name
            if name != oldname:
                clear_events()
            oldname = name
            name = add_event(name)
            secondaryreference = "SUSPECT"
            secondaryrefurl = "https://www.nhn.ou.edu/~suspect/"
            secondarysource = get_source(name, reference = secondaryreference, url = secondaryrefurl, secondary = True)
            eventspectra = next(os.walk('../sne-external-spectra/Suspect/'+folder+'/'+eventfolder))[2]
            for spectrum in eventspectra:
                sources = [secondarysource]
                bibcode = ''
                if spectrum in sourcedict:
                    bibcode = sourcedict[spectrum]
                elif name in sourcedict:
                    bibcode = sourcedict[name]
                if bibcode:
                    source = get_source(name, bibcode = bibcode)
                    sources += source
                sources = ','.join(sources)

                date = spectrum.split('_')[1]
                year = date[:4]
                month = date[4:6]
                day = date[6:]
                sig = get_sig_digits(day) + 5
                time = pretty_num(astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day)), sig = sig)

                with open('../sne-external-spectra/Suspect/'+folder+'/'+eventfolder+'/'+spectrum) as f:
                    specdata = list(csv.reader(f, delimiter=' ', skipinitialspace=True))
                    specdata = list(filter(None, specdata))
                    newspec = []
                    oldval = ''
                    for row in specdata:
                        if row[1] == oldval:
                            continue
                        newspec.append(row)
                        oldval = row[1]
                    specdata = newspec
                haserrors = len(specdata[0]) == 3 and specdata[0][2] and specdata[0][2] != 'NaN'
                specdata = [list(i) for i in zip(*specdata)]

                wavelengths = specdata[0]
                fluxes = specdata[1]
                errors = ''
                if haserrors:
                    errors = specdata[2]

                add_spectrum(name = name, timeunit = 'MJD', time = time, waveunit = 'Angstrom', fluxunit = 'Uncalibrated', wavelengths = wavelengths,
                    fluxes = fluxes, errors = errors, errorunit = 'Uncalibrated', source = sources)
            journal_events(clear = False)

if do_task('snfspectra'): 
    eventfolders = next(os.walk('../sne-external-spectra/SNFactory'))[1]
    bibcodes = {'SN2005gj':'2006ApJ...650..510A', 'SN2006D':'2007ApJ...654L..53T', 'SN2007if':'2010ApJ...713.1073S', 'SN2011fe':'2013A&A...554A..27P'}
    oldname = ''
    for eventfolder in eventfolders:
        name = eventfolder
        if name != oldname:
            clear_events()
        oldname = name
        name = add_event(name)
        secondaryreference = "Nearby Supernova Factory"
        secondaryrefurl = "http://snfactory.lbl.gov/"
        secondarysource = get_source(name, reference = secondaryreference, url = secondaryrefurl, secondary = True)
        bibcode = bibcodes[name]
        source = get_source(name, bibcode = bibcode)
        sources = ','.join([source,secondarysource])
        eventspectra = glob.glob('../sne-external-spectra/SNFactory/'+eventfolder+'/*.dat')
        for spectrum in eventspectra:
            with open(spectrum) as f:
                specdata = list(csv.reader(f, delimiter=' ', skipinitialspace=True))
            specdata = list(filter(None, specdata))
            newspec = []
            time = ''
            if 'Keck_20060202_R' in spectrum:
                time = '53768.23469'
            elif 'Spectrum05_276' in spectrum:
                time = pretty_num(astrotime('2005-10-03').mjd, sig = 5)
            elif 'Spectrum05_329' in spectrum:
                time = pretty_num(astrotime('2005-11-25').mjd, sig = 5)
            elif 'Spectrum05_336' in spectrum:
                time = pretty_num(astrotime('2005-12-02').mjd, sig = 5)
            for row in specdata:
                if not time:
                    if row[0] == '#MJD-OBS':
                        time = row[2].strip("'")
                    elif len(row) >= 2:
                        if row[1] == 'JD':
                            time = str(jd_to_mjd(Decimal(row[3])))
                        elif row[1] == 'MJD':
                            time = row[3]
                        elif row[1] == 'MJD-OBS':
                            time = row[3].strip("'")
                if row[0][0] == '#':
                    continue
                newspec.append(row)
            if not time:
                print(spectrum)
                sys.exit()
            specdata = newspec
            haserrors = len(specdata[0]) == 3 and specdata[0][2] and specdata[0][2] != 'NaN'
            specdata = [list(i) for i in zip(*specdata)]

            wavelengths = specdata[0]
            fluxes = specdata[1]
            errors = ''
            if haserrors:
                errors = specdata[2]

            add_spectrum(name = name, timeunit = 'MJD', time = time, waveunit = 'Angstrom', fluxunit = 'erg/s/cm^2/Angstrom', wavelengths = wavelengths,
                fluxes = fluxes, errors = errors, errorunit = ('Variance' if name == 'SN2011fe' else 'erg/s/cm^2/Angstrom'), source = sources)
        journal_events(clear = False)

if do_task('writeevents'): 
    files = []
    for rep in repfolders:
        files += glob.glob('../' + rep + "/*.json")

    for fi in files:
        events = OrderedDict()
        name = os.path.basename(os.path.splitext(fi)[0])
        name = add_event(name)
        derive_and_sanitize()
        write_all_events(empty = True)

print("Memory used (MBs on Mac, GBs on Linux): " + "{:,}".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024./1024.))

