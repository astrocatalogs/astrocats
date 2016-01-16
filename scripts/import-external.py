#!/usr/local/bin/python3.5

import csv
import glob
import os
import re
import urllib.request, urllib.error, urllib.parse
import calendar
import sys
import subprocess
import json
import codecs
import numpy
from cdecimal import Decimal
from astroquery.vizier import Vizier
from astropy.time import Time as astrotime
from astropy.cosmology import Planck15 as cosmo
from collections import OrderedDict
from math import log10, floor, sqrt, isnan
from bs4 import BeautifulSoup, SoupStrainer
from operator import itemgetter
from string import ascii_letters

clight = 29979245800.

eventnames = []

dovizier =        True
dosuspect =       True
docfa =           True
doucb =           True
dosdss =          True
dogaia =          True
docsp =           True
doitep =          True
doasiago =        True
dorochester =     True
dofirstmax =      True
dolennarz =       True
docfaiaspectra =  True
docfaibcspectra = True
writeevents =     True
printextra =      False

events = OrderedDict()

repfolders = [
    'sne-pre-1990',
    'sne-1990-1999',
    'sne-2000-2004',
    'sne-2005-2009',
    'sne-2010-2014',
    'sne-2015-2019'
]

repyears = [int(repfolders[x][-4:]) for x in range(len(repfolders))]

typereps = {
    'II P':  ['II pec', 'IIpec', 'II Pec', 'IIPec', 'IIP', 'IIp', 'II p', 'II-pec', 'II P pec'],
    'Ib/c':  ['Ibc'],
    'Ia P':  ['Ia pec', 'Ia-pec', 'Iapec'],
    'Ic P':  ['Ic pec', 'Ic-pec'],
    'IIn P': ['IIn pec', 'IIn-pec'],
    'II L':  ['IIL'],
    'I P':   ['I pec', 'I-pec', 'I Pec', 'I-Pec'],
    'Ib P':  ['Ib pec', 'Ib-pec']
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

def add_event(name):
    if name not in events:
        for event in events:
            if len(events[event]['aliases']) > 1 and name in events[event]['aliases']:
                return event
        print(name)
        events[name] = OrderedDict()
        events[name]['name'] = name
        add_alias(name, name)
        events[name]['sources'] = []
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

def round_sig(x, sig=4):
    if x == 0.0:
        return 0.0
    return round(x, sig-int(floor(log10(abs(x))))-1)

def pretty_num(x, sig=4):
    return str('%g'%(round_sig(x, sig)))

def get_source(name, reference, secondary = ''):
    nsources = len(events[name]['sources'])
    if len(events[name]['sources']) == 0 or reference not in [events[name]['sources'][x]['name'] for x in range(nsources)]:
        source = str(nsources + 1)
        newsource = OrderedDict()
        newsource['name'] = reference
        newsource['alias'] =  source
        if secondary:
            newsource['secondary'] = True
        events[name].setdefault('sources',[]).append(newsource)
    else:
        sourcexs = range(nsources)
        source = [events[name]['sources'][x]['alias'] for x in sourcexs][
            [events[name]['sources'][x]['name'] for x in sourcexs].index(reference)]
    return source

def add_photometry(name, timeunit = "MJD", time = "", instrument = "", band = "", abmag = "", aberr = "", source = ""):
    if not time or not abmag:
        print('Warning: Time or AB mag not specified when adding photometry.\n')
        print('Name : "' + name + '", Time: "' + time + '", Band: "' + band + '", AB mag: "' + abmag + '"')
        return

    if not is_number(time) or not is_number(abmag):
        print('Warning: Time or AB mag not numerical.\n')
        print('Name : "' + name + '", Time: "' + time + '", Band: "' + band + '", AB mag: "' + abmag + '"')
        return

    if aberr and not is_number(aberr):
        print('Warning: AB error not numerical.\n')
        print('Name : "' + name + '", Time: "' + time + '", Band: "' + band + '", AB err: "' + aberr + '"')
        return

    # Look for duplicate data and don't add if duplicate
    if 'photometry' in events[name]:
        for photo in events[name]['photometry']:
            if (photo['timeunit'] == timeunit and photo['band'] == band and
                Decimal(photo['time']) == Decimal(time) and
                Decimal(photo['abmag']) == Decimal(abmag) and
                (('aberr' not in photo and not aberr) or ('aberr' in photo and aberr and Decimal(photo['aberr']) == Decimal(aberr)) or
                ('aberr' in photo and not aberr))):
                return

    photoentry = OrderedDict()
    photoentry['timeunit'] = timeunit
    photoentry['time'] = str(time)
    photoentry['band'] = band
    photoentry['abmag'] = str(abmag)
    if instrument:
        photoentry['instrument'] = instrument
    if aberr:
        photoentry['aberr'] = str(aberr)
    if source:
        photoentry['source'] = source
    events[name].setdefault('photometry',[]).append(photoentry)

def add_spectrum(name, waveunit, fluxunit, wavelengths, fluxes, timeunit, time, instrument = "",
    deredshifted = False, dereddened = False, errorunit = "", errors = "", source = ""):
    if not waveunit:
        'Warning: No error unit specified, not adding spectrum.'
        return
    if not fluxunit:
        'Warning: No flux unit specified, not adding spectrum.'
        return
    spectrumentry = OrderedDict()
    spectrumentry['deredshifted'] = deredshifted
    spectrumentry['dereddened'] = dereddened
    spectrumentry['instrument'] = instrument
    spectrumentry['timeunit'] = timeunit
    spectrumentry['time'] = time
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

def add_claimed_type(name, inputtype, source = ""):
    claimedtype = inputtype.strip()
    for rep in typereps:
        if inputtype in typereps[rep]:
            claimedtype = rep
            break

    if 'claimedtype' in events[name]:
        for i, ct in enumerate(events[name]['claimedtype']):
            if ct['type'] == claimedtype:
                if source:
                    events[name]['claimedtype'][i]['source'] += ',' + source
                return

    claimedtypeentry = OrderedDict()
    claimedtypeentry['type'] = claimedtype
    if source:
        claimedtypeentry['source'] = source
    events[name].setdefault('claimedtype',[]).append(claimedtypeentry)

def get_max_light(name):
    if 'photometry' not in events[name]:
        return (None, None)

    eventphoto = [Decimal(events[name]['photometry'][x]['abmag']) for x in range(len(events[name]['photometry']))]
    mlmag = min(eventphoto)
    mlindex = eventphoto.index(mlmag)
    mlmjd = float(events[name]['photometry'][mlindex]['time'])
    return (astrotime(mlmjd, format='mjd').datetime, mlmag)

def get_first_light(name):
    if 'photometry' not in events[name]:
        return None

    eventtime = [events[name]['photometry'][x]['time'] for x in range(len(events[name]['photometry']))]
    flindex = eventtime.index(min(eventtime))
    flmjd = float(events[name]['photometry'][flindex]['time'])
    return astrotime(flmjd, format='mjd').datetime

def set_first_max_light(name):
    (mldt, mlmag) = get_max_light(name)
    if mldt:
        events[name]['maxyear'] = pretty_num(mldt.year)
        events[name]['maxmonth'] = pretty_num(mldt.month)
        events[name]['maxday'] = pretty_num(mldt.day)
        events[name]['maxappmag'] = pretty_num(mlmag)

    fldt = get_first_light(name)
    if fldt:
        events[name]['discoveryear'] = pretty_num(fldt.year)
        events[name]['discovermonth'] = pretty_num(fldt.month)
        events[name]['discoverday'] = pretty_num(fldt.day)

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

catalog = OrderedDict()
def convert_aq_output(row):
    return OrderedDict([(x, '%g'%Decimal(float(row[x])) if is_number(row[x]) else utf8(row[x])) for x in row.colnames])

# Import primary data sources from Vizier
if dovizier:
    Vizier.ROW_LIMIT = -1
    result = Vizier.get_catalogs("VII/272/snrs")
    table = result[list(result.keys())[0]]

    for row in table:
        row = convert_aq_output(row)
        name = ''
        if row["Names"]:
            names = row["Names"].split(',')
            for nam in names:
                if nam.strip()[:2] == 'SN':
                    name = nam.strip()
            if not name:
                for nam in names:
                    if nam.strip('()') == nam:
                        name = nam.strip()
                        break
        if not name:
            name = row["SNR"]

        name = add_event(name)

        if row["Names"]:
            names = row["Names"].split(',')
            for nam in names:
                if nam.strip()[:2] == 'SN':
                    events[name]['discoveryear'] = nam.strip()[2:]

        events[name]['snra'] = row['RAJ2000']
        events[name]['sndec'] = row['DEJ2000']

    reference = "<a href='http://adsabs.harvard.edu/abs/2014MNRAS.442..844F'>2014MNRAS.442..844F</a>"
    result = Vizier.get_catalogs("J/MNRAS/442/844/table1")
    table = result[list(result.keys())[0]]
    for row in table:
        row = convert_aq_output(row)
        name = 'SN' + row['SN']
        name = add_event(name)
        events[name]['redshift'] = row['zhost']
        events[name]['ebv'] = row['E_B-V_']

    result = Vizier.get_catalogs("J/MNRAS/442/844/table2")
    table = result[list(result.keys())[0]]
    for row in table:
        row = convert_aq_output(row)
        name = 'SN' + str(row['SN'])
        source = get_source(name, reference)
        if 'Bmag' in row and is_number(row['Bmag']) and not isnan(float(row['Bmag'])):
            add_photometry(name, time = row['MJD'], band = 'B', abmag = row['Bmag'], aberr = row['e_Bmag'], source = source)
        if 'Vmag' in row and is_number(row['Vmag']) and not isnan(float(row['Vmag'])):
            add_photometry(name, time = row['MJD'], band = 'V', abmag = row['Vmag'], aberr = row['e_Vmag'], source = source)
        if 'Rmag' in row and is_number(row['Rmag']) and not isnan(float(row['Rmag'])):
            add_photometry(name, time = row['MJD'], band = 'R', abmag = row['Rmag'], aberr = row['e_Rmag'], source = source)
        if 'Imag' in row and is_number(row['Imag']) and not isnan(float(row['Imag'])):
            add_photometry(name, time = row['MJD'], band = 'I', abmag = row['Imag'], aberr = row['e_Imag'], source = source)

    reference = "<a href='http://adsabs.harvard.edu/abs/2015ApJS..219...13W'>2015ApJS..219...13W</a>"
    result = Vizier.get_catalogs("J/ApJS/219/13/table3")
    table = result[list(result.keys())[0]]
    for row in table:
        row = convert_aq_output(row)
        name = u'LSQ' + str(row['LSQ'])
        name = add_event(name)
        events[name]['snra'] = row['RAJ2000']
        events[name]['sndec'] = row['DEJ2000']
        events[name]['redshift'] = row['z']
        events[name]['ebv'] = row['E_B-V_']
    result = Vizier.get_catalogs("J/ApJS/219/13/table2")
    table = result[list(result.keys())[0]]
    for row in table:
        row = convert_aq_output(row)
        name = 'LSQ' + row['LSQ']
        source = get_source(name, reference)
        add_photometry(name, time = str(jd_to_mjd(Decimal(row['JD']))), instrument = 'La Silla-QUEST', band = row['Filt'], abmag = row['mag'], aberr = row['e_mag'], source = source)

# Suspect catalog
if dosuspect:
    response = urllib.request.urlopen('http://www.nhn.ou.edu/cgi-bin/cgiwrap/~suspect/snindex.cgi')

    soup = BeautifulSoup(response.read(), "html5lib")
    i = 0
    for a in soup.findAll('a'):
        if 'phot=yes' in a['href'] and not 'spec=yes' in a['href']:
            if int(a.contents[0]) > 0:
                i = i + 1
                photlink = 'http://www.nhn.ou.edu/cgi-bin/cgiwrap/~suspect/' + a['href']
                eventresp = urllib.request.urlopen(photlink)
                eventsoup = BeautifulSoup(eventresp, "html5lib")
                ei = 0
                for ea in eventsoup.findAll('a'):
                    if ea.contents[0] == 'I':
                        ei = ei + 1
                        bandlink = 'http://www.nhn.ou.edu/cgi-bin/cgiwrap/~suspect/' + ea['href']
                        bandresp = urllib.request.urlopen(bandlink)
                        bandsoup = BeautifulSoup(bandresp, "html5lib")
                        bandtable = bandsoup.find('table')
                        if ei == 1:
                            names = bandsoup.body.findAll(text=re.compile("Name"))
                            name = 'SN' + names[0].split(':')[1].strip()
                            name = add_event(name)
                            year = re.findall(r'\d+', name)[0]
                            events[name]['discoveryear'] = year
                            events[name]['host'] = names[1].split(':')[1].strip()
                            redshifts = bandsoup.body.findAll(text=re.compile("Redshift"))
                            if redshifts:
                                events[name]['redshift'] = redshifts[0].split(':')[1].strip()
                            hvels = bandsoup.body.findAll(text=re.compile("Heliocentric Velocity"))
                            if hvels:
                                events[name]['hvel'] = hvels[0].split(':')[1].strip().split(' ')[0]
                            types = bandsoup.body.findAll(text=re.compile("Type"))

                            reference = ''
                            for link in bandsoup.body.findAll('a'):
                                if 'adsabs' in link['href']:
                                    reference = str(link).replace('"', "'")
                            source = get_source(name, reference)

                            add_claimed_type(name, types[0].split(':')[1].strip().split(' ')[0], source)

                        bands = bandsoup.body.findAll(text=re.compile("^Band"))
                        band = bands[0].split(':')[1].strip()

                        secondaryreference = "<a href='https://www.nhn.ou.edu/~suspect/'>SUSPECT</a>"
                        secondarysource = get_source(name, secondaryreference, 1)

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
                            aberr = col[4].contents[0]
                            if aberr.isspace():
                                aberr = ''
                            else:
                                aberr = str(aberr)
                            add_photometry(name, time = mjd, band = band, abmag = mag, aberr = aberr, source = secondarysource + ',' + source)

# CfA data
if docfa:
    for file in sorted(glob.glob("../external/cfa-input/*.dat"), key=lambda s: s.lower()):
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

        year = re.findall(r'\d+', name)[0]
        events[name]['discoveryear'] = year

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
                            print(bibcode)
                            secondaryreference = "<a href='https://www.cfa.harvard.edu/supernova/SNarchive.html'>CfA Supernova Archive</a>"
                            secondarysource = get_source(name, secondaryreference, 1)
                            reference = "<a href='http://adsabs.harvard.edu/abs/" + bibcode + "'>" + refstr + "</a>"
                            source = get_source(name, reference)

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
                            add_photometry(name, timeunit = tuout, time = mjd, band = eventbands[(v-1)//2], abmag = row[v], aberr = row[v+1], source = secondarysource + ',' + source)
        f.close()

    # Hicken 2012
    f = open("../external/hicken-2012-standard.dat", 'r')
    tsvin = csv.reader(f, delimiter='|', skipinitialspace=True)
    for r, row in enumerate(tsvin):
        if r <= 47:
            continue

        if row[0][:2] != 'sn':
            name = 'SN' + row[0].strip()
        else:
            name = row[0].strip()

        name = add_event(name)

        reference = "<a href='http://adsabs.harvard.edu/abs/2012ApJS..200...12H'>Hicken et al. 2012</a>"
        source = get_source(name, reference)
        add_photometry(name, timeunit = 'MJD', time = row[2].strip(), band = row[1].strip(),
            abmag = row[6].strip(), aberr = row[7].strip(), source = source)
    
    # Bianco 2014
    tsvin = open("../external/bianco-2014-standard.dat", 'r')
    tsvin = csv.reader(tsvin, delimiter=' ', skipinitialspace=True)
    for row in tsvin:
        name = 'SN' + row[0]
        name = add_event(name)

        reference = "<a href='http://adsabs.harvard.edu/abs/2014ApJS..213...19B'>Bianco et al. 2014</a>"
        source = get_source(name, reference)
        add_photometry(name, timeunit = 'MJD', time = row[2], band = row[1], abmag = row[3], aberr = row[4], instrument = row[5], source = source)
    f.close()

# Now import the UCB SNDB
if doucb:
    for file in sorted(glob.glob("../external/SNDB/*.dat"), key=lambda s: s.lower()):
        f = open(file,'r')
        tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)

        eventname = os.path.basename(os.path.splitext(file)[0])

        eventparts = eventname.split('.')

        name = snname(eventparts[0])
        name = add_event(name)

        year = re.findall(r'\d+', name)[0]
        events[name]['discoveryear'] = year

        reference = "<a href='http://heracles.astro.berkeley.edu/sndb/info'>UCB Filippenko Group's Supernova Database (SNDB)</a>"
        source = get_source(name, reference)

        for r, row in enumerate(tsvin):
            if len(row) > 0 and row[0] == "#":
                continue
            mjd = row[0]
            abmag = row[1]
            aberr = row[2]
            band = row[4]
            instrument = row[5]
            add_photometry(name, time = mjd, instrument = instrument, band = band, abmag = abmag, aberr = aberr, source = source)
        f.close()
    
# Import SDSS
if dosdss:
    sdssbands = ['u', 'g', 'r', 'i', 'z']
    for file in sorted(glob.glob("../external/SDSS/*.sum"), key=lambda s: s.lower()):
        f = open(file,'r')
        tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)

        for r, row in enumerate(tsvin):
            if r == 0:
                if row[5] == "RA:":
                    name = "SDSS" + row[3]
                else:
                    name = "SN" + row[5]
                name = add_event(name)

                if row[5] != "RA:":
                    year = re.findall(r'\d+', name)[0]
                    events[name]['discoveryear'] = year

                events[name]['snra'] = row[-4]
                events[name]['sndec'] = row[-2]

                reference = "<a href='http://classic.sdss.org/supernova/lightcurves.html'>SDSS Supernova Survey</a>"
                source = get_source(name, reference)
            if r == 1:
                events[name]['redshift'] = row[2]
            if r >= 19:
                # Skip bad measurements
                if int(row[0]) > 1024:
                    continue

                mjd = row[1]
                band = sdssbands[int(row[2])]
                abmag = row[3]
                aberr = row[4]
                instrument = "SDSS"
                add_photometry(name, time = mjd, instrument = instrument, band = band, abmag = abmag, aberr = aberr, source = source)
        f.close()

#Import GAIA
if dogaia:
    #response = urllib2.urlopen('https://gaia.ac.uk/selected-gaia-science-alerts')
    path = os.path.abspath('../external/selected-gaia-science-alerts')
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

        year = '20' + re.findall(r'\d+', name)[0]
        events[name]['discoveryear'] = year

        reference = "<a href='https://gaia.ac.uk/selected-gaia-science-alerts'>Gaia Photometric Science Alerts</a>"
        source = get_source(name, reference)

        events[name]['snra'] = col[2].contents[0].strip()
        events[name]['sndec'] = col[3].contents[0].strip()
        add_claimed_type(name, classname.replace('SN', '').strip(), source)

        photlink = 'http://gsaweb.ast.cam.ac.uk/alerts/alert/' + name + '/lightcurve.csv/'
        photresp = urllib.request.urlopen(photlink)
        photsoup = BeautifulSoup(photresp, "html5lib")
        photodata = str(photsoup.contents[0]).split('\n')[2:-1]
        for ph in photodata:
            photo = ph.split(',')
            mjd = str(jd_to_mjd(Decimal(photo[1].strip())))
            abmag = photo[2].strip()
            aberr = 0.
            instrument = 'GAIA'
            band = 'G'
            add_photometry(name, time = mjd, instrument = instrument, band = band, abmag = abmag, aberr = aberr, source = source)

# Import CSP
if docsp:
    cspbands = ['u', 'B', 'V', 'g', 'r', 'i', 'Y', 'J', 'H', 'K']
    for file in sorted(glob.glob("../external/CSP/*.dat"), key=lambda s: s.lower()):
        f = open(file,'r')
        tsvin = csv.reader(f, delimiter='\t', skipinitialspace=True)

        eventname = os.path.basename(os.path.splitext(file)[0])

        eventparts = eventname.split('opt+')

        name = snname(eventparts[0])
        name = add_event(name)

        year = re.findall(r'\d+', name)[0]
        events[name]['discoveryear'] = year

        reference = "<a href='http://csp.obs.carnegiescience.edu/data'>Carnegie Supernova Project</a>"
        source = get_source(name, reference)

        for r, row in enumerate(tsvin):
            if len(row) > 0 and row[0][0] == "#":
                if r == 2:
                    events[name]['redshift'] = row[0].split(' ')[-1]
                    events[name]['snra'] = row[1].split(' ')[-1]
                    events[name]['sndec'] = row[2].split(' ')[-1]
                continue
            for v, val in enumerate(row):
                if v == 0:
                    mjd = val
                elif v % 2 != 0:
                    if float(row[v]) < 90.0:
                        add_photometry(name, time = mjd, instrument = 'CSP', band = cspbands[(v-1)//2], abmag = row[v], aberr = row[v+1], source = source)
        f.close()

# Import ITEP
if doitep:
    f = open("../external/itep-refs.txt",'r')
    refrep = f.read().splitlines()
    f.close()
    refrepf = dict(list(zip(refrep[1::2], refrep[::2])))
    f = open("../external/itep-lc-cat-28dec2015.txt",'r')
    tsvin = csv.reader(f, delimiter='|', skipinitialspace=True)
    curname = ''
    for r, row in enumerate(tsvin):
        if r <= 1 or len(row) < 7:
            continue
        name = 'SN' + row[0].strip()
        mjd = str(jd_to_mjd(Decimal(row[1].strip())))
        band = row[2].strip()
        abmag = row[3].strip()
        aberr = row[4].strip()
        reference = row[6].strip().strip(',')
        if reference in refrepf:
            reference = refrepf[reference]
        if curname != name:
            curname = name
            name = add_event(name)
            year = re.findall(r'\d+', name)[0]
            events[name]['discoveryear'] = year

            secondaryreference = "<a href='http://dau.itep.ru/sn/node/72'>Sternberg Astronomical Institute Supernova Light Curve Catalogue</a>"
            secondarysource = get_source(name, secondaryreference, 1)

        source = get_source(name, reference) if reference else ''

        add_photometry(name, time = mjd, band = band, abmag = abmag, aberr = aberr, source = secondarysource + ',' + source)
    f.close()


# Now import the Asiago catalog
if doasiago:
    response = urllib.request.urlopen('http://graspa.oapd.inaf.it/cgi-bin/sncat.php')
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

            reference = 'http://graspa.oapd.inaf.it/cgi-bin/sncat.php'
            source = get_source(name, reference, secondary = True)

            year = re.findall(r'\d+', name)[0]
            events[name]['discoveryear'] = year

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
                    events[name][daykey] = dayarr[0]
                monthstr = ''.join(re.findall("[a-zA-Z]+", datestr))
                events[name][monthkey] = list(calendar.month_abbr).index(monthstr)

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
                events[name]['host'] = hostname
            if (claimedtype != ''):
                add_claimed_type(name, claimedtype, source)
            if (redshift != ''):
                events[name]['redshift'] = redshift
            if (hvel != ''):
                events[name]['hvel'] = hvel
            if (galra != ''):
                events[name]['galra'] = galra
            if (galdec != ''):
                events[name]['galdec'] = galdec
            if (snra != ''):
                events[name]['snra'] = snra
            if (sndec != ''):
                events[name]['sndec'] = sndec
            if (discoverer != ''):
                events[name]['discoverer'] = discoverer

if dorochester:
    path = os.path.abspath('../external/snredshiftall.html')
    response = urllib.request.urlopen('file://' + path)
    #response = urllib2.urlopen('http://www.rochesterastronomy.org/snimages/snredshiftall.html')
    html = response.read()

    soup = BeautifulSoup(html, "html5lib")
    rows = soup.findAll('tr')
    secondaryreference = "<a href='http://www.rochesterastronomy.org/snimages/snredshiftall.html'>Latest Supernovae</a>"
    for r, row in enumerate(rows):
        if r == 0:
            continue
        cols = row.findAll('td')
        if not len(cols):
            continue
        name = re.sub('<[^<]+?>', '', str(cols[0].contents[0])).strip()
        if name[:4].isdigit():
            name = 'SN' + name
        name = add_event(name)

        secondarysource = get_source(name, secondaryreference, secondary = 1)
        if str(cols[1].contents[0]).strip() != 'unk':
            add_claimed_type(name, str(cols[1].contents[0]).strip(), secondarysource)
        if str(cols[2].contents[0]).strip() != 'anonymous':
            events[name]['host'] = str(cols[2].contents[0]).strip()
        events[name]['snra'] = str(cols[3].contents[0]).strip()
        events[name]['sndec'] = str(cols[4].contents[0]).strip()
        if str(cols[6].contents[0]).strip() not in ['2440587', '2440587.292']:
            astrot = astrotime(float(str(cols[6].contents[0]).strip()), format='jd')
            events[name]['discoverday'] = str(astrot.datetime.day)
            events[name]['discovermonth'] = str(astrot.datetime.month)
            events[name]['discoveryear'] = str(astrot.datetime.year)
        if str(cols[7].contents[0]).strip() not in ['2440587', '2440587.292']:
            astrot = astrotime(float(str(cols[7].contents[0]).strip()), format='jd')
            source = get_source(name, str(cols[12].contents[0]).strip().replace('"', "'"))
            if float(str(cols[8].contents[0]).strip()) <= 90.0:
                add_photometry(name, time = str(astrot.mjd), abmag = str(cols[8].contents[0]).strip(), source = ','.join([source, secondarysource]))
        if cols[11].contents[0] != 'n/a':
            events[name]['redshift'] = str(cols[11].contents[0]).strip()
        events[name]['discoverer'] = str(cols[13].contents[0]).strip()
        if cols[14].contents:
            add_alias(name, str(cols[14].contents[0]).strip())

    vsnetfiles = ["latestsne.dat"]
    for vsnetfile in vsnetfiles:
        f = open("../external/" + vsnetfile,'r',encoding='latin1')
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
            abmag = row[2].rstrip(ascii_letters)
            if not is_number(abmag):
                continue
            if abmag.isdigit():
                if int(abmag) > 100:
                    abmag = abmag[:2] + '.' + abmag[2:]
            secondaryreference = "<a href='http://www.rochesterastronomy.org/snimages/snredshiftall.html'>Latest Supernovae</a>"
            secondarysource = get_source(name, secondaryreference, secondary = 1)
            band = row[2].lstrip('1234567890.')
            if len(row) >= 4:
                if is_number(row[3]):
                    aberr = row[3]
                    refind = 4
                else:
                    aberr = ''
                    refind = 3
                reference = ' '.join(row[refind:])
                source = get_source(name, reference)
                sources = ','.join([source,secondarysource])
            else:
                sources = secondarysource
            add_photometry(name, time = mjd, band = band, abmag = abmag, aberr = aberr, source = sources)
        f.close()

if dofirstmax:
    for name in events:
        set_first_max_light(name)

if dolennarz:
    Vizier.ROW_LIMIT = -1
    result = Vizier.get_catalogs("J/A+A/538/A120/usc")
    table = result[list(result.keys())[0]]

    reference = "<a href='http://adsabs.harvard.edu/abs/2012A%26A...538A.120L'>2012A&A...538A.120L</a>"
    for row in table:
        row = convert_aq_output(row)
        name = 'SN' + row['SN']
        name = add_event(name)

        source = get_source(name, reference)

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
                    add_photometry(name, time = mjd, band = row['Dband'], abmag = row['Dmag'], source = source)
            if 'discoveryear' not in events[name] and 'discovermonth' not in events[name] and 'discoverday' not in events[name]:
                events[name]['discoveryear'] = str(astrot.datetime.year)
                if len(dateparts) >= 2:
                    events[name]['discovermonth'] = str(astrot.datetime.month)
                if len(dateparts) == 3:
                    events[name]['discoverday'] = str(astrot.datetime.day)
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
                    add_photometry(name, time = mjd, band = row['Mband'], abmag = row['Mmag'], source = source)
            if 'maxyear' not in events[name] and 'maxmonth' not in events[name] and 'maxday' not in events[name]:
                events[name]['maxyear'] = str(astrot.datetime.year)
                if len(dateparts) >= 2:
                    events[name]['maxmonth'] = str(astrot.datetime.month)
                if len(dateparts) == 3:
                    events[name]['maxday'] = str(astrot.datetime.day)

if docfaiaspectra:
    for name in sorted(next(os.walk("../sne-external-spectra/CfA_SNIa"))[1], key=lambda s: s.lower()):
        fullpath = "../sne-external-spectra/CfA_SNIa/" + name
        if name[:2] == 'sn' and is_number(name[2:6]):
            name = 'SN' + name[2:]
        name = add_event(name)
        reference = "<a href='https://www.cfa.harvard.edu/supernova/SNarchive.html'>CfA Supernova Archive</a>"
        source = get_source(name, reference, secondary = True)
        for file in sorted(glob.glob(fullpath + '/*'), key=lambda s: s.lower()):
            fileparts = os.path.basename(file).split('-')
            if name[:2] == "SN":
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
                errorunit = "ergs/s/cm^2/Angstrom", errors = errors, source = source)

if docfaibcspectra:
    for name in sorted(next(os.walk("../sne-external-spectra/CfA_SNIbc"))[1], key=lambda s: s.lower()):
        fullpath = "../sne-external-spectra/CfA_SNIbc/" + name
        if name[:2] == 'sn' and is_number(name[2:6]):
            name = 'SN' + name[2:]
        name = add_event(name)
        reference = "<a href='https://www.cfa.harvard.edu/supernova/SNarchive.html'>CfA Supernova Archive</a>"
        source = get_source(name, reference, secondary = True)
        for file in sorted(glob.glob(fullpath + '/*'), key=lambda s: s.lower()):
            if os.path.basename(file) == 'sn1993J-19941128.flm':
                print ('Warning: Need to fix sn1993J-19941128.flm!')
                continue
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
            add_spectrum(name = name, waveunit = 'Angstrom', fluxunit = 'erg/s/cm^2/Angstrom', wavelengths = wavelengths,
                fluxes = fluxes, timeunit = 'MJD', time = time, instrument = instrument, source = source)

if writeevents:
    # Calculate some columns based on imported data, sanitize some fields
    for name in events:
        if 'claimedtype' in events[name]:
            events[name]['claimedtype'][:] = [ct for ct in events[name]['claimedtype'] if (ct != '?' and ct != '-')]
        if 'redshift' in events[name] and 'hvel' not in events[name]:
            z = float(events[name]['redshift'])
            sigdigits = len(events[name]['redshift'].strip('0').strip('.'))
            events[name]['hvel'] = pretty_num(clight/1.e5*((z + 1.)**2. - 1.)/
                ((z + 1.)**2. + 1.), sig = sigdigits)
        elif 'hvel' in events[name] and 'z' not in events[name]:
            voc = float(events[name]['hvel'])*1.e5/clight
            sigdigits = len(events[name]['hvel'].strip('0').strip('.'))
            events[name]['redshift'] = pretty_num(sqrt((1. + voc)/(1. - voc)) - 1., sig = sigdigits)
        if 'redshift' in events[name] and float(events[name]['redshift']) > 0.0:
            dl = cosmo.luminosity_distance(float(events[name]['redshift']))
            sigdigits = len(events[name]['redshift'].strip('0').strip('.'))
            if 'lumdist' not in events[name]:
                events[name]['lumdist'] = pretty_num(dl.value, sig = sigdigits)
            if 'maxabsmag' not in events[name] and 'maxappmag' in events[name]:
                events[name]['maxabsmag'] = pretty_num(float(events[name]['maxappmag']) - 5.0*(log10(dl.to('pc').value) - 1.0), sig = sigdigits)
        if 'redshift' in events[name]:
            events[name]['redshift'] = pretty_num(Decimal(events[name]['redshift']))
        if 'hvel' in events[name]:
            events[name]['hvel'] = pretty_num(Decimal(events[name]['hvel']))
        if 'host' in events[name]:
            events[name]['host'] = events[name]['host'].replace("NGC", "NGC ")
            events[name]['host'] = events[name]['host'].replace("UGC", "UGC ")
            events[name]['host'] = events[name]['host'].replace("IC", "IC ")
            events[name]['host'] = ' '.join(events[name]['host'].split())
        if 'photometry' in events[name]:
            events[name]['photometry'].sort(key=lambda x: float(x['time']))
        if 'spectra' in events[name]:
            events[name]['spectra'].sort(key=lambda x: float(x['time']))
        events[name] = OrderedDict(sorted(events[name].items(), key=lambda key: event_attr_priority(key[0])))

    # Write it all out!
    for name in events:
        print('Writing ' + name)
        filename = event_filename(name)

        jsonstring = json.dumps({name:events[name]}, indent=4, separators=(',', ': '), ensure_ascii=False)

        outdir = '../'
        if 'discoveryear' in events[name]:
            for r, year in enumerate(repyears):
                if int(events[name]['discoveryear']) <= year:
                    outdir += repfolders[r]
                    break
        else:
            outdir += str(repfolders[0])

        f = codecs.open(outdir + '/' + filename + '.json', 'w', encoding='utf8')
        f.write(jsonstring)
        f.close()

# Print some useful facts
if printextra:
    print('Printing events without any photometry:')
    for name in events:
        if 'photometry' not in events[name]:
            print(name)
