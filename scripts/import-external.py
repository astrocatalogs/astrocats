#!/usr/local/bin/python2.7

import csv
import glob
import os
import re
import urllib2
import calendar
import sys
import subprocess
from astroquery.vizier import Vizier
from astropy.time import Time as astrotime
from collections import OrderedDict
from sortedcontainers import SortedDict
from math import log10, floor, sqrt
from BeautifulSoup import BeautifulSoup, SoupStrainer
from operator import itemgetter

clight = 29979245800.

outdir = '../data/'

eventnames = []

dovizier =       True
dosuspect =      True
docfa =          True
doucb =          True
dosdss =         True
dogaia =         True
docsp =          True
doitep =         True
doasiago =       True
dorochester =    True
dofirstmax =     True
dolennarz =      True
writeevents =    True
compressevents = True
printextra =     False

photometrykeys = [
    'timeunit',
    'time',
    'band',
    'instrument',
    'abmag',
    'aberr',
    'source'
]

events = OrderedDict()
eventphotometry = OrderedDict()
eventsources = OrderedDict()

def add_event(name):
    if name not in events:
        for event in events:
            if len(events[event]['aliases']) > 1 and name in events[event]['aliases']:
                return event
        print name
        events[name] = SortedDict()
        eventphotometry[name] = []
        eventsources[name] = []
        add_alias(name, name)
        return name
    else:
        return name

def event_filename(name):
    return(name.replace('/', '_'))

def add_alias(name, alias):
    if 'aliases' in events[name]:
        if alias not in events[name]['aliases']:
            events[name]['aliases'].append(alias)
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

def round_sig(x, sig=2):
    if x == 0.0:
        return 0.0
    return round(x, sig-int(floor(log10(abs(x))))-1)

def get_source(name, reference, secondary = ''):
    if len(eventsources[name]) == 0 or reference not in [eventsources[name][es][2] for es in xrange(len(eventsources[name]))]:
        source = str(len(eventsources[name]) + 1)
        eventsources[name].append(['source', 'name', reference, 'alias', source] + (['secondary', 1] if secondary else []))
    else:
        source = [eventsources[name][es][4] for es in xrange(len(eventsources[name]))][
            [eventsources[name][es][2] for es in xrange(len(eventsources[name]))].index(reference)]
    return source

def add_photometry(name, timeunit = "MJD", time = "", instrument = "", band = "", abmag = "", aberr = "", source = ""):
    if not time or not abmag:
        print 'Error: Time or AB mag not specified when adding photometry.\n'
        print 'Name : "' + name + '", Time: "' + time + '", Band: "' + band + '", AB mag: "' + abmag + '"'
        sys.exit()

    # Look for duplicate data and don't add if duplicate
    for photo in eventphotometry[name]:
        if (photo['timeunit'] == timeunit and float(photo['time']) == float(time) and
            photo['band'] == band and float(photo['abmag']) == float(abmag) and
            ((not photo['aberr'] and not aberr) or (photo['aberr'] and aberr and float(photo['aberr']) == float(aberr)) or
            (photo['aberr'] and not aberr))):
            return

    photoentry = OrderedDict.fromkeys(photometrykeys)
    photoentry['timeunit'] = timeunit
    photoentry['time'] = time
    photoentry['band'] = band
    photoentry['abmag'] = abmag
    if instrument:
        photoentry['instrument'] = instrument
    if aberr:
        photoentry['aberr'] = aberr
    if source:
        photoentry['source'] = source
    eventphotometry[name].append(photoentry)

def get_max_light(name):
    if not eventphotometry[name]:
        return None

    eventphoto = [eventphotometry[name][x]['abmag'] for x in xrange(len(eventphotometry[name]))]
    mlindex = eventphoto.index(min(eventphoto))
    mlmjd = float(eventphotometry[name][mlindex]['time'])
    return astrotime(mlmjd, format='mjd').datetime

def get_first_light(name):
    if not eventphotometry[name]:
        return None

    eventtime = [eventphotometry[name][x]['time'] for x in xrange(len(eventphotometry[name]))]
    flindex = eventtime.index(min(eventtime))
    flmjd = float(eventphotometry[name][flindex]['time'])
    return astrotime(flmjd, format='mjd').datetime

def set_first_max_light(name):
    mldt = get_max_light(name)
    if mldt:
        events[name]['maxyear'] = mldt.year
        events[name]['maxmonth'] = mldt.month
        events[name]['maxday'] = mldt.day

    fldt = get_first_light(name)
    if fldt:
        events[name]['discoveryear'] = fldt.year
        events[name]['discovermonth'] = fldt.month
        events[name]['discoverday'] = fldt.day

def jd_to_mjd(jd):
    return jd - 2400000.5

# Import primary data sources from Vizier
if dovizier:
    Vizier.ROW_LIMIT = -1
    result = Vizier.get_catalogs("VII/272/snrs")
    table = result[result.keys()[0]]

    for row in table:
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
    table = result[result.keys()[0]]
    for row in table:
        name = 'SN' + row['SN']
        name = add_event(name)
        events[name]['redshift'] = row['zhost']
        events[name]['ebv'] = row['E_B-V_']
    result = Vizier.get_catalogs("J/MNRAS/442/844/table2")
    table = result[result.keys()[0]]
    for row in table:
        name = 'SN' + row['SN']
        source = get_source(name, reference)
        if row['Bmag']:
            add_photometry(name, time = row['MJD'], band = 'B', abmag = row['Bmag'], aberr = row['e_Bmag'], source = source)
        if row['Vmag']:
            add_photometry(name, time = row['MJD'], band = 'V', abmag = row['Vmag'], aberr = row['e_Vmag'], source = source)
        if row['Rmag']:
            add_photometry(name, time = row['MJD'], band = 'R', abmag = row['Rmag'], aberr = row['e_Rmag'], source = source)
        if row['Imag']:
            add_photometry(name, time = row['MJD'], band = 'I', abmag = row['Imag'], aberr = row['e_Imag'], source = source)

    reference = "<a href='http://adsabs.harvard.edu/abs/2015ApJS..219...13W'>2015ApJS..219...13W</a>"
    result = Vizier.get_catalogs("J/ApJS/219/13/table3")
    table = result[result.keys()[0]]
    for row in table:
        name = 'LSQ' + row['LSQ']
        name = add_event(name)
        events[name]['snra'] = row['RAJ2000']
        events[name]['sndec'] = row['DEJ2000']
        events[name]['redshift'] = row['z']
        events[name]['ebv'] = row['E_B-V_']
    result = Vizier.get_catalogs("J/ApJS/219/13/table2")
    table = result[result.keys()[0]]
    for row in table:
        name = 'LSQ' + row['LSQ']
        source = get_source(name, reference)
        add_photometry(name, time = jd_to_mjd(row['JD']), instrument = 'La Silla-QUEST', band = row['Filt'], abmag = row['mag'], aberr = row['e_mag'], source = source)

# Suspect catalog
if dosuspect:
    response = urllib2.urlopen('http://www.nhn.ou.edu/cgi-bin/cgiwrap/~suspect/snindex.cgi')

    soup = BeautifulSoup(response)
    i = 0
    for a in soup.findAll('a'):
        if 'phot=yes' in a['href'] and not 'spec=yes' in a['href']:
            if int(a.contents[0]) > 0:
                i = i + 1
                photlink = 'http://www.nhn.ou.edu/cgi-bin/cgiwrap/~suspect/' + a['href']
                eventresp = urllib2.urlopen(photlink)
                eventsoup = BeautifulSoup(eventresp)
                ei = 0
                for ea in eventsoup.findAll('a'):
                    if ea.contents[0] == 'I':
                        ei = ei + 1
                        bandlink = 'http://www.nhn.ou.edu/cgi-bin/cgiwrap/~suspect/' + ea['href']
                        bandresp = urllib2.urlopen(bandlink)
                        bandsoup = BeautifulSoup(bandresp)
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
                            events[name]['claimedtype'] = types[0].split(':')[1].strip().split(' ')[0]

                        bands = bandsoup.body.findAll(text=re.compile("^Band"))
                        band = bands[0].split(':')[1].strip()

                        secondaryreference = "<a href='https://www.nhn.ou.edu/~suspect/'>SUSPECT</a>"
                        secondarysource = get_source(name, secondaryreference, 1)

                        reference = ''
                        for link in bandsoup.body.findAll('a'):
                            if 'adsabs' in link['href']:
                                reference = str(link).replace('"', "'")
                        source = get_source(name, reference)

                        for r, row in enumerate(bandtable.findAll('tr')):
                            if r == 0:
                                continue
                            col = row.findAll('td')
                            mjd = str(jd_to_mjd(float(col[0].renderContents())))
                            mag = col[3].renderContents()
                            if mag.isspace():
                                mag = ''
                            else:
                                mag = str(float(mag))
                            aberr = col[4].renderContents()
                            if aberr.isspace():
                                aberr = ''
                            else:
                                aberr = str(float(aberr))
                            add_photometry(name, time = mjd, band = band, abmag = mag, aberr = aberr, source = secondarysource + ',' + source)

# CfA data
if docfa:
    for file in sorted(glob.glob("../external/cfa-input/*.dat"), key=lambda s: s.lower()):
        f = open(file,'rb')
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
            csv_data[r] = filter(None, csv_data[r])

        eventname = os.path.basename(os.path.splitext(file)[0])

        eventparts = eventname.split('_')

        name = snname(eventparts[0])
        name = add_event(name)

        year = re.findall(r'\d+', name)[0]
        events[name]['discoveryear'] = year

        eventbands = list(eventparts[1])

        tu = 'MJD'
        jdoffset = 0.
        for rc, row in enumerate(csv_data):
            if len(row) > 0 and row[0][0] == "#":
                if len(row[0]) > 2 and row[0][:3] == "#JD":
                    tu = 'JD'
                    rowparts = row[0].split('-')
                    jdoffset = float(rowparts[1])
                elif len(row[0]) > 6 and row[0][:7] == "#Julian":
                    tu = 'JD'
                    jdoffset = 0.
                elif len(row) > 1 and row[1].lower() == "photometry":
                    for ci, col in enumerate(row[2:]):
                        if col[0] == "(":
                            refstr = ' '.join(row[2+ci:])
                            refstr = refstr.translate(None, '()')
                            bibcode = refstr
                            print bibcode
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
                            mjd = str(jd_to_mjd(float(val) + jdoffset))
                            tuout = 'MJD'
                        elif tu == 'HJD':
                            mjd = str(jd_to_mjd(float(val)))
                            tuout = 'MJD'
                        else:
                            mjd = val
                            tuout = tu
                    elif v % 2 != 0:
                        if float(row[v]) < 90.0:
                            add_photometry(name, timeunit = tuout, time = mjd, band = eventbands[(v-1)/2], abmag = row[v], aberr = row[v+1], source = secondarysource + ',' + source)
        f.close()

    # Hicken 2012
    f = open("../external/hicken-2012-standard.dat", 'rb')
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
    tsvin = open("../external/bianco-2014-standard.dat", 'rb')
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
        f = open(file,'rb')
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
        f = open(file,'rb')
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
    response = urllib2.urlopen('file:///var/www/html/sne/sne/external/selected-gaia-science-alerts')
    html = response.read()

    soup = BeautifulSoup(html)
    table = soup.findAll("table")[1]
    for r, row in enumerate(table.findAll('tr')):
        if r == 0:
            continue

        col = row.findAll('td')
        classname = col[7].renderContents()

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

        events[name]['snra'] = col[2].renderContents().strip()
        events[name]['sndec'] = col[3].renderContents().strip()
        events[name]['claimedtype'] = classname.replace('SN', '').strip()

        photlink = 'http://gsaweb.ast.cam.ac.uk/alerts/alert/' + name + '/lightcurve.csv/'
        photresp = urllib2.urlopen(photlink)
        photsoup = BeautifulSoup(photresp)
        photodata = photsoup.renderContents().split('\n')[2:-1]
        for ph in photodata:
            photo = ph.split(',')
            mjd = str(jd_to_mjd(float(photo[1].strip())))
            abmag = photo[2].strip()
            aberr = 0.
            instrument = 'GAIA'
            band = 'G'
            add_photometry(name, time = mjd, instrument = instrument, band = band, abmag = abmag, aberr = aberr, source = source)

# Import CSP
if docsp:
    cspbands = ['u', 'B', 'V', 'g', 'r', 'i', 'Y', 'J', 'H', 'K']
    for file in sorted(glob.glob("../external/CSP/*.dat"), key=lambda s: s.lower()):
        f = open(file,'rb')
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
                    redshift = float(events[name]['redshift'])
                    events[name]['hvel'] = round(round_sig(clight/1.e5*((redshift + 1.)**2. - 1.)/
                        ((redshift + 1.)**2. + 1.), sig = 3))
                    events[name]['snra'] = row[1].split(' ')[-1]
                    events[name]['sndec'] = row[2].split(' ')[-1]
                continue
            for v, val in enumerate(row):
                if v == 0:
                    mjd = val
                elif v % 2 != 0:
                    if float(row[v]) < 90.0:
                        add_photometry(name, time = mjd, instrument = 'CSP', band = cspbands[(v-1)/2], abmag = row[v], aberr = row[v+1], source = source)
        f.close()

# Import ITEP
if doitep:
    f = open("../external/itep-lc-cat-28dec2015.txt",'rb')
    tsvin = csv.reader(f, delimiter='|', skipinitialspace=True)
    curname = ''
    for r, row in enumerate(tsvin):
        if r <= 1 or len(row) < 7:
            continue
        name = 'SN' + row[0].strip()
        mjd = str(jd_to_mjd(float(row[1].strip())))
        band = row[2].strip()
        abmag = row[3].strip()
        aberr = row[4].strip()
        reference = row[6].strip().strip(',')
        if curname != name:
            curname = name
            name = add_event(name)
            year = re.findall(r'\d+', name)[0]
            events[name]['discoveryear'] = year

            secondaryreference = "<a href='http://dau.itep.ru/sn/node/72'>Sternberg Astronomical Institute Supernova Light Curve Catalogue</a>"
            secondarysource = get_source(name, secondaryreference, 1)

        source = get_source(name, reference) if reference else ''

        add_photometry(name, time = mjd, band = band, abmag = abmag, aberr = aberr, source = secondarysource + ',' + source)
    f.close


# Now import the Asiago catalog
if doasiago:
    response = urllib2.urlopen('http://graspa.oapd.inaf.it/cgi-bin/sncat.php')
    html = response.read()
    html = html.replace('\r', '')

    soup = BeautifulSoup(html)
    table = soup.find("table")

    records = []
    for r, row in enumerate(table.findAll('tr')):
        if r == 0:
            continue

        col = row.findAll('td')
        records.append([x.renderContents() for x in col])

    for record in records:
        if len(record) > 1 and record[1] != '':
            name = snname("SN" + record[1])
            name = add_event(name)

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
                    voc = float(hvel)*1.e5/clight
                    redshift = round_sig(sqrt((1. + voc)/(1. - voc)) - 1., sig = 3)
                else:
                    redshift = float(redvel)
                    hvel = round(round_sig(clight/1.e5*((redshift + 1.)**2. - 1.)/((redshift + 1.)**2. + 1.), sig = 3))
                redshift = str(redshift)
                hvel = str(hvel)

            claimedtype = record[17].strip(':')

            if (hostname != ''):
                events[name]['host'] = hostname
            if (claimedtype != ''):
                events[name]['claimedtype'] = claimedtype
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
    #response = urllib2.urlopen('file:///var/www/html/sne/sne/external/snredshiftall.html')
    response = urllib2.urlopen('http://www.rochesterastronomy.org/snimages/snredshiftall.html')
    html = response.read()

    soup = BeautifulSoup(html)
    rows = soup.findAll('tr')
    secondaryreference = "<a href='http://www.rochesterastronomy.org/snimages/snredshiftall.html'>Latest Supernovae</a>"
    for r, row in enumerate(rows):
        if r == 0:
            continue
        cols = row.findAll('td')
        name = re.sub('<[^<]+?>', '', str(cols[0].contents[0])).strip()
        if name[:4].isdigit():
            name = 'SN' + name
        name = add_event(name)

        secondarysource = get_source(name, secondaryreference, secondary = 1)

        if str(cols[1].contents[0]).strip() != 'unk':
            events[name]['claimedtype'] = str(cols[1].contents[0]).strip()
        if str(cols[2].contents[0]).strip() != 'anonymous':
            events[name]['host'] = str(cols[2].contents[0]).strip()
        events[name]['snra'] = str(cols[3].contents[0]).strip()
        events[name]['sndec'] = str(cols[4].contents[0]).strip()
        if str(cols[6].contents[0]).strip() not in ['2440587', '2440587.292']:
            astrot = astrotime(float(str(cols[6].contents[0]).strip()), format='jd')
            events[name]['discoverday'] = astrot.datetime.day
            events[name]['discovermonth'] = astrot.datetime.month
            events[name]['discoveryear'] = astrot.datetime.year
        if str(cols[7].contents[0]).strip() not in ['2440587', '2440587.292']:
            astrot = astrotime(float(str(cols[7].contents[0]).strip()), format='jd')
            events[name]['maxday'] = astrot.datetime.day
            events[name]['maxmonth'] = astrot.datetime.month
            events[name]['maxyear'] = astrot.datetime.year
            source = get_source(name, str(cols[12].contents[0]).strip().replace('"', "'"))
            if float(str(cols[8].contents[0]).strip()) <= 90.0:
                add_photometry(name, time = astrot.mjd, abmag = float(str(cols[8].contents[0]).strip()), source = ','.join([source, secondarysource]))
        if cols[11].contents[0] != 'n/a':
            events[name]['redshift'] = float(str(cols[11].contents[0]).strip())
        events[name]['discoverer'] = str(cols[13].contents[0]).strip()
        if cols[14].contents:
            add_alias(name, str(cols[14].contents[0]).strip())

if dofirstmax:
    for name in events:
        set_first_max_light(name)

if dolennarz:
    Vizier.ROW_LIMIT = -1
    result = Vizier.get_catalogs("J/A+A/538/A120/usc")
    table = result[result.keys()[0]]

    reference = "<a href='http://adsabs.harvard.edu/abs/2012A%26A...538A.120L'>2012A&A...538A.120L</a>"
    for row in table:
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

            if not eventphotometry[name]:
                if row['Dmag']:
                    mjd = astrot.mjd
                    add_photometry(name, time = mjd, band = row['Dband'], abmag = row['Dmag'], source = source)
            if 'discoveryear' not in events[name] and 'discovermonth' not in events[name] and 'discoverday' not in events[name]:
                events[name]['discoveryear'] = astrot.datetime.year
                if len(dateparts) >= 2:
                    events[name]['discovermonth'] = astrot.datetime.month
                if len(dateparts) == 3:
                    events[name]['discoverday'] = astrot.datetime.day
        if row['Mdate']:
            dateparts = row['Mdate'].split('-')
            if len(dateparts) == 3:
                astrot = astrotime(row['Mdate'], scale='utc')
            elif len(dateparts) == 2:
                astrot = astrotime(row['Mdate'] + '-01', scale='utc')
            else:
                astrot = astrotime(row['Mdate'] + '-01-01', scale='utc')

            if not eventphotometry[name]:
                if row['Mmag']:
                    mjd = astrot.mjd
                    add_photometry(name, time = mjd, band = row['Mband'], abmag = row['Mmag'], source = source)
            if 'maxyear' not in events[name] and 'maxmonth' not in events[name] and 'maxday' not in events[name]:
                events[name]['maxyear'] = astrot.datetime.year
                if len(dateparts) >= 2:
                    events[name]['maxmonth'] = astrot.datetime.month
                if len(dateparts) == 3:
                    events[name]['maxday'] = astrot.datetime.day

if writeevents:
    # Calculate some columns based on imported data, sanitize some fields
    for name in events:
        if 'claimedtype' in events[name] and events[name]['claimedtype'] == '?':
            del events[name]['claimedtype']
        if 'hvel' in events[name]:
            events[name]['hvel'] = '%g'%(float(events[name]['hvel']))
        if 'host' in events[name]:
            events[name]['host'] = events[name]['host'].replace("NGC", "NGC ")
            events[name]['host'] = events[name]['host'].replace("UGC", "UGC ")
            events[name]['host'] = events[name]['host'].replace("IC", "IC ")
            events[name]['host'] = ' '.join(events[name]['host'].split())

    # Write it all out!
    for name in events:
        print 'Writing ' + name
        filename = event_filename(name)
        outfile = open(outdir + filename + '.dat', 'wb')
        csvout = csv.writer(outfile, quotechar='"', quoting=csv.QUOTE_ALL, delimiter="\t")
        
        csvout.writerow(['name', name])

        for row in eventsources[name]:
            csvout.writerow(row)
        
        for key in events[name]:
            if events[name][key]:
                if isinstance(events[name][key], list) and not isinstance(events[name][key], basestring):
                    events[name][key] = ",".join(events[name][key])
                csvout.writerow([key, events[name][key]])

        sortedphotometry = sorted(eventphotometry[name], key=itemgetter('time'))
        for entry in sortedphotometry:
            row = ['photometry'] + list(sum(((k, v) for k, v in entry.items() if v), ()))
            csvout.writerow(row)

        outfile.close()

# Compress the output
if compressevents:
    print 'Compressing output...'
    print 'bzip2 -k -f ' + outdir + '*.dat'
    os.system('bzip2 -k -f ' + outdir + '*.dat')

# Print some useful facts
if printextra:
    print 'Printing events without any photometry:'
    for name in events:
        if not eventphotometry[name]:
            print name
