#!/usr/local/bin/python2.7

import csv
import glob
import os
import re
import urllib
import urllib2
import calendar
from collections import OrderedDict
from math import log10, floor, sqrt
from BeautifulSoup import BeautifulSoup, SoupStrainer
from operator import itemgetter
from subprocess import Popen

clight = 29979245800.

outdir = '../data/'

eventnames = []

doitep =		False

docfa = 		True
dosuspect = 	True
doucb = 		True
dosdss = 		True
dogaia =		True
docsp =			True
doasiago = 		True

writeevents = 	True

columnkey = [
	"host",
	"citations",
	"instruments",
	"redshift",
	"hvel",
	"claimedtype",
	"notes",
	"galra",
	"galdec",
	"snra",
	"sndec",
	"discoverer",
	"year",
	"discovermonth",
	"discoverday",
	"maxmonth",
	"maxday"
	]

columnkey.sort(key=str.lower)

events = {}
eventphotometry = {}
eventsources = {}

def newevent(name):
	print name
	events[name] = OrderedDict.fromkeys(columnkey, '')
	eventphotometry[name] = []
	eventsources[name] = []

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
				refc = 0
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
							if name not in events:
								newevent(name)
							year = re.findall(r'\d+', name)[0]
							events[name]['year'] = year
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

						reference = ''
						for link in bandsoup.body.findAll('a'):
							if 'adsabs' in link['href']:
								reference = str(link).replace('"', "'")

						if reference:
							if len(eventsources[name]) == 0 or reference not in [eventsources[name][es][2] for es in xrange(len(eventsources[name]))]:
								refc = refc + 1
								eventsources[name].append(['source', 'name', reference, 'alias', refc])
								alias = refc
							else:
								alias = [eventsources[name][es][4] for es in xrange(len(eventsources[name]))][
									[eventsources[name][es][2] for es in xrange(len(eventsources[name]))].index(reference)]
						else:
							refc = refc + 1
							eventsources[name].append(['source', 'name', 'SUSPECT', 'alias', refc])
							alias = refc

						for r, row in enumerate(bandtable.findAll('tr')):
							if r == 0:
								continue
							col = row.findAll('td')
							mjd = str(float(col[0].renderContents()) - 2400000.5)
							mag = col[3].renderContents()
							if mag.isspace():
								mag = ''
							err = col[4].renderContents()
							if err.isspace():
								err = ''
							photometryrow = ['photometry', 'timeunit', 'MJD', 'time', mjd, 'band', band, 'abmag', mag] + (['aberr', err] if err else []) + ['source', alias]
							eventphotometry[name].append(photometryrow)


# CfA data
if docfa:
	for file in sorted(glob.glob("../external/cfa-input/*.dat"), key=lambda s: s.lower()):
		tsvin = open(file,'rb')
		tsvin = csv.reader(tsvin, delimiter=' ', skipinitialspace=True)
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
		if name not in events:
			newevent(name)

		year = re.findall(r'\d+', name)[0]
		events[name]['year'] = year

		eventbands = list(eventparts[1])

		tu = 'MJD'
		jdoffset = 0.
		for rc, row in enumerate(csv_data):
			if len(row) > 0 and row[0][0] == "#":
				if len(row[0]) > 2 and row[0][:3] == "#JD":
					tu = 'JD'
					rowparts = row[0].split('-')
					jdoffset = float(rowparts[1])
				elif len(row) > 1 and row[1] == "HJD":
					tu = "HJD"
				continue
			elif len(row) > 0:
				mjd = row[0]
				for v, val in enumerate(row):
					if v == 0:
						if tu == 'JD':
							mjd = str(float(val) + jdoffset - 2400000.5)
							tuout = 'MJD'
						elif tu == 'HJD':
							mjd = str(float(val) - 2400000.5)
							tuout = 'MJD'
						else:
							mjd = val
							tuout = tu
					elif v % 2 != 0:
						if float(row[v]) < 90.0:
							eventphotometry[name].append(['photometry', 'timeunit', tuout, 'time', mjd, 'band', eventbands[(v-1)/2], 'abmag', row[v], 'aberr', row[v+1], 'upperlimit', 0])

# Now import the UCB SNDB
if doucb:
	for file in sorted(glob.glob("../external/SNDB/*.dat"), key=lambda s: s.lower()):
		tsvin = open(file,'rb')
		tsvin = csv.reader(tsvin, delimiter=' ', skipinitialspace=True)

		eventname = os.path.basename(os.path.splitext(file)[0])

		eventparts = eventname.split('.')

		name = snname(eventparts[0])
		if name not in events:
			newevent(name)

		year = re.findall(r'\d+', name)[0]
		events[name]['year'] = year

		for r, row in enumerate(tsvin):
			if len(row) > 0 and row[0] == "#":
				continue
			mjd = row[0]
			abmag = row[1]
			aberr = row[2]
			band = row[4]
			instrument = row[5]
			eventphotometry[name].append(['photometry', 'timeunit', 'MJD', 'time', mjd, 'band', band, 'instrument', instrument, 'abmag', abmag, 'aberr', aberr, 'upperlimit', 0])
	
# Import SDSS
sdssbands = ['u', 'g', 'r', 'i', 'z']

if dosdss:
	for file in sorted(glob.glob("../external/SDSS/*.sum"), key=lambda s: s.lower()):
		tsvin = open(file,'rb')
		tsvin = csv.reader(tsvin, delimiter=' ', skipinitialspace=True)

		for r, row in enumerate(tsvin):
			if r == 0:
				if row[5] == "RA:":
					name = "SDSS" + row[3]
				else:
					name = "SN" + row[5]
				if name not in events:
					newevent(name)

				if row[5] != "RA:":
					year = re.findall(r'\d+', name)[0]
					events[name]['year'] = year

				events[name]['snra'] = row[-4]
				events[name]['sndec'] = row[-2]
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
				eventphotometry[name].append(['photometry', 'timeunit', 'MJD', 'time', mjd, 'band', band, 'instrument', instrument, 'abmag', abmag, 'aberr', aberr, 'upperlimit', 0])

#Import GAIA
if dogaia:
	response = urllib2.urlopen('https://gaia.ac.uk/selected-gaia-science-alerts')
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

		if name not in events:
			newevent(name)

		year = '20' + re.findall(r'\d+', name)[0]
		events[name]['year'] = year

		events[name]['snra'] = col[2].renderContents().strip()
		events[name]['sndec'] = col[3].renderContents().strip()
		events[name]['claimedtype'] = classname.replace('SN', '').strip()

		photlink = 'http://gsaweb.ast.cam.ac.uk/alerts/alert/' + name + '/lightcurve.csv/'
		photresp = urllib2.urlopen(photlink)
		photsoup = BeautifulSoup(photresp)
		photodata = photsoup.renderContents().split('\n')[2:-1]
		for ph in photodata:
			photo = ph.split(',')
			mjd = str(float(photo[1].strip()) - 2400000.5)
			abmag = photo[2].strip()
			aberr = 0.
			instrument = 'GAIA'
			band = 'G'
			eventphotometry[name].append(['photometry', 'timeunit', 'MJD', 'time', mjd, 'band', band, 'instrument', instrument, 'abmag', abmag, 'aberr', aberr, 'upperlimit', 0])

# Import ITEP, not currently working.
if doitep:
	url = 'http://dau.itep.ru/sn/lctheory/draw'
	values = {'edit-snnameobs' : '1993J' }

	data = urllib.urlencode(values)
	req = urllib2.Request(url, data)
	response = urllib2.urlopen(req)
	the_page = response.read()

# Import CSP
cspbands = ['u', 'B', 'V', 'g', 'r', 'i', 'Y', 'J', 'H', 'K']

if docsp:
	for file in sorted(glob.glob("../external/CSP/*.dat"), key=lambda s: s.lower()):
		tsvin = open(file,'rb')
		tsvin = csv.reader(tsvin, delimiter='\t', skipinitialspace=True)

		eventname = os.path.basename(os.path.splitext(file)[0])

		eventparts = eventname.split('opt+')

		name = snname(eventparts[0])
		if name not in events:
			newevent(name)

		year = re.findall(r'\d+', name)[0]
		events[name]['year'] = year

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
						eventphotometry[name].append(['photometry', 'timeunit', 'MJD', 'time', mjd, 'band', cspbands[(v-1)/2], 'instrument', 'CSP', 'abmag', row[v], 'aberr', row[v+1], 'upperlimit', 0])

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
			if name not in events:
				newevent(name)

			year = re.findall(r'\d+', name)[0]
			events[name]['year'] = year

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

# Write it all out!
if writeevents:
	for name in events:
		print 'writing ' + name
		outfile = open(outdir + name + '.dat', 'wb')
		csvout = csv.writer(outfile, quotechar='"', quoting=csv.QUOTE_ALL, delimiter="\t")
		
		csvout.writerow(['name', name])

		for key in events[name]:
			if events[name][key]:
				csvout.writerow([key, events[name][key]])

		for row in eventsources[name]:
			csvout.writerow(row)
		
		sortedphotometry = sorted(eventphotometry[name], key=itemgetter(2))
		for row in sortedphotometry:
			csvout.writerow(row)

		outfile.close()

		# Compress the output
		Popen(["bzip2", '-f', outdir + name + '.dat'])

