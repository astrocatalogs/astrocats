#!/usr/local/bin/python2.7

import csv
import glob
import os
import re
import urllib2
from collections import OrderedDict
from pprint import pprint
from math import log10, floor, sqrt
from BeautifulSoup import BeautifulSoup, SoupStrainer
from operator import itemgetter

outdir = '../'

eventnames = []

docfa = 		True
dosuspect = 	True
doucb = 		True
dosdss = 		True
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
	"discoverer"
	]

columnkey.sort(key=str.lower)

events = {}
eventphotometry = {}

def newevent(name):
	print name
	events[name] = OrderedDict.fromkeys(columnkey, '')
	eventphotometry[name] = []

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
							photometryrow = ['photometry', 'MJD', mjd, band, '', mag, err, 0]
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
							eventphotometry[name].append(['photometry', tuout, mjd, eventbands[(v-1)/2], '', row[v], row[v+1], 0])

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

		for r, row in enumerate(tsvin):
			if len(row) > 0 and row[0] == "#":
				continue
			mjd = row[0]
			abmag = row[1]
			aberr = row[2]
			band = row[4]
			instrument = row[5]
			eventphotometry[name].append(['photometry', 'MJD', mjd, band, instrument, abmag, aberr, 0])
	
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
				events[name]['snra'] = row[-4]
				events[name]['sndec'] = row[-2]
			if r == 1:
				events[name]['redshift'] = row[1]
			if r >= 19:
				mjd = row[1]
				band = sdssbands[int(row[2])]
				abmag = row[3]
				aberr = row[4]
				instrument = "SDSS"
				eventphotometry[name].append(['photometry', 'MJD', mjd, band, instrument, abmag, aberr, 0])


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
			hostname = record[2]
			galra = record[3]
			galdec = record[4]
			snra = record[5]
			sndec = record[6]
			redvel = record[11].strip(':')
			discoverer = record[19]

			hvel = ''
			redshift = ''
			if redvel != '':
				c = 29979245800.
				if round(float(redvel)) == float(redvel):
					hvel = int(redvel)
					voc = float(hvel)*1.e5/c
					redshift = round_sig(sqrt((1. + voc)/(1. - voc)) - 1., sig = 3)
				else:
					redshift = float(redvel)
					hvel = round(round_sig(c/1.e5*((redshift + 1.)**2. - 1.)/((redshift + 1.)**2. + 1.), sig = 3))
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
			if not events[name][key].isspace():
				csvout.writerow([key, events[name][key]])
		
		sortedphotometry = sorted(eventphotometry[name], key=itemgetter(2))
		for row in sortedphotometry:
			csvout.writerow(row)

		outfile.close()

