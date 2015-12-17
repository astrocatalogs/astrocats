#!/usr/local/bin/python2.7

import csv
import glob
import os
import re
import urllib2
from BeautifulSoup import BeautifulSoup

outdir = '../'

eventnames = []

def snname(string):
	newstring = string.replace(' ', '').upper()
	if (newstring[:2] == "SN"):
		head = newstring[:6]
		tail = newstring[6:]
		if len(tail) >= 2 and tail[1] != '?':
			tail = tail.lower()
		newstring = head + tail

	return newstring

# First import the CfA data.
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
	eventname = snname(eventparts[0])
	eventbands = list(eventparts[1])

	if (eventname in eventnames):
		outfile = open(outdir + eventname + '.dat', 'a')
		csvout = csv.writer(outfile, quotechar='"', quoting=csv.QUOTE_ALL, delimiter="\t")
	else:
		eventnames.append(eventname)
		outfile = open(outdir + eventname + '.dat', 'wb')
		csvout = csv.writer(outfile, quotechar='"', quoting=csv.QUOTE_ALL, delimiter="\t")

		csvout.writerow(['name', eventname])

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
						mjd = float(val) + jdoffset - 2400000.5
						tuout = 'MJD'
					elif tu == 'HJD':
						mjd = float(val) - 2400000.5
						tuout = 'MJD'
					else:
						mjd = val
						tuout = tu
				elif v % 2 != 0:
					if float(row[v]) < 90.0:
						csvout.writerow(['photometry', tuout, mjd, eventbands[(v-1)/2], '', row[v], row[v+1], 0])

	outfile.close()

# Now import the UCB SNDB
for file in sorted(glob.glob("../external/SNDB/*.dat"), key=lambda s: s.lower()):
	tsvin = open(file,'rb')
	tsvin = csv.reader(tsvin, delimiter=' ', skipinitialspace=True)

	eventname = os.path.basename(os.path.splitext(file)[0])

	eventparts = eventname.split('.')
	eventname = snname(eventparts[0])

	if (eventname in eventnames):
		outfile = open(outdir + eventname + '.dat', 'a')
		csvout = csv.writer(outfile, quotechar='"', quoting=csv.QUOTE_ALL, delimiter="\t")
	else:
		eventnames.append(eventname)
		outfile = open(outdir + eventname + '.dat', 'wb')
		csvout = csv.writer(outfile, quotechar='"', quoting=csv.QUOTE_ALL, delimiter="\t")

		csvout.writerow(['name', eventname])

	for r, row in enumerate(tsvin):
		if len(row) > 0 and row[0] == "#":
			continue
		mjd = row[0]
		abmag = row[1]
		aberr = row[2]
		band = row[4]
		instrument = row[5]
		csvout.writerow(['photometry', 'MJD', mjd, band, instrument, abmag, aberr, 0])
	
	outfile.close()

# Now import the Asiago catalog
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
		eventname = snname("SN" + record[1])
		hostname = record[2]
		redvel = record[11].strip(':')
		hvel = ''
		redshift = ''
		if redvel != '':
			if round(float(redvel)) == float(redvel):
				hvel = int(redvel)
			else:
				redshift = float(redvel)

		claimedtype = record[17].strip(':')

		if (eventname in eventnames):
			outfile = open(outdir + eventname + '.dat', 'a')
			csvout = csv.writer(outfile, quotechar='"', quoting=csv.QUOTE_ALL, delimiter="\t")
		else:
			eventnames.append(eventname)
			outfile = open(outdir + eventname + '.dat', 'wb')
			csvout = csv.writer(outfile, quotechar='"', quoting=csv.QUOTE_ALL, delimiter="\t")

			csvout.writerow(['name', eventname])

		if (hostname != ''):
			csvout.writerow(['host', hostname])
		if (claimedtype != ''):
			csvout.writerow(['claimedtype', claimedtype])
		if (redshift != ''):
			csvout.writerow(['redshift', redshift])
		if (hvel != ''):
			csvout.writerow(['hvel', hvel])

		outfile.close()
