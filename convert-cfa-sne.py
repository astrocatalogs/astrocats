#!/usr/local/bin/python2.7

import csv
import glob
import os

for file in sorted(glob.glob("cfa-input/*.dat"), key=lambda s: s.lower()):
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
	eventname = (eventparts[0]).upper()
	eventbands = list(eventparts[1])
	print eventname
	print eventbands

	csvout = open(eventname + '.dat', 'wb')
	csvout = csv.writer(csvout, quotechar='"', quoting=csv.QUOTE_ALL, delimiter="\t")

	csvout.writerow(['name', eventname])

	for rc, row in enumerate(csv_data):
		tu = 'MJD'
		if len(row) > 0 and row[0][0] == "#":
			if len(row) > 1 and row[1] == 'HJD':
				tu = 'HJD'
			continue
		elif len(row) > 0:
			mjd = row[0]
			for v, val in enumerate(row):
				if v == 0:
					mjd = val
				elif v % 2 != 0:
					csvout.writerow(['photometry', tu, mjd, eventbands[(v-1)/2], '', row[v], row[v+1], 0])
