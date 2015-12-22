#!/usr/local/bin/python2.7

import csv
import glob
import sys
import os
import re
import bz2
from bokeh.io import hplot, vplot
from bokeh.plotting import figure, show, save
from bokeh.resources import CDN
from bokeh.embed import file_html

indir = "../data/"
outdir = "../"

header = [
	"Names",
	"Host Names",
#	"Publications",
	"Instruments/Surveys",
	"<em>z</em>",
	r"<em>v</em><sub>Helio</sub>",
#	"$N_{\\rm h}$",
	"Claimed Type",
#	"Notes",
	"Plots",
	"Data"
	]

columnkey = [
	"name",
	"host",
#	"citations",
	"instruments",
	"redshift",
	"hvel",
#	"nh",
	"claimedtype",
#	"notes",
	"plot",
	"data"
	]

footer = [
	"Note: IAU name preferred",
	"*&nbsp;Uncertain",
#	"* discovery\n&Dagger; sne type first proposed",
	"",
	"",
	"",
#	"Line of sight H column",
	"",
#	"",
	"",
	""
	]

if (len(columnkey) != len(header) or len(columnkey) != len(footer)):
	print 'Error: Header and/or footer not same length as key list.'
	sys.exit(0)

dataavaillink = "<a href='https://bitbucket.org/Guillochon/sne'>Y</a>";

header = dict(zip(columnkey,header))
headerrow = [header[x] for x in columnkey]

bandcolors = [
	"indigo",
	"firebrick",
	"forestgreen",
	"red",
	"crimson",
	"indigo",
	"darkblue",
	"mediumvioletred",
	"pink",
	"#d82930",
	"orangered",
	"mediumvioletred",
	"mediumspringgreen",
	"orange",
	"chocolate",
	"darkorange",
]

bandcodes = [
	"u",
	"g",
	"r",
	"i",
	"z",
	"U",
	"B",
	"V",
	"R",
	"I",
	"G",
	"Y",
	"J",
	"H",
	"K"
]

bandnames = [
	"u",
	"g",
	"r",
	"i",
	"z",
	"U",
	"B",
	"V",
	"R",
	"I",
	"G",
	"Y",
	"J",
	"H",
	"K"
]

bandcolordict = dict(zip(bandcodes,bandcolors))
bandnamedict = dict(zip(bandcodes,bandnames))

def bandcolorf(color):
	if (color in bandcolordict):
		return bandcolordict[color]
	else:
		return 'black'

def bandnamef(code):
	if (code in bandnamedict):
		return bandnamedict[code]
	else:
		return code

catalogrows = []
for file in (sorted(glob.glob(indir + "*.bz2"), key=lambda s: s.lower()) + sorted(glob.glob(indir + "*.dat"), key=lambda s: s.lower())):
	print file
	filehead, ext = os.path.splitext(file)
	if ext == ".dat":
		tsvin = open(file,'rb')
	elif ext == ".bz2":
		tsvin = bz2.BZ2File(file,'rb')
	else:
		print "illegal file extension"
	tsvin = csv.reader(tsvin, delimiter='\t')

	catalog = dict(zip(columnkey,['' for _ in xrange(len(columnkey))]))

	table = []

	phototu = []
	phototime = []
	photoband = []
	photoinstrument = []
	photoAB = []
	photoerr = []
	phototype = []

	eventname = os.path.basename(file).split('.')[0]

	plotavail = False;
	for row in tsvin:
		if row[0] == 'photometry':
			plotavail = True;
			plotlink = "<a href='https://sne.space/sne/" + eventname + ".html' target='_blank'><img alt='plot' width='32' height='32' src='https://sne.space/light-curve-icon.png'></a>";
			catalog['plot'] = plotlink

			phototu.append(row[1])
			phototime.append(float(row[2]))
			photoband.append(row[4])
			photoinstrument.append(row[6].strip())
			photoAB.append(float(row[8]))
			if not row[10]:
				photoerr.append(0.)
			else:
				photoerr.append(float(row[10]))
			phototype.append(int(row[12]))

		if row[0] in columnkey:
			table.append(row)
			catalog[row[0]] = row[1]

	catalog['data'] = r'<a href="https://sne.space/sne/' + eventname + r'.dat">Download</a>'
	
	instrulist = sorted(filter(None, list(set(photoinstrument))))
	instruments = ", ".join(instrulist)
	if len(instrulist) > 0:
		catalog['instruments'] = instruments

	catalogrows.append([catalog[x] for x in columnkey])

	tools = "pan,wheel_zoom,box_zoom,save,crosshair,hover,reset,resize"

	if plotavail:
		x_data = phototime

		x_buffer = 0.1*(max(x_data) - min(x_data))
		x_range = [-x_buffer + min(x_data), x_buffer + max(x_data)]

		p1 = figure(title='Photometry for ' + eventname, x_axis_label='Time (' + phototu[0] + ')',
			y_axis_label='AB Magnitude', x_range = x_range, tools = tools,
			y_range = (0.5 + max([x + y for x, y in zip(photoAB, photoerr)]), -0.5 + min([x - y for x, y in zip(photoAB, photoerr)])))

		err_xs = []
		err_ys = []

		for x, y, yerr in zip(phototime, photoAB, photoerr):
			err_xs.append((x, x))
			err_ys.append((y - yerr, y + yerr))

		bandset = set(photoband)
		bandset = [i for (j, i) in sorted(zip(map(bandnamef, bandset), bandset))]

		for band in bandset:
			bandname = bandnamef(band)
			indb = [i for i, j in enumerate(photoband) if j == band]
			indt = [i for i, j in enumerate(phototype) if j == 0]
			ind = set(indb).intersection(indt)

			p1.circle([phototime[x] for x in ind], [photoAB[x] for x in ind],
				color=bandcolorf(band), legend=bandname, size=5)
			p1.multi_line([err_xs[x] for x in ind], [err_ys[x] for x in ind], color=bandcolorf(band))

			upplimlegend = bandname if len(ind) == 0 else ''

			indt = [i for i, j in enumerate(phototype) if j == 1]
			ind = set(indb).intersection(indt)
			p1.inverted_triangle([phototime[x] for x in ind], [photoAB[x] for x in ind],
				color=bandcolorf(band), legend=upplimlegend, size=8)

		p = p1

		#save(p)
		html = file_html(p, CDN, eventname)
		returnlink = r'    <a href="https://sne.space"><< Return to supernova catalog</a>';
		html = re.sub(r'(\<body\>)', r'\1\n    '+returnlink, html)
		html = re.sub(r'(\<\/body\>)', r'<a href="https://sne.space/sne/data/' + eventname + r'.dat.bz2">Download datafile</a><br><br>\n	  \1', html)
		html = re.sub(r'(\<\/body\>)', returnlink+r'\n	  \1', html)
		print outdir + eventname + ".html"
		with open(outdir + eventname + ".html", "w") as f:
			f.write(html)

# Write it all out at the end
csvout = open(outdir + 'sne-catalog.csv', 'wb')
csvout = csv.writer(csvout, quotechar='"', quoting=csv.QUOTE_ALL)

csvout.writerow(headerrow)

for row in catalogrows:
	csvout.writerow(row)

footer = dict(zip(columnkey,footer))
csvout.writerow([footer[x] for x in columnkey])
