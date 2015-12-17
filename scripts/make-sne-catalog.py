#!/usr/local/bin/python2.7

import csv
import glob
import sys
import os
import re
from bokeh.io import hplot, vplot
from bokeh.plotting import figure, show, save
from bokeh.resources import CDN
from bokeh.embed import file_html

outdir = "../"

header = [
	"Name(s)",
	"Host Name(s)",
	"Publications",
	"Instruments",
	"$z$",
	r"$v_{\rm Helio}$",
#	"$N_{\\rm h}$",
	"Claimed Type",
	"Notes",
	"Plot(s)"
	]

columnkey = [
	"name",
	"host",
	"citations",
	"instruments",
	"redshift",
	"hvel",
#	"nh",
	"claimedtype",
	"notes",
	"plot"
	]

footer = [
	"Note: Name of transient preferred, if unnamed host galaxy is used",
	"*&nbsp;Uncertain\nNote: Usually host galaxy name(s)",
	"* discovery\n&Dagger; sne type first proposed",
	"&Dagger; upper limit\n* host only\n** host and flare",
	"*&nbsp;Uncertain",
	"",
#	"Line of sight H column",
	"Note: Papers claiming a sne without specifying disruptor/disruptee are listed as \"sne.\"",
	"",
	""
	]

if (len(columnkey) != len(header) or len(columnkey) != len(footer)):
	print 'Error: Header and/or footer not same length as key list.'
	sys.exit(0)

dataavaillink = "<a href='https://bitbucket.org/Guillochon/sne'>Y</a>";

csvout = open(outdir + 'sne-catalog.csv', 'wb')
csvout = csv.writer(csvout, quotechar='"', quoting=csv.QUOTE_ALL)

header = dict(zip(columnkey,header))
csvout.writerow([header[x] for x in columnkey])

bandcolors = [
	"indigo",
	"darkblue",
	"mediumvioletred",
	"pink",
	"#d82930"
]

bandcodes = [
	"U",
	"B",
	"V",
	"R",
	"I"
]

bandnames = [
	"U",
	"B",
	"V",
	"R",
	"I"
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

for file in sorted(glob.glob(outdir + "*.dat"), key=lambda s: s.lower()):
	tsvin = open(file,'rb')
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

	plotavail = False;
	for row in tsvin:
		if row[0] == 'photometry':
			plotavail = True;
			eventname = os.path.splitext(file)[0]
			plotlink = "<a href='https://sne.space/sne/" + eventname + ".html'><img src='https://sne.space/light-curve-icon.png'></a>";
			catalog['plot'] = plotlink

			phototu.append(row[1])
			phototime.append(float(row[2]))
			photoband.append(row[3])
			photoinstrument.append(row[4].strip())
			photoAB.append(float(row[5]))
			photoerr.append(float(row[6]))
			phototype.append(int(row[7]))

		if row[0] in columnkey:
			table.append(row)
			catalog[row[0]] = row[1]
	
	instrulist = sorted(filter(None, list(set(photoinstrument))))
	instruments = ", ".join(instrulist)
	if len(instrulist) > 0:
		catalog['instruments'] = instruments

	csvout.writerow([catalog[x] for x in columnkey])

	tools = "pan,wheel_zoom,box_zoom,save,crosshair,hover,reset,resize"

	if plotavail:
		eventname = os.path.splitext(file)[0]

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
		html = re.sub(r'(\<\/body\>)', r'<a href="https://sne.space/sne/' + eventname + r'.dat">Download datafile</a><br><br>\n	  \1', html)
		html = re.sub(r'(\<\/body\>)', returnlink+r'\n	  \1', html)
		with open(outdir + eventname + ".html", "w") as f:
			f.write(html)

footer = dict(zip(columnkey,footer))
csvout.writerow([footer[x] for x in columnkey])

