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

header = [
	"Name(s)",
	"Host Name(s)",
	"Publications",
	"Instruments",
	"$z$",
	"$N_{\\rm h}$",
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
	"nh",
	"claimedtype",
	"notes",
	"plot"
	]

footer = [
	"Note: Name of transient preferred, if unnamed host galaxy is used",
	"*&nbsp;Uncertain\nNote: Usually host galaxy name(s)",
	"* discovery\n&Dagger; TDE type first proposed",
	"&Dagger; upper limit\n* host only\n** host and flare",
	"*&nbsp;Uncertain",
	"Line of sight H column",
	"Note: Papers claiming a TDE without specifying disruptor/disruptee are listed as \"TDE.\"",
	"",
	""
	]

if (len(columnkey) != len(header) or len(columnkey) != len(footer)):
	print 'Error: Header and/or footer not same length as key list.'
	sys.exit(0)

dataavaillink = "<a href='https://bitbucket.org/Guillochon/tde_events'>Y</a>";

csvout = open('tde-catalog.csv', 'wb')
csvout = csv.writer(csvout, quotechar='"', quoting=csv.QUOTE_ALL)

header = dict(zip(columnkey,header))
csvout.writerow([header[x] for x in columnkey])

bandcolors = [
	"pink",
	"darkorchid",
	"darkseagreen",
	"darkmagenta",
	"darkolivegreen",
	"#d82930",
	"blueviolet",
	"darkgreen",
	"darkred",
	"brown",
	"purple",
	"mediumvioletred",
	"indigo",
	"forestgreen",
	"red",
	"crimson",
	"firebrick",
	"forestgreen",
	"red",
	"crimson",
	"firebrick",
	"cornflowerblue",
	"limegreen",
	"darkblue",
	"darkslateblue",
	"midnightblue",
	"mediumvioletred",
	"mediumspringgreen",
	"orange",
	"chocolate",
	"darkorange",
	"orangered"
]

bandcodes = [
	"Cr",
	"4g",
	"4r",
	"6g",
	"6r",
	"6i",
	"Mg",
	"Mr",
	"Mi",
	"RO",
	"GN",
	"GF",
	"Sg",
	"Su",
	"Sr",
	"Si",
	"Sz",
	"Pg",
	"Pr",
	"Pi",
	"Pz",
	"Ub",
	"Uv",
	"Um",
	"U1",
	"U2",
	"Uu",
	"bV",
	"X1",
	"X2",
	"Xs",
	"Xm",
	"Xh"
]

bandnames = [
	"CSS V",
	"P48 g",
	"P48 r",
	"P60 g",
	"P60 r",
	"P60 i",
	"MegaCam g",
	"MegaCam r",
	"MegaCam i",
	"ROTSE",
	"GALEX NUV",
	"GALEX FUV",
	"SDSS g",
	"SDSS u",
	"SDSS r",
	"SDSS i",
	"SDSS z",
	"PAN-STARRS g",
	"PAN-STARRS r",
	"PAN-STARRS i",
	"PAN-STARRS z",
	"UVOT B",
	"UVOT V",
	"UVOT M2",
	"UVOT W1",
	"UVOT W2",
	"UVOT U",
	"ASAS-SN V",
	"X1",
	"X2",
	"XRT (0.2 - 0.5 keV)",
	"XRT (0.5 - 1.5 keV)",
	"Xh"
]

xraybands = [
	"X1",
	"X2",
	"Xs",
	"Xm",
	"Xh"
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

for file in sorted(glob.glob("*.dat"), key=lambda s: s.lower()):
	tsvin = open(file,'rb')
	tsvin = csv.reader(tsvin, delimiter='\t')

	catalog = dict(zip(columnkey,['' for _ in xrange(len(columnkey))]))

	table = []

	phototu = []
	phototime = []
	photoband = []
	photoAB = []
	photoerr = []
	phototype = []

	xraytu = []
	xraytime = []
	xrayband = []
	xraycts = []
	xrayerr = []
	xraytype = []

	plotavail = False;
	xrayavail = False;
	for row in tsvin:
		if row[0] == 'photometry':
			plotavail = True;
			eventname = os.path.splitext(file)[0]
			plotlink = "<a href='https://tde.space/tde_events/" + eventname + ".html'><img src='https://tde.space/plot_icon.png'></a>";
			catalog['plot'] = plotlink
			if row[3] in xraybands:
				xrayavail = True;
				xraytu.append(row[1])
				xraytime.append(float(row[2]))
				xrayband.append(row[3])
				xraycts.append(float(row[4]))
				xrayerr.append(float(row[5]))
				xraytype.append(int(row[6]))
			else:
				phototu.append(row[1])
				phototime.append(float(row[2]))
				photoband.append(row[3])
				photoAB.append(float(row[4]))
				photoerr.append(float(row[5]))
				phototype.append(int(row[6]))

		if row[0] in columnkey:
			table.append(row)
			catalog[row[0]] = row[1]
	
	csvout.writerow([catalog[x] for x in columnkey])

	tools = "pan,wheel_zoom,box_zoom,save,crosshair,hover,reset,resize"

	if plotavail:
		eventname = os.path.splitext(file)[0]

		if xrayavail:
			x_data = phototime + xraytime
		else:
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

		if xrayavail:
			p2 = figure(title='X-ray counts for ' + eventname, x_axis_label='Time (' + xraytu[0] + ')',
				y_axis_label='Log X-ray Counts', x_range = p1.x_range, tools = tools,
				y_range = (-0.5 + min([x + y for x, y in zip(xraycts, xrayerr)]), 0.5 + max([x - y for x, y in zip(xraycts, xrayerr)])))

			err_xs = []
			err_ys = []

			for x, y, yerr in zip(xraytime, xraycts, xrayerr):
				err_xs.append((x, x))
				err_ys.append((y - yerr, y + yerr))

			bandset = set(xrayband)
			bandset = [i for (j, i) in sorted(zip(map(bandnamef, bandset), bandset))]

			for band in bandset:
				bandname = bandnamef(band)
				indb = [i for i, j in enumerate(xrayband) if j == band]
				indt = [i for i, j in enumerate(xraytype) if j == 0]
				ind = set(indb).intersection(indt)

				p2.circle([xraytime[x] for x in ind], [xraycts[x] for x in ind],
					color=bandcolorf(band), legend=bandname, size=5)
				p2.multi_line([err_xs[x] for x in ind], [err_ys[x] for x in ind], color=bandcolorf(band))
				
				upplimlegend = bandname if len(ind) == 0 else ''

				indt = [i for i, j in enumerate(xraytype) if j == 1]
				ind = set(indb).intersection(indt)
				p2.inverted_triangle([xraytime[x] for x in ind], [xraycts[x] for x in ind],
					color=bandcolorf(band), legend=upplimlegend, size=8)

			p = vplot(hplot(p1, p2))
		else:
			p = p1

		#save(p)
		html = file_html(p, CDN, eventname)
		returnlink = r'    <a href="https://tde.space"><< Return to TDE catalog</a>';
		html = re.sub(r'(\<body\>)', r'\1\n    '+returnlink, html)
		html = re.sub(r'(\<\/body\>)', r'<a href="https://tde.space/tde_events/' + eventname + r'.dat">Download datafile</a><br><br>\n    \1', html)
		html = re.sub(r'(\<\/body\>)', returnlink+r'\n    \1', html)
		with open(eventname + ".html", "w") as f:
			f.write(html)

footer = dict(zip(columnkey,footer))
csvout.writerow([footer[x] for x in columnkey])

