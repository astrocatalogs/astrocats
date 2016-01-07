#!/usr/local/bin/python2.7

import csv
import glob
import sys
import os
import re
import bz2
import operator
import json
from colorpy.ciexyz import xyz_from_wavelength
from colorpy.colormodels import irgb_string_from_xyz
from random import shuffle, seed
from collections import OrderedDict
from bokeh.plotting import figure, show, save, ColumnDataSource
from bokeh.models import HoverTool
from bokeh.resources import CDN
from bokeh.embed import file_html

writecatalog = True
indir = "../data/"
outdir = "../"

columnkey = [
    "num",
    "name",
    "aliases",
    "discoveryear",
    "discovermonth",
    "discoverday",
    "discoverdate",
    "maxyear",
    "maxmonth",
    "maxday",
    "maxdate",
    "host",
    "citations",
    "instruments",
    "redshift",
    "hvel",
    "claimedtype",
    "notes",
    "numphoto",
    "plot",
    "data"
]

header = [
    "",
    "Name",
    "Other Names",
    "Discovery Year",
    "Discovery Month",
    "Discovery Day",
    "Discovery Date",
    "Year of Maximum",
    "Month of Maximum",
    "Day of Maximum",
    "Date of Maximum",
    "Host Name",
    "Publications",
    "Instruments/Bands",
    "<em>z</em>",
    r"<em>v</em><sub>Helio</sub>",
    "Claimed Type",
    "Notes",
    "# Photo Pts.",
    "Plot",
    "Data"
]

showcols = [
    True,
    True,
    True,
    False,
    False,
    False,
    True,
    False,
    False,
    False,
    True,
    True,
    False,
    True,
    True,
    True,
    True,
    False,
    False,
    False,
    True
]

photokeys = [
    'timeunit',
    'time',
    'band',
    'instrument',
    'abmag',
    'aberr',
    'upperlimit',
    'source'
]

sourcekeys = [
    'name',
    'alias',
    'secondary'
]

if len(columnkey) != len(header):
    print 'Error: Header not same length as key list.'
    sys.exit(0)

dataavaillink = "<a href='https://bitbucket.org/Guillochon/sne'>Y</a>";

header = dict(zip(columnkey,header))

bandcodes = [
    "u",
    "g",
    "r",
    "i",
    "z",
    "u'",
    "g'",
    "r'",
    "i'",
    "z'",
    "u_SDSS",
    "g_SDSS",
    "r_SDSS",
    "i_SDSS",
    "z_SDSS",
    "U",
    "B",
    "V",
    "R",
    "I",
    "G",
    "Y",
    "J",
    "H",
    "K",
    "C"
]

bandaliases = {
    "u_SDSS" : "u (SDSS)",
    "g_SDSS" : "g (SDSS)",
    "r_SDSS" : "r (SDSS)",
    "i_SDSS" : "i (SDSS)",
    "z_SDSS" : "z (SDSS)"
}

bandshortaliases = {
    "u_SDSS" : "u",
    "g_SDSS" : "g",
    "r_SDSS" : "r",
    "i_SDSS" : "i",
    "z_SDSS" : "z",
    "G" : ""
}

bandwavelengths = {
    "u" : 354.,
    "g" : 475.,
    "r" : 622.,
    "i" : 763.,
    "z" : 905.,
    "u'" : 354.,
    "g'" : 475.,
    "r'" : 622.,
    "i'" : 763.,
    "z'" : 905.,
    "u_SDSS" : 354.3,
    "g_SDSS" : 477.0,
    "r_SDSS" : 623.1,
    "i_SDSS" : 762.5,
    "z_SDSS" : 913.4,
    "U" : 365.,
    "B" : 445.,
    "V" : 551.,
    "R" : 658.,
    "I" : 806.,
    "Y" : 1020.,
    "J" : 1220.,
    "H" : 1630.,
    "K" : 2190.
}

wavedict = dict(zip(bandcodes,bandwavelengths))

seed(101)
bandcolors = ["#%06x" % round(float(x)/float(len(bandcodes))*0xFFFEFF) for x in xrange(len(bandcodes))]
shuffle(bandcolors)

def event_filename(name):
    return(name.replace('/', '_'))

# Replace bands with real colors, if possible.
for b, code in enumerate(bandcodes):
    if (code in bandwavelengths):
        hexstr = irgb_string_from_xyz(xyz_from_wavelength(bandwavelengths[code]))
        if (hexstr != "#000000"):
            bandcolors[b] = hexstr

bandcolordict = dict(zip(bandcodes,bandcolors))

coldict = dict(zip(range(len(columnkey)),columnkey))

def bandcolorf(color):
    if (color in bandcolordict):
        return bandcolordict[color]
    return 'black'

def bandaliasf(code):
    if (code in bandaliases):
        return bandaliases[code]
    return code

def bandshortaliasf(code):
    if (code in bandshortaliases):
        return bandshortaliases[code]
    return code

def bandwavef(code):
    if (code in bandwavelengths):
        return bandwavelengths[code]
    return 0.

catalogrows = []
sourcerows = []
for fcnt, file in enumerate(sorted(glob.glob(indir + "*.dat"), key=lambda s: s.lower())):
    #sorted(glob.glob(indir + "*.bz2"), key=lambda s: s.lower()) + 
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
    photometry = []
    sources = []

    eventname = os.path.basename(file).split('.')[0]

    plotavail = False;
    for row in tsvin:
        photorow = OrderedDict.fromkeys(photokeys, '')
        sourcerow = OrderedDict.fromkeys(sourcekeys, '')
        if row[0] == 'photometry':
            plotavail = True;
            plotlink = "sne/" + eventname + ".html";
            catalog['plot'] = plotlink
            plotlink = "<a class='lci' href='" + plotlink + "' target='_blank'></a>";

            photodict = dict(zip(row[1:], row[2:]))

            for key in photorow:
                if key in photodict:
                    photorow[key] = photodict[key]

            photometry.append(photorow)

        elif row[0] == 'source':
            sourcedict = dict(zip(row[1:], row[2:]))

            for key in sourcerow:
                if key in sourcedict:
                    sourcerow[key] = sourcedict[key]

            sources.append(sourcerow)

        elif row[0] in columnkey:
            table.append(row)
            catalog[row[0]] = row[1]

    catalog['data'] = "<span class='ics'>"
    catalog['data'] += "<a class='dci' href='sne/data/" + eventname + ".dat.bz2'></a>"
    if plotavail:
        catalog['data'] += plotlink
    catalog['data'] += "</span>"
    
    catalog['numphoto'] = len(photometry)
    prange = xrange(len(photometry))
    instrulist = sorted(filter(None, list(set([photometry[x]['instrument'] for x in prange]))))
    if len(instrulist) > 0:
        instruments = ''
        for i, instru in enumerate(instrulist):
            instruments += instru
            bandlist = sorted(filter(None, list(set([bandshortaliasf(photometry[x]['band'])
                if photometry[x]['instrument'] == instru else "" for x in prange]))), key=lambda y: bandwavef(y))
            if bandlist:
                instruments += ' (' + ", ".join(bandlist) + ')'
            if i < len(instrulist) - 1:
                instruments += ', '

        catalog['instruments'] = instruments
    else:
        bandlist = sorted(filter(None, list(set([bandshortaliasf(photometry[x]['band']) for x in prange]))), key=lambda y: bandwavef(y))
        if len(bandlist) > 0:
            catalog['instruments'] = ", ".join(bandlist)

    catalogrows.append(catalog)
    sourcerows.append(sources)

    tools = "pan,wheel_zoom,box_zoom,save,crosshair,reset,resize"
    hover = HoverTool(
        tooltips=[
            ("MJD", "@x{1.11}"),
            ("Magnitude", "@y{1.111}"),
            ("Error", "@err{1.111}"),
            ("Instrument", "@instr"),
            ("Band", "@desc")
        ]
    )

    if plotavail:
        phototime = [float(photometry[x]['time']) for x in prange]
        photoAB = [float(photometry[x]['abmag']) for x in prange]
        photoerrs = [float(photometry[x]['aberr'] if photometry[x]['aberr'] else 0.) for x in prange]

        x_buffer = 0.1*(max(phototime) - min(phototime)) if len(phototime) > 1 else 1.0
        x_range = [-x_buffer + min(phototime), x_buffer + max(phototime)]

        p1 = figure(title='Photometry for ' + eventname, x_axis_label='Time (' + photometry[0]['timeunit'] + ')',
            y_axis_label='AB Magnitude', x_range = x_range, tools = tools,
            y_range = (0.5 + max([x + y for x, y in zip(photoAB, photoerrs)]), -0.5 + min([x - y for x, y in zip(photoAB, photoerrs)])))
        p1.add_tools(hover)

        err_xs = []
        err_ys = []

        for x, y, yerr in zip(phototime, photoAB, photoerrs):
            err_xs.append((x, x))
            err_ys.append((y - yerr, y + yerr))

        photoband = [photometry[x]['band'] for x in prange]
        photoinstru = [photometry[x]['instrument'] for x in prange]
        phototype = [int(photometry[x]['upperlimit']) if photometry[x]['upperlimit'] else 0 for x in prange]
        bandset = set(photoband)
        bandset = [i for (j, i) in sorted(zip(map(bandaliasf, bandset), bandset))]

        for band in bandset:
            bandname = bandaliasf(band)
            indb = [i for i, j in enumerate(photoband) if j == band]
            indt = [i for i, j in enumerate(phototype) if j == 0]
            ind = set(indb).intersection(indt)

            source = ColumnDataSource(
                data = dict(
                    x = [phototime[i] for i in ind],
                    y = [photoAB[i] for i in ind],
                    err = [photoerrs[i] for i in ind],
                    desc = [photoband[i] for i in ind],
                    instr = [photoinstru[i] for i in ind]
                )
            )
            p1.circle('x', 'y', source = source, color=bandcolorf(band), legend=bandname, size=5)
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
        #html = re.sub(r'(\<body\>)', r'\1\n    '+returnlink, html)
        html = re.sub(r'(\<\/body\>)', r'    <a href="data/' + eventname + r'.dat.bz2">Download datafile</a><br><br>\n        \1', html)
        if len(sources):
            html = re.sub(r'(\<\/body\>)', r'<em>Sources of data:</em><br>\n        \1', html)
            for source in sources:
                html = re.sub(r'(\<\/body\>)', source['name']+r'<br>\n        \1', html)
            html = re.sub(r'(\<\/body\>)', r'<br>\n    \1', html)
        html = re.sub(r'(\<\/body\>)', returnlink+r'\n    \1', html)
        print outdir + eventname + ".html"
        with open(outdir + eventname + ".html", "w") as f:
            f.write(html)
    #if fcnt > 100:
    #    break

# Construct the date
for r, row in enumerate(catalogrows):
    year = None
    month = None
    day = None
    
    discoverdatestr = ''
    if catalogrows[r]['discoveryear']:
        discoverdatestr += catalogrows[r]['discoveryear']
        if catalogrows[r]['discovermonth']:
            discoverdatestr += '-' + catalogrows[r]['discovermonth'].zfill(2)
            if catalogrows[r]['discoverday']:
                discoverdatestr += '-' + catalogrows[r]['discoverday'].zfill(2)
    catalogrows[r]['discoverdate'] = discoverdatestr

    maxdatestr = ''
    if catalogrows[r]['maxyear']:
        maxdatestr += catalogrows[r]['maxyear']
        if catalogrows[r]['maxmonth']:
            maxdatestr += '-' + catalogrows[r]['maxmonth'].zfill(2)
            if catalogrows[r]['maxday']:
                maxdatestr += '-' + catalogrows[r]['maxday'].zfill(2)
    catalogrows[r]['maxdate'] = maxdatestr

# Write it all out at the end
if writecatalog:
    f = open(outdir + 'sne-catalog.csv', 'wb')
    csvout = csv.writer(f, quotechar='"', quoting=csv.QUOTE_ALL)

    prunedheader = [header[coldict[i]] for (i, j) in enumerate(showcols) if j]
    csvout.writerow(prunedheader)

    catalogrows.sort(key=operator.itemgetter('discoverdate'), reverse=True)
    prunedrows = []
    for r, row in enumerate(catalogrows):
        prunedrow = [row[coldict[i]] for (i, j) in enumerate(showcols) if j]
        csvout.writerow(prunedrow)
        prunedrows.append(prunedrow)

    prunedfooter = [header[coldict[i]] for (i, j) in enumerate(showcols) if j]
    csvout.writerow(prunedfooter)
    f.close()

    jsonobj = dict.fromkeys(['data'])
    jsonobj['data'] = prunedrows
    jsonstring = json.dumps(jsonobj, indent=4, separators=(',', ': '))
    f = open(outdir + 'sne-catalog.json', 'wb')
    f.write(jsonstring)
    f.close()

    f = open(outdir + 'catalog.html', 'wb')
    f.write('<table id="example" class="display" cellspacing="0" width="100%">\n')
    f.write('\t<thead>\n')
    f.write('\t\t<tr>\n')
    for i, j in enumerate(showcols):
        if j:
            f.write('\t\t\t<th class="' + coldict[i] + '">' + header[coldict[i]] + '</th>\n')
    f.write('\t\t</tr>\n')
    f.write('\t</thead>\n')
    f.write('\t<tfoot>\n')
    f.write('\t\t<tr>\n')
    for i, j in enumerate(showcols):
        if j:
            f.write('\t\t\t<th>' + header[coldict[i]] + '</th>\n')
    f.write('\t\t</tr>\n')
    f.write('\t</thead>\n')
    f.write('</table>\n')
    f.close()

# Make a few small files for generating charts
f = open(outdir + 'snepages.csv', 'wb')
csvout = csv.writer(f, quotechar='"', quoting=csv.QUOTE_ALL)
for row in catalogrows:
    if row['plot']:
        csvout.writerow([row['aliases'], 'https://sne.space/' + row['plot']])
f.close()

f = open(outdir + 'sources.csv', 'wb')
sourcedict = dict()
for sources in sourcerows:
    for sourcerow in sources:
        strippedname = re.sub('<[^<]+?>', '', sourcerow['name'])
        if strippedname in sourcedict:
            sourcedict[strippedname] += 1
        else:
            sourcedict[strippedname] = 1

sortedsources = sorted(sourcedict.items(), key=operator.itemgetter(1), reverse=True)
csvout = csv.writer(f)
csvout.writerow(['Source','Number'])
for source in sortedsources:
    csvout.writerow(source)
f.close()

nophoto = sum(i < 3 for i in [x['numphoto'] for x in catalogrows])
hasphoto = len(catalogrows) - nophoto
f = open(outdir + 'pie.csv', 'wb')
csvout = csv.writer(f)
csvout.writerow(['Category','Number'])
csvout.writerow(['Has light curve', hasphoto])
csvout.writerow(['No light curve', nophoto])
f.close()
f = open(outdir + 'hasphoto.html', 'wb')
f.write(str(hasphoto))
f.close()

ctypedict = dict()
for row in catalogrows:
    cleanedtype = row['claimedtype'].strip('?* ')
    cleanedtype = cleanedtype.replace('Ibc', 'Ib/c')
    cleanedtype = cleanedtype.replace('IIP', 'II P')
    if not cleanedtype:
        cleanedtype = 'Unknown'
    if cleanedtype in ctypedict:
        ctypedict[cleanedtype] += 1
    else:
        ctypedict[cleanedtype] = 1
sortedctypes = sorted(ctypedict.items(), key=operator.itemgetter(1), reverse=True)
f = open(outdir + 'types.csv', 'wb')
csvout = csv.writer(f)
csvout.writerow(['Type','Number'])
for ctype in sortedctypes:
    csvout.writerow(ctype)
f.close()

years = filter(None, [int(x['discoveryear']) if x['discoveryear'] else '' for x in catalogrows])
yearrange = range(min(years), max(years))
f = open(outdir + 'area.csv', 'wb')
csvout = csv.writer(f)
csvout.writerow(['Year','Has light curve','No light curve'])
csvout.writerow(['date','number','number'])
for year in yearrange:
    yearind = [i for i, x in enumerate(catalogrows) if x['discoveryear'] and int(x['discoveryear']) == year]
    hasphoto = len(yearind) - sum(i < 3 for i in [x['numphoto'] for x in [catalogrows[y] for y in yearind]])
    nophoto = sum(i < 3 for i in [x['numphoto'] for x in [catalogrows[y] for y in yearind]])
    csvout.writerow([year, hasphoto, nophoto])
f.close()
