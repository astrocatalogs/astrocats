#!/usr/local/bin/python3.5

import csv
import glob
import sys
import os
import re
import operator
import json
import argparse
import hashlib
import numpy
import shutil
import gzip
import requests
import urllib.request
import urllib.parse
import filecmp
from photometry import *
from digits import *
from datetime import datetime
from astropy.time import Time as astrotime
from astropy.coordinates import SkyCoord as coord
from astropy import units as un
from copy import deepcopy
from collections import OrderedDict
from bokeh.plotting import Figure, show, save, reset_output
from bokeh.models import (HoverTool, CustomJS, Slider, ColumnDataSource,
                          HBox, VBox, Range1d, LinearAxis, DatetimeAxis)
from bokeh.resources import CDN, INLINE
from bokeh.embed import file_html, components
from palettable import cubehelix
from bs4 import BeautifulSoup, Tag, NavigableString
from math import isnan, floor, ceil

parser = argparse.ArgumentParser(description='Generate a catalog JSON file and plot HTML files from SNE data.')
parser.add_argument('--no-write-catalog', '-nwc', dest='writecatalog', help='Don\'t write catalog file',          default=True, action='store_false')
parser.add_argument('--no-write-html', '-nwh',    dest='writehtml',    help='Don\'t write html plot files',       default=True, action='store_false')
parser.add_argument('--no-collect-hosts', '-nch', dest='collecthosts', help='Don\'t collect host galaxy images',  default=True, action='store_false')
parser.add_argument('--force-html', '-fh',        dest='forcehtml',    help='Force write html plot files',        default=False, action='store_true')
parser.add_argument('--event-list', '-el',        dest='eventlist',    help='Process a list of events',           default=[], type=str, nargs='+')
parser.add_argument('--test', '-te',              dest='test',         help='Test this script',                   default=False, action='store_true')
parser.add_argument('--travis', '-tr',            dest='travis',       help='Set some options when using Travis', default=False, action='store_true')
args = parser.parse_args()

outdir = "../"

travislimit = 1000

googlepingurl = "http://www.google.com/webmasters/tools/ping?sitemap=https%3A%2F%2Fsne.space%2Fsitemap.xml"

linkdir = "https://sne.space/sne/"

testsuffix = '.test' if args.test else ''

mycolors = cubehelix.perceptual_rainbow_16.hex_colors[:14]

columnkey = [
    "check",
    "name",
    "aliases",
    "discoverdate",
    "maxdate",
    "maxappmag",
    "maxabsmag",
    "host",
    "ra",
    "dec",
    "instruments",
    "redshift",
    "velocity",
    "lumdist",
    "claimedtype",
    "photolink",
    "spectralink",
    "references",
    "download",
    "responsive"
]

eventignorekey = [
    "download"
]

header = [
    "",
    "Name",
    "Aliases",
    "Disc. Date",
    "Max Date",
    r"<em>m</em><sub>max</sub>",
    r"<em>M</em><sub>max</sub>",
    "Host Name",
    "R.A. (h:m:s)",
    "Dec. (d:m:s)",
    "Instruments/Bands",
    r"<em>z</em>",
    r"<em>v</em><sub>&#9737;</sub> (km/s)",
    r"<em>d</em><sub>L</sub> (Mpc)",
    "Claimed Type",
    "Phot.",
    "Spec.",
    "References",
    "Data",
    ""
]

eventpageheader = [
    "",
    "Name",
    "Aliases",
    "Discovery Date",
    "Maximum Date [band]",
    r"<em>m</em><sub>max</sub> [band]",
    r"<em>M</em><sub>max</sub> [band]",
    "Host Name",
    "R.A. (h:m:s)",
    "Dec. (d:m:s)",
    "Instruments/Bands",
    r"<em>z</em>",
    r"<em>v</em><sub>&#9737;</sub> (km/s)",
    r"<em>d</em><sub>L</sub> (Mpc)",
    "Claimed Type",
    "# Phot. Obs.",
    "# Spectra",
    "References",
    "Download",
    ""
]

titles = [
    "",
    "Name (IAU name preferred)",
    "Aliases",
    "Discovey Date (year-month-day)",
    "Date of Maximum (year-month-day)",
    "Maximum apparent AB magnitude",
    "Maximum absolute AB magnitude",
    "Host Name",
    "J2000 Right Ascension (h:m:s)",
    "J2000 Declination (d:m:s)",
    "List of Instruments and Bands",
    "Redshift",
    "Heliocentric velocity (km/s)",
    "Luminosity distance (Mpc)",
    "Claimed Type",
    "Photometry",
    "Spectra",
    "Bibcodes of references with most data on event",
    "Download and edit data",
    ""
]

photokeys = [
    'timeunit',
    'time',
    'band',
    'instrument',
    'magnitude',
    'aberr',
    'upperlimit',
    'source'
]

sourcekeys = [
    'name',
    'alias',
    'secondary'
]

newfiletemplate = (
'''{
\t"{0}":{
\t\t"name":"{0}",
\t\t"aliases":[
\t\t\t"{0}"
\t\t]
\t}
}'''
)

sitemaptemplate = (
'''<?xml version="1.0" encoding="UTF-8"?>
<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9"> 
  <url>
    <loc>https://sne.space</loc>
    <priority>1.0</priority>
    <changefreq>daily</changefreq>
  </url>
  <url>
    <loc>https://sne.space/about</loc>
    <priority>0.8</priority>
  </url>
  <url>
    <loc>https://sne.space/contribute</loc>
    <priority>0.7</priority>
  </url>
  <url>
    <loc>https://sne.space/derivations</loc>
    <priority>0.7</priority>
  </url>
  <url>
    <loc>https://sne.space/statistics</loc>
    <priority>0.7</priority>
  </url>
  <url>
    <loc>https://sne.space/download</loc>
    <priority>0.7</priority>
  </url>
  <url>
    <loc>https://sne.space/links</loc>
    <priority>0.7</priority>
  </url>
{0}</urlset>'''
)

with open('rep-folders.txt', 'r') as f:
    repfolders = f.read().splitlines()

repyears = [int(repfolders[x][-4:]) for x in range(len(repfolders))]
repyears[0] -= 1

if len(columnkey) != len(header):
    raise(ValueError('Header not same length as key list.'))
    sys.exit(0)

if len(columnkey) != len(eventpageheader):
    raise(ValueError('Event page header not same length as key list.'))
    sys.exit(0)

dataavaillink = "<a href='https://bitbucket.org/Guillochon/sne'>Y</a>"

header = OrderedDict(list(zip(columnkey,header)))
eventpageheader = OrderedDict(list(zip(columnkey,eventpageheader)))
titles = OrderedDict(list(zip(columnkey,titles)))

wavedict = dict(list(zip(bandcodes,bandwavelengths)))

def event_filename(name):
    return(name.replace('/', '_'))

# Replace bands with real colors, if possible.
#for b, code in enumerate(bandcodes):
#    if (code in bandwavelengths):
#        hexstr = irgb_string_from_xyz(xyz_from_wavelength(bandwavelengths[code]))
#        if (hexstr != "#000000"):
#            bandcolors[b] = hexstr

coldict = dict(list(zip(list(range(len(columnkey))),columnkey)))

def utf8(x):
    return str(x, 'utf-8')

def get_rep_folder(entry):
    if 'discoverdate' not in entry:
        return repfolders[0]
    if not is_number(entry['discoverdate'][0]['value'].split('/')[0]):
        raise(ValueError('Discovery year is not a number!'))
        sys.exit()
    for r, repyear in enumerate(repyears):
        if int(entry['discoverdate'][0]['value'].split('/')[0]) <= repyear:
            return repfolders[r]
    return repfolders[0]

def label_format(label):
    newlabel = label.replace('Angstrom', 'Å')
    newlabel = newlabel.replace('^2', '²')
    return newlabel

def is_valid_link(url):
    response = requests.get(url)
    try:
        response.raise_for_status()
    except:
        return False
    return True

def get_first_value(name, field):
    return catalog[name][field][0]['value'] if field in catalog[name] and catalog[name][field] else ''

def get_first_kind(name, field):
    return (catalog[name][field][0]['kind'] if field in catalog[name] and
        catalog[name][field] and 'kind' in catalog[name][field][0] else '')

catalog = OrderedDict()
catalogcopy = OrderedDict()
snepages = [["# name", "aliases", "max apparent mag", "max mag date", "claimed type", "redshift", "redshift kind",
    "ra", "dec", "# of photometric obs.", "URL"]]
sourcedict = {}
nophoto = []
lcspye = []
lcspno = []
lconly = []
sponly = []
hasalc = []
hasasp = []
totalphoto = 0
totalspectra = 0

hostimgs = []
if os.path.isfile(outdir + 'hostimgs.json'):
    with open(outdir + 'hostimgs.json', 'r') as f:
        filetext = f.read()
    oldhostimgs = json.loads(filetext)
    oldhostimgs = [list(i) for i in zip(*oldhostimgs)]
    hostimgdict = dict(list(zip(oldhostimgs[0], oldhostimgs[1])))
else:
    hostimgdict = {}

files = []
for rep in repfolders:
    files += glob.glob('../' + rep + "/*.json") + glob.glob('../' + rep + "/*.json.gz")

md5s = []
md5 = hashlib.md5
if os.path.isfile(outdir + 'md5s.json'):
    with open(outdir + 'md5s.json', 'r') as f:
        filetext = f.read()
    oldmd5s = json.loads(filetext)
    oldmd5s = [list(i) for i in zip(*oldmd5s)]
    md5dict = dict(list(zip(oldmd5s[0], oldmd5s[1])))
else:
    md5dict = {}

for fcnt, eventfile in enumerate(sorted(files, key=lambda s: s.lower())):
    fileeventname = os.path.splitext(os.path.basename(eventfile))[0].replace('.json','')
    if args.eventlist and fileeventname not in args.eventlist:
        continue

    if args.travis and fcnt >= travislimit:
        break

    checksum = md5(open(eventfile, 'rb').read()).hexdigest()
    md5s.append([eventfile, checksum])

    if eventfile.split('.')[-1] == 'gz':
        with gzip.open(eventfile, 'rt') as f:
            filetext = f.read()
    else:
        with open(eventfile, 'r') as f:
            filetext = f.read()

    catalog.update(json.loads(filetext, object_pairs_hook=OrderedDict))
    entry = next(reversed(catalog))

    eventname = entry

    if args.eventlist and eventname not in args.eventlist:
        continue

    print(eventfile + ' [' + checksum + ']')

    repfolder = get_rep_folder(catalog[entry])
    if os.path.isfile("../sne-internal/" + fileeventname + ".json"):
        catalog[entry]['download'] = 'e'
    else:
        catalog[entry]['download'] = ''
    if 'discoverdate' in catalog[entry]:
        for d, date in enumerate(catalog[entry]['discoverdate']):
            catalog[entry]['discoverdate'][d]['value'] = catalog[entry]['discoverdate'][d]['value'].split('.')[0]
    if 'maxdate' in catalog[entry]:
        for d, date in enumerate(catalog[entry]['maxdate']):
            catalog[entry]['maxdate'][d]['value'] = catalog[entry]['maxdate'][d]['value'].split('.')[0]

    hostmag = ''
    hosterr = ''
    if 'photometry' in catalog[entry]:
        for pi, photo in enumerate(catalog[entry]['photometry']):
            if 'host' in photo and ('upperlimit' not in photo or not photo['upperlimit']):
                hostmag = float(photo['magnitude'])
                hosterr = float(photo['e_magnitude']) if 'e_magnitude' in photo else 0.0

        # Delete the host magnitudes so they are not plotted as points
        catalog[entry]['photometry'][:] = [x for x in catalog[entry]['photometry']
            if 'host' not in x]

    photoavail = 'photometry' in catalog[entry]
    # Must be two sigma above host magnitude, if host magnitude known, to add to phot count.
    numphoto = len([x for x in catalog[entry]['photometry'] if 'upperlimit' not in x and
        (not hostmag or not 'includeshost' in x or float(x['magnitude']) <= (hostmag - 2.0*hosterr))]) if photoavail else 0
    catalog[entry]['numphoto'] = numphoto

    maxabsappoffset = 0.0
    if 'maxabsmag' in catalog[entry] and 'maxappmag' in catalog[entry]:
        maxabsappoffset = float(get_first_value(entry, 'maxabsmag')) - float(get_first_value(entry, 'maxappmag'))

    plotlink = "sne/" + fileeventname + "/"
    if photoavail:
        catalog[entry]['photoplot'] = plotlink
        catalog[entry]['photolink'] = str(numphoto)
    spectraavail = 'spectra' in catalog[entry]
    catalog[entry]['numspectra'] = len(catalog[entry]['spectra']) if spectraavail else 0
    if spectraavail:
        catalog[entry]['spectraplot'] = plotlink
        catalog[entry]['spectralink'] = str(len(catalog[entry]['spectra']))

    prange = list(range(len(catalog[entry]['photometry']))) if 'photometry' in catalog[entry] else []
    
    instrulist = sorted([_f for _f in list({catalog[entry]['photometry'][x]['instrument']
        if 'instrument' in catalog[entry]['photometry'][x] else None for x in prange}) if _f])
    if len(instrulist) > 0:
        instruments = ''
        for i, instru in enumerate(instrulist):
            instruments += instru
            bandlist = sorted([_f for _f in list({bandshortaliasf(catalog[entry]['photometry'][x]['band'] if 'band' in catalog[entry]['photometry'][x] else '')
                if 'instrument' in catalog[entry]['photometry'][x] and catalog[entry]['photometry'][x]['instrument'] == instru
                else "" for x in prange}) if _f], key=lambda y: (bandwavef(y), y))
            if bandlist:
                instruments += ' (' + ", ".join(bandlist) + ')'
            if i < len(instrulist) - 1:
                instruments += ', '

        catalog[entry]['instruments'] = instruments
    else:
        bandlist = sorted([_f for _f in list({bandshortaliasf(catalog[entry]['photometry'][x]['band']
            if 'band' in catalog[entry]['photometry'][x] else '') for x in prange}) if _f], key=lambda y: (bandwavef(y), y))
        if len(bandlist) > 0:
            catalog[entry]['instruments'] = ", ".join(bandlist)

    tools = "pan,wheel_zoom,box_zoom,save,crosshair,reset,resize"

    # Check file modification times before constructing .html files, which is expensive
    dohtml = True
    if not args.forcehtml:
        if os.path.isfile(outdir + fileeventname + ".html"):
            if eventfile in md5dict and checksum == md5dict[eventfile]:
                dohtml = False

    # Copy JSON files up a directory if they've changed
    if dohtml:
        shutil.copy2(eventfile, '../' + os.path.basename(eventfile))

    if photoavail and dohtml and args.writehtml:
        phototime = [float(x['time']) for x in catalog[entry]['photometry'] if 'magnitude' in x]
        phototimelowererrs = [float(x['e_lower_time']) if ('e_lower_time' in x and 'e_upper_time' in x)
            else (float(x['e_time']) if 'e_time' in x else 0.) for x in catalog[entry]['photometry'] if 'magnitude' in x]
        phototimeuppererrs = [float(x['e_upper_time']) if ('e_lower_time' in x and 'e_upper_time' in x) in x
            else (float(x['e_time']) if 'e_time' in x else 0.) for x in catalog[entry]['photometry'] if 'magnitude' in x]
        photoAB = [float(x['magnitude']) for x in catalog[entry]['photometry'] if 'magnitude' in x]
        photoABerrs = [(float(x['e_magnitude']) if 'e_magnitude' in x else 0.) for x in catalog[entry]['photometry'] if 'magnitude' in x]
        photoband = [(x['band'] if 'band' in x else '') for x in catalog[entry]['photometry'] if 'magnitude' in x]
        photoinstru = [(x['instrument'] if 'instrument' in x else '') for x in catalog[entry]['photometry'] if 'magnitude' in x]
        photosource = [', '.join(str(j) for j in sorted(int(i) for i in catalog[entry]['photometry'][x]['source'].split(','))) for x in prange]
        phototype = [(x['upperlimit'] if 'upperlimit' in x else False) for x in catalog[entry]['photometry'] if 'magnitude' in x]

        x_buffer = 0.1*(max(phototime) - min(phototime)) if len(phototime) > 1 else 1.0

        hastimeerrs = (len(list(filter(None, phototimelowererrs))) and len(list(filter(None, phototimeuppererrs))))
        hasABerrs = len(list(filter(None, photoABerrs)))
        tt = [  
                ("Source ID(s)", "@src"),
                ("Epoch (" + catalog[entry]['photometry'][0]['timeunit'] + ")",
                 "@x{1.11}" + ("<sub>-@xle{1}</sub><sup>+@xue{1}</sup>" if hastimeerrs else ""))
             ]
        tt += [("Apparent Magnitude", "@y{1.111}" + ("&nbsp;±&nbsp;@err{1.11}" if hasABerrs else ""))]
        if 'maxabsmag' in catalog[entry] and 'maxappmag' in catalog[entry]:
            tt += [("Absolute Magnitude", "@yabs{1.111}" + ("&nbsp;±&nbsp;@err{1.11}" if hasABerrs else ""))]
        if len(list(filter(None, photoband))):
            tt += [("Band", "@desc")]
        if len(list(filter(None, photoinstru))):
            tt += [("Instrument", "@instr")]
        hover = HoverTool(tooltips = tt)

        min_x_range = -x_buffer + min([x - y for x, y in list(zip(phototime, phototimeuppererrs))])
        max_x_range = x_buffer + max([x + y for x, y in list(zip(phototime, phototimelowererrs))])
        min_y_range = 0.5 + max([x + y for x, y in list(zip(photoAB, photoABerrs))])
        max_y_range = -0.5 + min([x - y for x, y in list(zip(photoAB, photoABerrs))])

        p1 = Figure(title='Photometry for ' + eventname, x_axis_label='Time (' + catalog[entry]['photometry'][0]['timeunit'] + ')',
            y_axis_label = 'Apparent Magnitude', tools = tools, plot_width = 485, plot_height = 485, #responsive = True,
            x_range = (min_x_range, max_x_range), y_range = (min_y_range, max_y_range),
            title_text_font_size='16pt', title_text_font = 'futura')
        p1.xaxis.axis_label_text_font = 'futura'
        p1.yaxis.axis_label_text_font = 'futura'
        p1.xaxis.major_label_text_font = 'futura'
        p1.yaxis.major_label_text_font = 'futura'
        p1.xaxis.axis_label_text_font_size = '12pt'
        p1.yaxis.axis_label_text_font_size = '12pt'
        p1.xaxis.major_label_text_font_size = '8pt'
        p1.yaxis.major_label_text_font_size = '8pt'

        min_x_date = astrotime(min_x_range, format='mjd').datetime
        max_x_date = astrotime(max_x_range, format='mjd').datetime


        p1.extra_x_ranges = {"gregorian date": Range1d(start=min_x_date, end=max_x_date)}
        p1.add_layout(DatetimeAxis(axis_label = "Time (Gregorian Date)", major_label_text_font_size = '8pt',
            major_label_text_font = 'futura', axis_label_text_font = 'futura',
            x_range_name = "gregorian date", axis_label_text_font_size = '12pt'), 'above')
        if 'maxabsmag' in catalog[entry] and 'maxappmag' in catalog[entry]:
            min_y_absmag = min_y_range + maxabsappoffset
            max_y_absmag = max_y_range + maxabsappoffset
            p1.extra_y_ranges = {"abs mag": Range1d(start=min_y_absmag, end=max_y_absmag)}
            p1.add_layout(LinearAxis(axis_label = "Absolute Magnitude", major_label_text_font_size = '8pt',
                major_label_text_font = 'futura', axis_label_text_font = 'futura',
                y_range_name = "abs mag", axis_label_text_font_size = '12pt'), 'right')
        p1.add_tools(hover)

        xs = []
        ys = []
        err_xs = []
        err_ys = []

        for x, y, xlowerr, xupperr, yerr in list(zip(phototime, photoAB, phototimelowererrs, phototimeuppererrs, photoABerrs)):
            xs.append(x)
            ys.append(y)
            err_xs.append((x - xlowerr, x + xupperr))
            err_ys.append((y - yerr, y + yerr))

        bandset = set(photoband)
        bandset = [i for (j, i) in sorted(list(zip(list(map(bandaliasf, bandset)), bandset)))]

        for band in bandset:
            bandname = bandaliasf(band)
            indb = [i for i, j in enumerate(photoband) if j == band]
            indt = [i for i, j in enumerate(phototype) if not j]
            # Should always have upper error if have lower error.
            indnex = [i for i, j in enumerate(phototimelowererrs) if j == 0.]
            indyex = [i for i, j in enumerate(phototimelowererrs) if j > 0.]
            indney = [i for i, j in enumerate(photoABerrs) if j == 0.]
            indyey = [i for i, j in enumerate(photoABerrs) if j > 0.]
            indne = set(indb).intersection(indt).intersection(indney).intersection(indnex)
            indye = set(indb).intersection(indt).intersection(set(indyey).union(indyex))

            noerrorlegend = bandname if len(indne) == 0 else ''

            data = dict(
                x = [phototime[i] for i in indne],
                y = [photoAB[i] for i in indne],
                err = [photoABerrs[i] for i in indne],
                desc = [photoband[i] for i in indne],
                instr = [photoinstru[i] for i in indne],
                src = [photosource[i] for i in indne]
            )
            if 'maxabsmag' in catalog[entry] and 'maxappmag' in catalog[entry]:
                data['yabs'] = [photoAB[i] + maxabsappoffset for i in indne]
            if hastimeerrs:
                data['xle'] = [phototimelowererrs[i] for i in indne]
                data['xue'] = [phototimeuppererrs[i] for i in indne]

            source = ColumnDataSource(data)
            p1.circle('x', 'y', source = source, color=bandcolorf(band), fill_color="white", legend=noerrorlegend, size=4)

            data = dict(
                x = [phototime[i] for i in indye],
                y = [photoAB[i] for i in indye],
                err = [photoABerrs[i] for i in indye],
                desc = [photoband[i] for i in indye],
                instr = [photoinstru[i] for i in indye],
                src = [photosource[i] for i in indye]
            )
            if 'maxabsmag' in catalog[entry] and 'maxappmag' in catalog[entry]:
                data['yabs'] = [photoAB[i] + maxabsappoffset for i in indye]
            if hastimeerrs:
                data['xle'] = [phototimelowererrs[i] for i in indye]
                data['xue'] = [phototimeuppererrs[i] for i in indye]

            source = ColumnDataSource(data)
            p1.multi_line([err_xs[x] for x in indye], [[ys[x], ys[x]] for x in indye], color=bandcolorf(band))
            p1.multi_line([[xs[x], xs[x]] for x in indye], [err_ys[x] for x in indye], color=bandcolorf(band))
            p1.circle('x', 'y', source = source, color=bandcolorf(band), legend=bandname, size=4)

            upplimlegend = bandname if len(indye) == 0 and len(indne) == 0 else ''

            indt = [i for i, j in enumerate(phototype) if j]
            ind = set(indb).intersection(indt)
            data = dict(
                x = [phototime[i] for i in ind],
                y = [photoAB[i] for i in ind],
                err = [photoABerrs[i] for i in ind],
                desc = [photoband[i] for i in ind],
                instr = [photoinstru[i] for i in ind],
                src = [photosource[i] for i in ind]
            )
            if 'maxabsmag' in catalog[entry] and 'maxappmag' in catalog[entry]:
                data['yabs'] = [photoAB[i] + maxabsappoffset for i in ind]
            if hastimeerrs:
                data['xle'] = [phototimelowererrs[i] for i in ind]
                data['xue'] = [phototimeuppererrs[i] for i in ind]

            source = ColumnDataSource(data)
            # Currently Bokeh doesn't support tooltips for inverted_triangle, so hide an invisible circle behind for the tooltip
            p1.circle('x', 'y', source = source, alpha=0.0, size=7)
            p1.inverted_triangle('x', 'y', source = source,
                color=bandcolorf(band), legend=upplimlegend, size=7)

        p1.legend.label_text_font = 'futura'
        p1.legend.label_text_font_size = '8pt'
        p1.legend.label_width = 20
        p1.legend.label_height = 14
        p1.legend.glyph_height = 14

    if spectraavail and dohtml and args.writehtml:
        spectrumwave = []
        spectrumflux = []
        spectrumerrs = []
        spectrummjdmax = []
        hasepoch = True
        hasmjdmax = False
        if 'redshift' in catalog[entry]:
            z = float(catalog[entry]['redshift'][0]['value'])
        for spectrum in catalog[entry]['spectra']:
            spectrumdata = deepcopy(spectrum['data'])
            oldlen = len(spectrumdata)
            specslice = ceil(float(len(spectrumdata))/10000)
            spectrumdata = spectrumdata[::specslice]
            spectrumdata = [x for x in spectrumdata if is_number(x[1]) and not isnan(float(x[1]))]
            specrange = range(len(spectrumdata))

            if 'deredshifted' in spectrum and spectrum['deredshifted'] and 'redshift' in catalog[entry]:
                spectrumwave.append([float(spectrumdata[x][0])*(1.0 + z) for x in specrange])
            else:
                spectrumwave.append([float(spectrumdata[x][0]) for x in specrange])

            spectrumflux.append([float(spectrumdata[x][1]) for x in specrange])
            if 'errorunit' in spectrum:
                spectrumerrs.append([float(spectrumdata[x][2]) for x in specrange])
                spectrumerrs[-1] = [x if is_number(x) and not isnan(float(x)) else 0. for x in spectrumerrs[-1]]

            if 'timeunit' not in spectrum or 'time' not in spectrum:
                hasepoch = False

            mjdmax = ''
            if 'timeunit' in spectrum and spectrum['timeunit'] == 'MJD' and 'redshift' in catalog[entry]:
                if 'maxdate' in catalog[entry]:
                    mjdmax = astrotime(catalog[entry]['maxdate'][0]['value'].replace('/', '-')).mjd
                if mjdmax:
                    hasmjdmax = True
                    mjdmax = (float(spectrum['time']) - mjdmax) / (1.0 + float(catalog[entry]['redshift'][0]['value']))
                    spectrummjdmax.append(mjdmax)

        nspec = len(catalog[entry]['spectra'])
        
        spectrumscaled = deepcopy(spectrumflux)
        for f, flux in enumerate(spectrumscaled):
            mean = numpy.std(flux)
            spectrumscaled[f] = [x/mean for x in flux]

        y_height = 0.
        y_offsets = [0. for x in range(nspec)]
        for i in reversed(range(nspec)):
            y_offsets[i] = y_height
            if (i-1 >= 0 and 'time' in catalog[entry]['spectra'][i] and 'time' in catalog[entry]['spectra'][i-1]
                and catalog[entry]['spectra'][i]['time'] == catalog[entry]['spectra'][i-1]['time']):
                    ydiff = 0
            else:
                ydiff = 0.8*(max(spectrumscaled[i]) - min(spectrumscaled[i]))
            spectrumscaled[i] = [j + y_height for j in spectrumscaled[i]]
            y_height += ydiff

        maxsw = max(map(max, spectrumwave))
        minsw = min(map(min, spectrumwave))
        maxfl = max(map(max, spectrumscaled))
        minfl = min(map(min, spectrumscaled))
        maxfldiff = max(map(operator.sub, list(map(max, spectrumscaled)), list(map(min, spectrumscaled))))
        x_buffer = 0.0 #0.1*(maxsw - minsw)
        x_range = [-x_buffer + minsw, x_buffer + maxsw]
        y_buffer = 0.1*maxfldiff
        y_range = [-y_buffer + minfl, y_buffer + maxfl]

        for f, flux in enumerate(spectrumscaled):
            spectrumscaled[f] = [x - y_offsets[f] for x in flux]

        tt2 = [ ("Source ID(s)", "@src") ]
        if 'redshift' in catalog[entry]:
            tt2 += [ ("λ (rest)", "@xrest{1.1} Å") ]
        tt2 += [
                ("λ (obs)", "@x{1.1} Å"),
                ("Flux", "@yorig"),
                ("Flux unit", "@fluxunit")
               ]

        if hasepoch:
            tt2 += [ ("Epoch (" + spectrum['timeunit'] + ")", "@epoch{1.11}") ]

        if hasmjdmax:
            tt2 += [ ("Rest days to max", "@mjdmax{1.11}") ]

        hover2 = HoverTool(tooltips = tt2)

        p2 = Figure(title='Spectra for ' + eventname, x_axis_label=label_format('Observed Wavelength (Å)'),
            y_axis_label=label_format('Flux (scaled)' + (' + offset'
            if (nspec > 1) else '')), x_range = x_range, tools = tools, #responsive = True,
            plot_width = 485, plot_height = 485, y_range = y_range, title_text_font_size='16pt',
            title_text_font = 'futura')
        p2.xaxis.axis_label_text_font = 'futura'
        p2.yaxis.axis_label_text_font = 'futura'
        p2.xaxis.major_label_text_font = 'futura'
        p2.yaxis.major_label_text_font = 'futura'
        p2.xaxis.axis_label_text_font_size = '12pt'
        p2.yaxis.axis_label_text_font_size = '12pt'
        p2.xaxis.major_label_text_font_size = '8pt'
        p2.yaxis.major_label_text_font_size = '8pt'
        p2.add_tools(hover2)

        sources = []
        for i in range(len(spectrumwave)):
            sl = len(spectrumscaled[i])

            data = dict(
                x0 = spectrumwave[i],
                y0 = spectrumscaled[i],
                yorig = spectrumflux[i],
                fluxunit = [label_format(catalog[entry]['spectra'][i]['fluxunit'])]*sl,
                x = spectrumwave[i],
                y = [y_offsets[i] + j for j in spectrumscaled[i]],
                src = [catalog[entry]['spectra'][i]['source']]*sl
            )
            if 'redshift' in catalog[entry]:
                data['xrest'] = [x/(1.0 + z) for x in spectrumwave[i]]
            if hasepoch:
                data['epoch'] = [catalog[entry]['spectra'][i]['time'] for j in spectrumscaled[i]]
            if hasmjdmax:
                data['mjdmax'] = [spectrummjdmax[i] for j in spectrumscaled[i]]
            sources.append(ColumnDataSource(data))
            p2.line('x', 'y', source=sources[i], color=mycolors[i % len(mycolors)], line_width=2, line_join='round')

        if 'redshift' in catalog[entry]:
            minredw = minsw/(1.0 + z)
            maxredw = maxsw/(1.0 + z)
            p2.extra_x_ranges = {"other wavelength": Range1d(start=minredw, end=maxredw)}
            p2.add_layout(LinearAxis(axis_label ="Restframe Wavelength (Å)",
                x_range_name="other wavelength", axis_label_text_font_size = '12pt',
                axis_label_text_font = 'futura',
                major_label_text_font_size = '8pt', major_label_text_font = 'futura'), 'above')

        sdicts = dict(zip(['s'+str(x) for x in range(len(sources))], sources))
        callback = CustomJS(args=sdicts, code="""
            var yoffs = [""" + ','.join([str(x) for x in y_offsets]) + """];
            for (s = 0; s < """ + str(len(sources)) + """; s++) {
                var data = eval('s'+s).get('data');
                var redshift = """ + str(z if 'redshift' in catalog[entry] else 0.) + """;
                if (!('binsize' in data)) {
                    data['binsize'] = 1.0
                }
                if (!('spacing' in data)) {
                    data['spacing'] = 1.0
                }
                if (cb_obj.get('title') == 'Spacing') {
                    data['spacing'] = cb_obj.get('value');
                } else {
                    data['binsize'] = cb_obj.get('value');
                }
                var f = data['binsize']
                var space = data['spacing']
                var x0 = data['x0'];
                var y0 = data['y0'];
                var dx0 = x0[1] - x0[0];
                var yoff = space*yoffs[s];
                data['x'] = [x0[0] - 0.5*Math.max(0., f - dx0)];
                data['xrest'] = [(x0[0] - 0.5*Math.max(0., f - dx0))/(1.0 + redshift)];
                data['y'] = [y0[0] + yoff];
                var xaccum = 0.;
                var yaccum = 0.;
                for (i = 0; i < x0.length; i++) {
                    var dx;
                    if (i == 0) {
                        dx = x0[i+1] - x0[i];
                    } else {
                        dx = x0[i] - x0[i-1];
                    }
                    xaccum += dx;
                    yaccum += y0[i]*dx;
                    if (xaccum >= f) {
                        data['x'].push(data['x'][data['x'].length-1] + xaccum);
                        data['xrest'].push(data['x'][data['x'].length-1]/(1.0 + redshift));
                        data['y'].push(yaccum/xaccum + yoff);
                        xaccum = 0.;
                        yaccum = 0.;
                    }
                }
                eval('s'+s).trigger('change');
            }
        """)

        binslider = Slider(start=0, end=20, value=1, step=0.5, title=label_format("Bin size (Angstrom)"), callback=callback)
        spacingslider = Slider(start=0, end=2, value=1, step=0.02, title=label_format("Spacing"), callback=callback)

    hasimage = False
    skyhtml = ''
    if 'ra' in catalog[entry] and 'dec' in catalog[entry] and args.collecthosts:
        snra = catalog[entry]['ra'][0]['value']
        sndec = catalog[entry]['dec'][0]['value']
        c = coord(ra=snra, dec=sndec, unit=(un.hourangle, un.deg))

        if 'lumdist' in catalog[entry] and float(catalog[entry]['lumdist'][0]['value']) > 0.:
            if 'host' in catalog[entry] and catalog[entry]['host'][0]['value'] == 'Milky Way':
                sdssimagescale = max(0.05,0.04125/float(catalog[entry]['lumdist'][0]['value']))
            else:
                sdssimagescale = max(0.05,20.6265/float(catalog[entry]['lumdist'][0]['value']))
        else:
            if 'host' in catalog[entry] and catalog[entry]['host'][0]['value'] == 'Milky Way':
                sdssimagescale = 0.0006
            else:
                sdssimagescale = 0.3
        dssimagescale = 0.13889*sdssimagescale
        #At the moment, no way to check if host is in SDSS footprint without comparing to empty image, which is only possible at fixed angular resolution.
        sdssimagescale = 0.3

        imgsrc = ''
        hasimage = True
        if eventname in hostimgdict:
            imgsrc = hostimgdict[eventname]
        else:
            try:
                response = urllib.request.urlopen('http://skyservice.pha.jhu.edu/DR12/ImgCutout/getjpeg.aspx?ra='
                    + str(c.ra.deg) + '&dec=' + str(c.dec.deg) + '&scale=' + sdssimagescale + '&width=500&height=500&opt=G')
            except:
                hasimage = False
            else:
                with open(outdir + fileeventname + '-host.jpg', 'wb') as f:
                    f.write(response.read())
                imgsrc = 'SDSS'

            if hasimage and filecmp.cmp(outdir + fileeventname + '-host.jpg', outdir + 'missing.jpg'):
                hasimage = False

            if not hasimage:
                hasimage = True
                url = ("http://skyview.gsfc.nasa.gov/current/cgi/runquery.pl?Position=" + str(urllib.parse.quote_plus(snra + " " + sndec)) +
                       "&coordinates=J2000&coordinates=&projection=Tan&pixels=500&size=" + str(dssimagescale) + "&float=on&scaling=Log&resolver=SIMBAD-NED" +
                       "&Sampler=_skip_&Deedger=_skip_&rotation=&Smooth=&lut=colortables%2Fb-w-linear.bin&PlotColor=&grid=_skip_&gridlabels=1" +
                       "&catalogurl=&CatalogIDs=on&RGB=1&survey=DSS2+IR&survey=DSS2+Red&survey=DSS2+Blue&IOSmooth=&contour=&contourSmooth=&ebins=null")

                response = urllib.request.urlopen(url)
                bandsoup = BeautifulSoup(response, "html5lib")
                images = bandsoup.findAll('img')
                imgname = ''
                for image in images:
                    if "Quicklook RGB image" in image.get('alt', ''):
                        imgname = image.get('src', '').split('/')[-1]

                if imgname:
                    try:
                        response = urllib.request.urlopen('http://skyview.gsfc.nasa.gov/tempspace/fits/' + imgname)
                    except:
                        hasimage = False
                    else:
                        with open(outdir + fileeventname + '-host.jpg', 'wb') as f:
                            f.write(response.read())
                        imgsrc = 'DSS'
                else:
                    hasimage = False

        if hasimage:
            if imgsrc == 'SDSS':
                hostimgs.append([eventname, 'SDSS'])
                skyhtml = ('<a href="http://skyserver.sdss.org/DR12/en/tools/chart/navi.aspx?opt=G&ra='
                    + str(c.ra.deg) + '&dec=' + str(c.dec.deg) + '&scale=0.15"><img style="margin:5px;" src="' + fileeventname + '-host.jpg" width=250></a>')
            elif imgsrc == 'DSS':
                hostimgs.append([eventname, 'DSS'])
                url = ("http://skyview.gsfc.nasa.gov/current/cgi/runquery.pl?Position=" + str(urllib.parse.quote_plus(snra + " " + sndec)) +
                       "&coordinates=J2000&coordinates=&projection=Tan&pixels=500&size=" + str(dssimagescale) + "float=on&scaling=Log&resolver=SIMBAD-NED" +
                       "&Sampler=_skip_&Deedger=_skip_&rotation=&Smooth=&lut=colortables%2Fb-w-linear.bin&PlotColor=&grid=_skip_&gridlabels=1" +
                       "&catalogurl=&CatalogIDs=on&RGB=1&survey=DSS2+IR&survey=DSS2+Red&survey=DSS2+Blue&IOSmooth=&contour=&contourSmooth=&ebins=null")
                skyhtml = ('<a href="' + url + '"><img style="margin:5px;" src="' + fileeventname + '-host.jpg" width=250></a>')
        else:
            hostimgs.append([eventname, 'None'])

    if dohtml and args.writehtml:
    #if (photoavail and spectraavail) and dohtml and args.writehtml:
        if photoavail and spectraavail:
            p = HBox(p1,VBox(p2,HBox(binslider,spacingslider)))
            #script, div = components(dict(p1=p1, p2=p2))#, binslider=binslider, spacingslider=spacingslider))
        elif photoavail:
            p = p1
            #script, div = components(dict(p1=p1))
        elif spectraavail:
            p = VBox(HBox(p2,VBox(binslider,spacingslider)), width=900)
            #script, div = components(dict(p2=p2, binslider=binslider, spacingslider=spacingslider))

        html = '<html><head><title>'+eventname+'</title>'
        if photoavail or spectraavail:
            html = file_html(p, CDN, eventname)
            #html = html + '''<link href="https://cdn.pydata.org/bokeh/release/bokeh-0.11.0.min.css" rel="stylesheet" type="text/css">
            #    <script src="https://cdn.pydata.org/bokeh/release/bokeh-0.11.0.min.js"></script>''' + script + '</head><body>'
        else:
            html = '<html><title></title><body></body></html>'

        #if photoavail and spectraavail:
        #    html = html + div['p1'] + div['p2']# + div['binslider'] + div['spacingslider']
        #elif photoavail:
        #    html = html + div['p1']
        #elif spectraavail:
        #    html = html + div['p2'] + div['binslider'] + div['spacingslider']

        #html = html + '</body></html>'

        html = html.replace('<body>',
            '''<body class='event-body'><div style="padding-bottom:8px;"><strong>Disclaimer:</strong> All data collected by the OSC was originally generated by others, if you intend to use this data in a publication, we ask that you please cite the linked sources and/or contact the sources of the data directly. Data sources are revealed by hovering over the data with your cursor.</div>''')
        html = re.sub(r'(\<\/title\>)', r'''\1\n
            <base target="_parent" />\n
            <link rel="stylesheet" href="event.css" type="text/css">\n
            <script type="text/javascript">\n
                if(top==self)\n
                this.location="''' + eventname + '''"\n
            </script>'''
            , html)

        repfolder = get_rep_folder(catalog[entry])
        html = re.sub(r'(\<\/body\>)', '<div style="width:100%; text-align:center;">' + r'<a class="event-download" href="' +
            linkdir + fileeventname + r'.json" download>' + r'Click to download all data for ' + eventname + ' in JSON format' +
            r'</a></div>\n\1', html)

        newhtml = r'<div class="event-tab-div"><h3 class="event-tab-title">Event metadata</h3><table class="event-table"><tr><th width=100px class="event-cell">Quantity</th><th class="event-cell">Value<sup>sources</sup></th></tr>\n'
        for key in columnkey:
            if key in catalog[entry] and key not in eventignorekey and len(catalog[entry][key]) > 0:
                keyhtml = ''
                if isinstance(catalog[entry][key], str):
                    keyhtml = keyhtml + re.sub('<[^<]+?>', '', catalog[entry][key])
                else:
                    for r, row in enumerate(catalog[entry][key]):
                        if 'value' in row and 'source' in row:
                            sources = [str(x) for x in sorted([x.strip() for x in row['source'].split(',')], key=lambda x: float(x) if is_number(x) else float("inf"))]
                            sourcehtml = ''
                            for s, source in enumerate(sources):
                                if source == 'D':
                                    sourcehtml = sourcehtml + (',' if s > 0 else '') + source
                                else:
                                    sourcehtml = sourcehtml + (',' if s > 0 else '') + r'<a href="#source' + source + r'">' + source + r'</a>'
                            keyhtml = keyhtml + (r'<br>' if r > 0 else '') + row['value']
                            if ((key == 'maxdate' or key == 'maxabsmag' or key == 'maxappmag') and 'maxband' in catalog[entry]
                                and catalog[entry]['maxband']):
                                keyhtml = keyhtml + r' [' + catalog[entry]['maxband'][0]['value'] + ']'
                            keyhtml = keyhtml + r'<sup>' + sourcehtml + r'</sup>'
                        elif isinstance(row, str):
                            keyhtml = keyhtml + (r'<br>' if r > 0 else '') + row.strip()
                if keyhtml:
                    newhtml = (newhtml + r'<tr><td class="event-cell">' + eventpageheader[key] +
                        r'</td><td width=250px class="event-cell">' + keyhtml)

                newhtml = newhtml + r'</td></tr>\n'
        newhtml = newhtml + r'</table><em>D = Derived value</em></div>\n\1'
        html = re.sub(r'(\<\/body\>)', newhtml, html)

        if 'sources' in catalog[entry] and len(catalog[entry]['sources']):
            newhtml = r'<div class="event-tab-div"><h3 class="event-tab-title">Sources of data</h3><table class="event-table"><tr><th width=30px class="event-cell">ID</th><th class="event-cell">Source</th></tr>\n'
            for source in catalog[entry]['sources']:
                hasurlnobib = ('url' in source and 'bibcode' not in source)
                newhtml = (newhtml + r'<tr><td class="event-cell" id="source' + source['alias'] + '">' + source['alias'] +
                    r'</td><td width=250px class="event-cell">' + (('<a href="' + source['url'] + '">') if hasurlnobib else '') +
                    source['name'].encode('ascii', 'xmlcharrefreplace').decode("utf-8") +
                    (r'</a>' if hasurlnobib else '') +
                    ((r'<br>\n' + (('<a href="' + source['url'] + '">') if 'url' in source else '') + source['bibcode'] +
                    (r'</a>' if 'url' in source else '')) if 'bibcode' in source and source['name'] != source['bibcode'] else '') +
                    r'</td></tr>\n')
            newhtml = newhtml + r'</table><em>Sources are presented in order of importation, not in order of importance.</em></div>\n'

            if hasimage:
                newhtml = newhtml + '<div class="event-host-div"><h3 class="event-host-title">Host Image</h3>' + skyhtml
                newhtml = newhtml + r'</table><em>Host images are taken from SDSS if available; if not, DSS is used.</em></div>\n'

        newhtml = newhtml + r'\n\1'

        html = re.sub(r'(\<\/body\>)', newhtml, html)

        with open(outdir + fileeventname + ".html", "w") as fff:
            fff.write(html)

    # Necessary to clear Bokeh state
    reset_output()

    #if spectraavail and dohtml:
    #    sys.exit()

    #if fcnt > 100:
    #    sys.exit()

    # Save this stuff because next line will delete it.
    if args.writecatalog:
        # Construct array for Bishop's webpage
        # Things David wants in this file: names (aliases), max mag, max mag date (gregorian), type, redshift (helio), redshift (host), r.a., dec., # obs., link
        snepages.append([entry, ",".join(catalog[entry]['aliases']), get_first_value(entry, 'maxappmag'), get_first_value(entry, 'maxdate'),
            get_first_value(entry, 'claimedtype'), get_first_value(entry, 'redshift'), get_first_kind(entry, 'redshift'),
            get_first_value(entry, 'ra'), get_first_value(entry, 'dec'), catalog[entry]['numphoto'], 'https://sne.space/' + plotlink])

        if 'sources' in catalog[entry]:
            lsourcedict = {}
            for sourcerow in catalog[entry]['sources']:
                strippedname = re.sub('<[^<]+?>', '', sourcerow['name'].encode('ascii','xmlcharrefreplace').decode("utf-8"))
                alias = sourcerow['alias']
                if 'bibcode' in sourcerow and 'secondary' not in sourcerow:
                    lsourcedict[alias] = {'bibcode':sourcerow['bibcode'], 'count':0}
                if strippedname in sourcedict:
                    sourcedict[strippedname] += 1
                else:
                    sourcedict[strippedname] = 1

            for key in catalog[entry].keys():
                if isinstance(catalog[entry][key], list):
                    for row in catalog[entry][key]:
                        if 'source' in row:
                            for lsource in lsourcedict:
                                if lsource in row['source'].split(','):
                                    if key == 'spectra':
                                        lsourcedict[lsource]['count'] += 10
                                    else:
                                        lsourcedict[lsource]['count'] += 1

            ssources = sorted(list(lsourcedict.values()), key=lambda x: x['count'], reverse=True)
            if ssources:
                seemorelink = ''
                catalog[entry]['references'] = ','.join([y['bibcode'] for y in ssources[:5]])

        lcspye.append(catalog[entry]['numphoto'] >= 5 and catalog[entry]['numspectra'] >  0)
        lconly.append(catalog[entry]['numphoto'] >= 5 and catalog[entry]['numspectra'] == 0)
        sponly.append(catalog[entry]['numphoto'] <  5 and catalog[entry]['numspectra'] >  0)
        lcspno.append(catalog[entry]['numphoto'] <  5 and catalog[entry]['numspectra'] == 0)

        hasalc.append(catalog[entry]['numphoto'] >= 5)
        hasasp.append(catalog[entry]['numspectra'] >  0)

        totalphoto += catalog[entry]['numphoto']
        totalspectra += catalog[entry]['numspectra']

        # Delete unneeded data from catalog, add blank entries when data missing.
        catalogcopy[entry] = OrderedDict()
        for col in columnkey:
            if col in catalog[entry]:
                catalogcopy[entry][col] = catalog[entry][col]
            else:
                catalogcopy[entry][col] = None

        del catalog[entry]

    if args.test and spectraavail and photoavail:
        break

# Write it all out at the end
if args.writecatalog and not args.eventlist:
    catalog = catalogcopy

    #Write the MD5 checksums
    jsonstring = json.dumps(md5s, separators=(',',':'))
    with open(outdir + 'md5s.json' + testsuffix, 'w') as f:
        f.write(jsonstring)

    #Write the host image info
    jsonstring = json.dumps(hostimgs, separators=(',',':'))
    with open(outdir + 'hostimgs.json' + testsuffix, 'w') as f:
        f.write(jsonstring)

    # Things David wants in this file: names (aliases), max mag, max mag date (gregorian), type, redshift, r.a., dec., # obs., link
    with open(outdir + 'snepages.csv' + testsuffix, 'w') as f:
        csvout = csv.writer(f, quotechar='"', quoting=csv.QUOTE_ALL)
        for row in snepages:
            csvout.writerow(row)

    # Make a few small files for generating charts
    with open(outdir + 'sources.csv' + testsuffix, 'w') as f:
        sortedsources = sorted(list(sourcedict.items()), key=operator.itemgetter(1), reverse=True)
        csvout = csv.writer(f)
        csvout.writerow(['Source','Number'])
        for source in sortedsources:
            csvout.writerow(source)

    with open(outdir + 'pie.csv' + testsuffix, 'w') as f:
        csvout = csv.writer(f)
        csvout.writerow(['Category','Number'])
        csvout.writerow(['Has light curve and spectra', sum(lcspye)])
        csvout.writerow(['Has light curve only', sum(lconly)])
        csvout.writerow(['Has spectra only', sum(sponly)])
        csvout.writerow(['No light curve or spectra', sum(lcspno)])

    with open(outdir + 'hasphoto.html' + testsuffix, 'w') as f:
        f.write("{:,}".format(sum(hasalc)))
    with open(outdir + 'hasspectra.html' + testsuffix, 'w') as f:
        f.write("{:,}".format(sum(hasasp)))
    with open(outdir + 'snecount.html' + testsuffix, 'w') as f:
        f.write("{:,}".format(len(catalog)))
    with open(outdir + 'photocount.html' + testsuffix, 'w') as f:
        f.write("{:,}".format(totalphoto))
    with open(outdir + 'spectracount.html' + testsuffix, 'w') as f:
        f.write("{:,}".format(totalspectra))

    ctypedict = dict()
    for entry in catalog:
        cleanedtype = ''
        if 'claimedtype' in catalog[entry] and catalog[entry]['claimedtype']:
            maxsources = 0
            for ct in catalog[entry]['claimedtype']:
                sourcecount = len(ct['source'].split(','))
                if sourcecount > maxsources:
                    maxsources = sourcecount
                    cleanedtype = ct['value'].strip('?* ')
        if not cleanedtype:
            cleanedtype = 'Unknown'
        if cleanedtype in ctypedict:
            ctypedict[cleanedtype] += 1
        else:
            ctypedict[cleanedtype] = 1
    sortedctypes = sorted(list(ctypedict.items()), key=operator.itemgetter(1), reverse=True)
    with open(outdir + 'types.csv' + testsuffix, 'w') as f:
        csvout = csv.writer(f)
        csvout.writerow(['Type','Number'])
        for ctype in sortedctypes:
            csvout.writerow(ctype)

    with open('../../sitemap.xml', 'w') as f:
        sitemapxml = sitemaptemplate
        sitemaplocs = ''
        for key in catalog.keys():
            sitemaplocs = sitemaplocs + "  <url>\n    <loc>https://sne.space/sne/" + key + "</loc>\n  </url>\n"
        sitemapxml = sitemapxml.replace('{0}', sitemaplocs)
        f.write(sitemapxml)

    # Ping Google to let them know sitemap has been updated
    #response = urllib.request.urlopen(googlepingurl)

    # Convert to array since that's what datatables expects
    catalog = list(catalog.values())

    jsonstring = json.dumps(catalog, separators=(',',':'))
    with open(outdir + 'catalog.min.json' + testsuffix, 'w') as f:
        f.write(jsonstring)

    jsonstring = json.dumps(catalog, indent='\t', separators=(',',':'))
    with open(outdir + 'catalog.json' + testsuffix, 'w') as f:
        f.write(jsonstring)

    with open(outdir + 'catalog.html' + testsuffix, 'w') as f:
        f.write('<table id="example" class="display" cellspacing="0" width="100%">\n')
        f.write('\t<thead>\n')
        f.write('\t\t<tr>\n')
        for h in header:
            f.write('\t\t\t<th class="' + h + '" title="' + titles[h] + '">' + header[h] + '</th>\n')
        f.write('\t\t</tr>\n')
        f.write('\t</thead>\n')
        f.write('\t<tfoot>\n')
        f.write('\t\t<tr>\n')
        for h in header:
            f.write('\t\t\t<th class="' + h + '" title="' + titles[h] + '">' + header[h] + '</th>\n')
        f.write('\t\t</tr>\n')
        f.write('\t</thead>\n')
        f.write('</table>\n')

    with open(outdir + 'catalog.min.json', 'rb') as f_in, gzip.open(outdir + 'catalog.min.json.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
