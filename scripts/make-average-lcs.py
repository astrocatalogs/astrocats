#!/usr/local/bin/python3.5

import json
import re
import sys
import glob
from photometry import bandcolorf, bandaliasf
from bokeh.plotting import Figure, show, save, reset_output
from bokeh.models import (HoverTool, CustomJS, Slider, ColumnDataSource,
                          HBox, VBox, Range1d, LinearAxis, DatetimeAxis)
from bokeh.resources import CDN, INLINE
from bokeh.embed import file_html, components
from random import randint, uniform
from math import log10, floor
from collections import OrderedDict
from astropy.time import Time as astrotime

phototime = []
phototimelowererrs = []
phototimeuppererrs = []
photoAB = []
photoABerrs = []
photoband = []
photoinstru = []
photoevent = []
phototype = []

tools = "pan,wheel_zoom,box_zoom,save,crosshair,reset,resize"

outdir = "../"

averagetype = "II"

def get_sig_digits(x):
    return len((''.join(x.split('.'))).strip('0'))

def round_sig(x, sig=4):
    if x == 0.0:
        return 0.0
    return round(x, sig-int(floor(log10(abs(x))))-1)

with open('rep-folders.txt', 'r') as f:
    repfolders = f.read().splitlines()

files = []
for rep in repfolders:
    files += glob.glob('../' + rep + "/*.json")

for fcnt, eventfile in enumerate(sorted(files, key=lambda s: s.lower())):
    #if fcnt > 2000:
    #    break

    with open(eventfile, 'r') as f:
        filetext = f.read()

    thisevent = json.loads(filetext, object_pairs_hook=OrderedDict)
    thisevent = thisevent[list(thisevent.keys())[0]]

    if ('photometry' not in thisevent or 'maxdate' not in thisevent or 'maxabsmag' not in thisevent or
        'maxappmag' not in thisevent or 'claimedtype' not in thisevent or
        len(thisevent['maxdate'][0]['value'].split('/')) < 3 or 'discoverdate' not in thisevent):
        continue

    foundtype = False
    for ct in thisevent['claimedtype']:
        if ct['value'] == averagetype:
            foundtype = True
            break

    if not foundtype:
        continue

    maxdate = astrotime(thisevent['maxdate'][0]['value'].replace('/', '-')).mjd
    discoverdate = astrotime(thisevent['discoverdate'][0]['value'].replace('/', '-')).mjd

    if maxdate == discoverdate:
        continue

    distmod = float(thisevent['maxappmag'][0]['value']) - float(thisevent['maxabsmag'][0]['value'])

    print(thisevent['name'])

    prange = list(range(len(thisevent['photometry']))) if 'photometry' in thisevent else []

    #phototime += [float(x['time']) - maxdate for x in thisevent['photometry'] if 'magnitude' in x]
    phototime += [float(x['time'][:-1] + str(randint(0,9)) if x['time'][-1] != '.' else x['time'] + '.' + str(randint(0,9))) - maxdate
                  for x in thisevent['photometry'] if 'magnitude' in x]
    #phototime += [(float(x['time']) + uniform(-1,1)*10.**(log10(abs(float(x['time']) - round_sig(float(x['time']), get_sig_digits(x['time'])-1)))) - maxdate)
    #              for x in thisevent['photometry'] if 'magnitude' in x]
    phototimelowererrs += [float(x['e_lower_time']) if ('e_lower_time' in x and 'e_upper_time' in x)
        else (float(x['e_time']) if 'e_time' in x else 0.) for x in thisevent['photometry'] if 'magnitude' in x]
    phototimeuppererrs += [float(x['e_upper_time']) if ('e_lower_time' in x and 'e_upper_time' in x) in x
        else (float(x['e_time']) if 'e_time' in x else 0.) for x in thisevent['photometry'] if 'magnitude' in x]
    photoAB += [float(x['magnitude'] + str(randint(0,9)) if '.' in x['magnitude'] else x['magnitude'] + '.' + str(randint(0,9))) - distmod for x in thisevent['photometry'] if 'magnitude' in x]
    photoABerrs += [(float(x['e_magnitude']) if 'e_magnitude' in x else 0.) for x in thisevent['photometry'] if 'magnitude' in x]
    photoband += [(x['band'] if 'band' in x else '') for x in thisevent['photometry'] if 'magnitude' in x]
    photoinstru += [(x['instrument'] if 'instrument' in x else '') for x in thisevent['photometry'] if 'magnitude' in x]
    photoevent += [thisevent['name'] for x in prange]
    phototype += [(x['upperlimit'] if 'upperlimit' in x else False) for x in thisevent['photometry'] if 'magnitude' in x]

bandset = set(photoband)
bandset = [i for (j, i) in sorted(list(zip(list(map(bandaliasf, bandset)), bandset)))]

x_buffer = 0.1*(max(phototime) - min(phototime)) if len(phototime) > 1 else 1.0

tt = [  
        ("Event", "@src"),
        ("Epoch (MJD)", "@x{1.11}"),
        ("Absolute Magnitude", "@y{1.111}")
     ]
if len(list(filter(None, photoABerrs))):
    tt += [("Error", "@err{1.111}")]
if len(list(filter(None, photoband))):
    tt += [("Band", "@desc")]
if len(list(filter(None, photoinstru))):
    tt += [("Instrument", "@instr")]
hover = HoverTool(tooltips = tt)

min_x_range = -x_buffer + min([x - y for x, y in list(zip(phototime, phototimeuppererrs))])
max_x_range = x_buffer + max([x + y for x, y in list(zip(phototime, phototimelowererrs))])

p1 = Figure(title='Average Photometry for Type ' + averagetype + ' SNe', x_axis_label='Time (MJD)',
    y_axis_label='Absolute Magnitude', tools = tools, plot_width = 1000, plot_height = 1000, #responsive = True,
    x_range = (min_x_range, max_x_range),
    y_range = (0.5 + max([x + y for x, y in list(zip(photoAB, photoABerrs))]),
               -0.5 + min([x - y for x, y in list(zip(photoAB, photoABerrs))])),
    title_text_font_size='20pt', webgl = True)
p1.xaxis.axis_label_text_font_size = '16pt'
p1.yaxis.axis_label_text_font_size = '16pt'
p1.xaxis.major_label_text_font_size = '12pt'
p1.yaxis.major_label_text_font_size = '12pt'

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

    source = ColumnDataSource(
        data = dict(
            x = [phototime[i] for i in indne],
            y = [photoAB[i] for i in indne],
            err = [photoABerrs[i] for i in indne],
            desc = [photoband[i] for i in indne],
            instr = [photoinstru[i] for i in indne],
            src = [photoevent[i] for i in indne]
        )
    )
    p1.circle('x', 'y', source = source, color=bandcolorf(band), legend='', size=2, line_alpha=0.75, fill_alpha=0.75)

    source = ColumnDataSource(
        data = dict(
            x = [phototime[i] for i in indye],
            y = [photoAB[i] for i in indye],
            err = [photoABerrs[i] for i in indye],
            desc = [photoband[i] for i in indye],
            instr = [photoinstru[i] for i in indye],
            src = [photoevent[i] for i in indye]
        )
    )
    #p1.multi_line([err_xs[x] for x in indye], [[ys[x], ys[x]] for x in indye], color=bandcolorf(band))
    #p1.multi_line([[xs[x], xs[x]] for x in indye], [err_ys[x] for x in indye], color=bandcolorf(band))
    p1.circle('x', 'y', source = source, color=bandcolorf(band), legend=bandname, size=2, line_alpha=0.75, fill_alpha=0.75)

    #upplimlegend = bandname if len(indye) == 0 and len(indne) == 0 else ''

    #indt = [i for i, j in enumerate(phototype) if j]
    #ind = set(indb).intersection(indt)
    #p1.inverted_triangle([phototime[x] for x in ind], [photoAB[x] for x in ind],
    #    color=bandcolorf(band), legend=upplimlegend, size=7)

p1.legend.label_text_font_size = '8pt'
p1.legend.label_width = 20
p1.legend.label_height = 14
p1.legend.glyph_height = 14

html = file_html(p1, CDN, 'Average ' + averagetype)

with open(outdir + "average-" + averagetype.lower().replace(' ', '_').replace('/', '-') + ".html", "w") as f:
    f.write(html)
