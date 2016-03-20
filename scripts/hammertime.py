#!/usr/local/bin/python3.5

import json
import sys
import glob
from bokeh.plotting import Figure, show, save, reset_output
from bokeh.models import (HoverTool, CustomJS, Slider, ColumnDataSource,
                          HBox, VBox, Range1d, LinearAxis, DatetimeAxis)
from bokeh.resources import CDN
from bokeh.embed import file_html, components
from collections import OrderedDict
from astropy.coordinates import SkyCoord as coord
from astropy import units as un
from palettable import cubehelix
from math import *
from random import shuffle, seed, randint, uniform

tools = "pan,wheel_zoom,box_zoom,save,crosshair,reset,resize"

outdir = "../"

snhxs = []
snhys = []
snras = []
sndecs = []
sntypes = []
snnames = []

seed(12483)
colors = (cubehelix.cubehelix1_16.hex_colors[2:13] +
          cubehelix.cubehelix2_16.hex_colors[2:13] +
          cubehelix.cubehelix3_16.hex_colors[2:13] +
          cubehelix.jim_special_16.hex_colors[2:13] +
          cubehelix.purple_16.hex_colors[2:13] +
          cubehelix.perceptual_rainbow_16.hex_colors)
shuffle(colors)

with open('rep-folders.txt', 'r') as f:
    repfolders = f.read().splitlines()

files = []
for rep in repfolders:
    files += glob.glob('../' + rep + "/*.json")

for fcnt, eventfile in enumerate(sorted(files, key=lambda s: s.lower())):
    #if fcnt > 20:
    #    break

    with open(eventfile, 'r') as f:
        filetext = f.read()

    thisevent = json.loads(filetext, object_pairs_hook=OrderedDict)
    thisevent = thisevent[list(thisevent.keys())[0]]

    if 'ra' in thisevent and 'dec' in thisevent:
        if 'claimedtype' in thisevent and thisevent['claimedtype']:
            thistype = thisevent['claimedtype'][0]['value'].replace('?', '').replace('*', '')
            if thistype in ('', 'QSO', 'AGN', 'Nova', 'Galaxy', 'CV'):
                continue
            elif thistype in ('Other', 'not Ia', 'SN', 'unconf', 'Radio', 'CC', 'CCSN', 'Candidate'):
                sntypes.append('Unknown')
            else:
                sntypes.append(thistype)
        else:
            sntypes.append('Unknown')

        print(thisevent['name'])
        c = coord(ra=thisevent['ra'][0]['value'], dec=thisevent['dec'][0]['value'], unit=(un.hourangle, un.deg))
        snnames.append(thisevent['name'])
        rarad = c.ra.radian - pi
        decrad = c.dec.radian
        snhx = 2.0**1.5*cos(decrad)*sin(rarad/2.0)/sqrt(1.0 + cos(decrad)*cos(rarad/2.0))
        snhy = sqrt(2.0)*sin(decrad)/sqrt(1.0 + cos(decrad)*cos(rarad/2.0))
        snras.append(c.ra.deg)
        sndecs.append(c.dec.deg)
        snhxs.append(snhx)
        snhys.append(snhy)

rangepts = 100
raseps = 24
decseps = 18
rarange = [-pi + i*2.0*pi/rangepts for i in range(0, rangepts+1)]
decrange = [-pi/2.0 + i*pi/rangepts for i in range(0, rangepts+1)]

ragrid = [-pi + i*2.0*pi/raseps for i in range(0, raseps+1)]
decgrid = [-pi/2.0 + i*pi/decseps for i in range(0, decseps+1)]

tt = [  
        ("Event", "@event"),
        ("R.A. (deg)", "@ra{1.111}"),
        ("Dec. (deg)", "@dec{1.111}"),
        ("Claimed type", "@claimedtype")
     ]
hover = HoverTool(tooltips = tt)

p1 = Figure(title='Supernova Positions', x_axis_label='Right Ascension (deg)',
    y_axis_label='Declination (deg)', tools = tools, plot_width = 1200, plot_height = 600, #responsive = True,
    x_range = (-1.2*(2.0**1.5), 1.2*2.0**1.5), y_range = (-1.2*sqrt(2.0), 1.2*sqrt(2.0)),
    title_text_font_size='20pt', webgl = True)
p1.axis.visible = None
p1.outline_line_color = None
p1.xgrid.grid_line_color = None
p1.ygrid.grid_line_color = None

raxs = []
rays = []
for rg in ragrid:
    raxs.append([2.0**1.5*cos(x)*sin(rg/2.0)/sqrt(1.0 + cos(x)*cos(rg/2.0)) for x in decrange])
    rays.append([sqrt(2.0)*sin(x)/sqrt(1.0 + cos(x)*cos(rg/2.0)) for x in decrange])

decxs = []
decys = []
for dg in decgrid:
    decxs.append([2.0**1.5*cos(dg)*sin(x/2.0)/sqrt(1.0 + cos(dg)*cos(x/2.0)) for x in rarange])
    decys.append([sqrt(2.0)*sin(dg)/sqrt(1.0 + cos(dg)*cos(x/2.0)) for x in rarange])

p1.add_tools(hover)
p1.multi_line(raxs, rays, color='#bbbbbb')
p1.multi_line(decxs, decys, color='#bbbbbb')
xs = []
ys = []

claimedtypes = sorted(list(set(sntypes)))

print(claimedtypes)

for ci, ct in enumerate(claimedtypes):
    ind = [i for i, t in enumerate(sntypes) if t == ct]

    source = ColumnDataSource(
        data = dict(
            x = [snhxs[i] for i in ind],
            y = [snhys[i] for i in ind],
            ra = [snras[i] for i in ind],
            dec = [sndecs[i] for i in ind],
            event = [snnames[i] for i in ind],
            claimedtype = [sntypes[i] for i in ind]
        )
    )
    if ct == 'Unknown':
        tcolor = 'black'
        falpha = 0.0
    else:
        tcolor = colors[ci]
        falpha = 1.0
    p1.circle('x', 'y', source = source, color=tcolor, fill_alpha=falpha, legend=ct, size=3)

p1.legend.label_text_font_size = '7pt'
p1.legend.label_width = 20
p1.legend.label_height = 10
p1.legend.glyph_height = 10

html = file_html(p1, CDN, 'Supernova locations')

with open(outdir + "sne-locations.html", "w") as f:
    f.write(html)

