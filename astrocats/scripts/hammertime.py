#!/usr/local/bin/python3.5

import argparse
import json
import warnings
import numpy as np
from collections import OrderedDict
from math import cos, pi, sin, sqrt
from random import seed, shuffle

from astropy import units as un
from astropy.coordinates import SkyCoord as coord
from bokeh.embed import file_html
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import Figure
from bokeh.resources import CDN
from palettable import cubehelix

from astrocats.catalog.utils import tprint, tq
from astrocats.scripts.events import get_event_text
from astrocats.scripts.repos import repo_file_list, get_rep_folders

parser = argparse.ArgumentParser(
    description='Generate a sky location map AstroCats data.'
)
parser.add_argument(
    '--catalog',
    '-c',
    dest='catalog',
    help='Select which catalog to generate',
    default='sne',
    type=str)
args = parser.parse_args()

if args.catalog == 'tde':
    moduledir = 'tidaldisruptions'
    modulename = 'tde'
    moduleurl = 'tde.space'
    moduletitle = 'TDE'
elif args.catalog == 'sne':
    moduledir = 'supernovae'
    modulename = 'sne'
    moduleurl = 'sne.space'
    moduletitle = 'Supernova'
else:
    raise ValueError('Unknown catalog!')

tools = "pan,wheel_zoom,box_zoom,save,crosshair,reset,resize"

outdir = "astrocats/" + moduledir + "/output/html/"

evhxs = []
evhys = []
evras = []
evdecs = []
evtypes = []
evnames = []

seed(12483)
colors = (cubehelix.cubehelix1_16.hex_colors[2:13] +
          cubehelix.cubehelix2_16.hex_colors[2:13] +
          cubehelix.cubehelix3_16.hex_colors[2:13] +
          cubehelix.jim_special_16.hex_colors[2:13] +
          cubehelix.purple_16.hex_colors[2:13] +
          cubehelix.purple_16.hex_colors[2:13] +
          cubehelix.purple_16.hex_colors[2:13] +
          cubehelix.purple_16.hex_colors[2:13] +
          cubehelix.perceptual_rainbow_16.hex_colors)
shuffle(colors)

repofolders = get_rep_folders(moduledir)
files = repo_file_list(moduledir, repofolders, normal=True, bones=False)

with open('astrocats/' + moduledir + '/input/non-' + modulename + '-types.json', 'r') as f:
    nontypes = json.loads(f.read(), object_pairs_hook=OrderedDict)
    nontypes = [x.upper() for x in nontypes]

for fcnt, eventfile in enumerate(tq(sorted(files, key=lambda s: s.lower()),
                                    "Collecting positions")):
    # if fcnt > 5000:
    #    break

    filetext = get_event_text(eventfile)

    thisevent = json.loads(filetext, object_pairs_hook=OrderedDict)
    thisevent = thisevent[list(thisevent.keys())[0]]

    if 'ra' in thisevent and 'dec' in thisevent:
        if 'claimedtype' in thisevent and thisevent['claimedtype']:
            for ct in [x['value'] for x in thisevent['claimedtype']]:
                thistype = ct.replace('?', '').replace('*', '')
                if thistype.upper() in nontypes:
                    continue
                elif thistype in ('Other', 'not Ia', 'SN', 'unconf', 'Radio',
                                  'CC', 'CCSN', 'Candidate', 'nIa'):
                    evtypes.append('Unknown')
                    break
                else:
                    evtypes.append(thistype)
                    break
        else:
            evtypes.append('Unknown')

        tprint(thisevent['name'])
        try:
            c = coord(ra=thisevent['ra'][0]['value'], dec=thisevent[
                      'dec'][0]['value'], unit=(un.hourangle, un.deg))
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            warnings.warn('Mangled coordinate, skipping')
            continue
        else:
            evnames.append(thisevent['name'])
            rarad = c.ra.radian - pi
            decrad = c.dec.radian
            snhx = 2.0**1.5 * cos(decrad) * sin(rarad / 2.0) / \
                sqrt(1.0 + cos(decrad) * cos(rarad / 2.0))
            snhy = sqrt(2.0) * sin(decrad) / \
                sqrt(1.0 + cos(decrad) * cos(rarad / 2.0))
            evras.append(c.ra.deg)
            evdecs.append(c.dec.deg)
            evhxs.append(snhx)
            evhys.append(snhy)

rangepts = 100
raseps = 24
decseps = 18
rarange = [-pi + i * 2.0 * pi / rangepts for i in range(0, rangepts + 1)]
decrange = [-pi / 2.0 + i * pi / rangepts for i in range(0, rangepts + 1)]
ragrid = [-pi + i * 2.0 * pi / raseps for i in range(0, raseps + 1)]
decgrid = [-pi / 2.0 + i * pi / decseps for i in range(0, decseps + 1)]

tt = [
    ("Event", "@event"),
    ("R.A. (deg)", "@ra{1.111}"),
    ("Dec. (deg)", "@dec{1.111}"),
    ("Claimed type", "@claimedtype")
]
p1 = Figure(title=moduletitle + ' Positions',
            # responsive = True,
            tools=tools, plot_width=990,
            x_range=(-1.05 * (2.0**1.5), 1.05 * 2.0**1.5),
            y_range=(-1.4 * sqrt(2.0), 1.0 * sqrt(2.0)))
p1.axis.visible = False
p1.outline_line_color = None
p1.xgrid.grid_line_color = None
p1.ygrid.grid_line_color = None
p1.title.text_font_size = '20pt'
p1.title.align = 'center'

raxs = []
rays = []
for rg in ragrid:
    raxs.append([2.0**1.5 * cos(x) * sin(rg / 2.0) /
                 sqrt(1.0 + cos(x) * cos(rg / 2.0)) for x in decrange])
    rays.append([sqrt(2.0) * sin(x) / sqrt(1.0 + cos(x) * cos(rg / 2.0))
                 for x in decrange])

decxs = []
decys = []
for dg in decgrid:
    decxs.append([2.0**1.5 * cos(dg) * sin(x / 2.0) /
                  sqrt(1.0 + cos(dg) * cos(x / 2.0)) for x in rarange])
    decys.append([sqrt(2.0) * sin(dg) / sqrt(1.0 + cos(dg) * cos(x / 2.0))
                  for x in rarange])

p1.multi_line(raxs, rays, color='#bbbbbb')
p1.multi_line(decxs, decys, color='#bbbbbb')

claimedtypes = sorted(list(set(evtypes)))

glyphs = []
glsize = max(2.5, 7.0 - np.log10(len(evtypes)))
for ci, ct in enumerate(claimedtypes):
    ind = [i for i, t in enumerate(evtypes) if t == ct]

    source = ColumnDataSource(
        data=dict(
            x=[evhxs[i] for i in ind],
            y=[evhys[i] for i in ind],
            ra=[evras[i] for i in ind],
            dec=[evdecs[i] for i in ind],
            event=[evnames[i] for i in ind],
            claimedtype=[evtypes[i] for i in ind]
        )
    )
    if ct == 'Unknown':
        tcolor = 'black'
        falpha = 0.0
    else:
        tcolor = colors[ci]
        falpha = 1.0
    glyphs.append(p1.circle('x', 'y', source=source, color=tcolor,
              fill_alpha=falpha, legend=ct, size=glsize))

hover = HoverTool(tooltips=tt, renderers=glyphs)
p1.add_tools(hover)

p1.legend.location = "bottom_center"
p1.legend.orientation = "horizontal"
p1.legend.label_text_font_size = '7pt'
p1.legend.label_width = 20
p1.legend.label_height = 8
p1.legend.glyph_height = 8
p1.legend.spacing = 0

html = file_html(p1, CDN, 'Supernova locations').replace('width: 90%;', 'width: inherit;')

with open(outdir + modulename + "-locations.html", "w") as f:
    f.write(html)
