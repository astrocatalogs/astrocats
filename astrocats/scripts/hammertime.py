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
from bokeh.models import ColumnDataSource, HoverTool, Label
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

moduletype = 'claimedtype'
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
elif args.catalog == 'kne':
    moduledir = 'kilonovae'
    modulename = 'kne'
    moduleurl = 'kilonova.space'
    moduletitle = 'Kilonova'
elif args.catalog == 'hvs':
    moduledir = 'faststars'
    modulename = 'hvs'
    moduleurl = 'faststars.space'
    moduletitle = 'Fast Stars'
    moduletype = 'spectraltype'
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
evbps = []
evpmrs = []
evpmds = []
evrvs = []
evgvs = []

seed(12483)
if moduletype == 'spectraltype':
    colors = list(reversed([
        'lightslategrey', 'firebrick', 'darkorange', 'orange',
        'gold', 'paleturquoise', 'deepskyblue', 'mediumpurple']))
    spectral_ordering = list('OBAFGKM')
    untype = 'Unclassified'
else:
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
    untype = 'Unknown'

repofolders = get_rep_folders(moduledir)
files = repo_file_list(moduledir, repofolders, normal=True, bones=True)

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

    # Code for Boubert 2018.
    # if 'discoverdate' in thisevent:
    #     if int(thisevent['discoverdate'][0]['value']) >= 2018:
    #         continue

    # if 'boundprobability' not in thisevent or float(thisevent['boundprobability'][0]['value']) > 0.5:
    #     continue

    if 'ra' in thisevent and 'dec' in thisevent:
        if moduletype in thisevent and thisevent[moduletype]:
            levtypes = []
            for ct in [x['value'] for x in thisevent[moduletype]]:
                thistype = ct.replace('?', '').replace('*', '')
                if thistype.upper() in nontypes:
                    continue

                if moduletype == 'claimedtype':
                    if thistype in ('Other', 'not Ia', 'SN', 'unconf', 'Radio',
                                      'CC', 'CCSN', 'Candidate', 'nIa'):
                        continue
                elif moduletype == 'spectraltype':
                    thistype = thistype[0]
                    if thistype not in spectral_ordering:
                        continue
                levtypes.append(thistype)

            if not len(levtypes):
                evtype = untype
            elif moduletype == 'spectraltype':
                evtype = [x for x in spectral_ordering if x in levtypes][0]
            else:
                evtype = levtypes[0]
        else:
            evtype = untype

        try:
            c = coord(ra=thisevent['ra'][0]['value'], dec=thisevent[
                      'dec'][0]['value'], unit=(un.hourangle, un.deg))
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            warnings.warn('Mangled coordinate, skipping {}.'.format(thisevent['name']))
            continue
        else:
            evtypes.append(evtype)
            evnames.append(thisevent['name'])
            rarad = c.ra.radian - pi
            decrad = c.dec.radian
            snhx = 2.0**1.5 * cos(decrad) * sin(rarad / 2.0) / \
                sqrt(1.0 + cos(decrad) * cos(rarad / 2.0))
            snhy = sqrt(2.0) * sin(decrad) / \
                sqrt(1.0 + cos(decrad) * cos(rarad / 2.0))
            rastr = str(c.ra.deg)
            decstr = str(c.dec.deg)
            if 'e_value' in thisevent['ra'][0]:
                rastr += ' ± ' + thisevent['ra'][0]['e_value'] + ' ' + thisevent['ra'][0].get('u_e_value', '')
            if 'e_value' in thisevent['dec'][0]:
                decstr += ' ± ' + thisevent['dec'][0]['e_value'] + ' ' + thisevent['dec'][0].get('u_e_value', '')
            evras.append(rastr)
            evdecs.append(decstr)
            evhxs.append(snhx)
            evhys.append(snhy)
        if args.catalog == 'hvs':
            if thisevent.get('boundprobability'):
                bpstr = str(np.round(100.0 * float(thisevent['boundprobability'][0]['value']), 3)) + '%'
                if 'upperlimit' in thisevent['boundprobability'][0]:
                    bpstr = '<' + bpstr
                evbps.append(bpstr)
            else:
                evbps.append('?')
            if thisevent.get('propermotionra'):
                bpstr = str(np.round(float(thisevent['propermotionra'][0]['value']), 4))
                if 'upperlimit' in thisevent['propermotionra'][0]:
                    bpstr = '<' + bpstr
                if 'e_value' in thisevent['propermotiondec'][0]:
                    bpstr += ' ± ' + thisevent['propermotiondec'][0]['e_value']
                evpmrs.append(bpstr)
            else:
                evpmrs.append('?')
            if thisevent.get('propermotiondec'):
                bpstr = str(np.round(float(thisevent['propermotiondec'][0]['value']), 4))
                if 'upperlimit' in thisevent['propermotiondec'][0]:
                    bpstr = '<' + bpstr
                if 'e_value' in thisevent['propermotionra'][0]:
                    bpstr += ' ± ' + thisevent['propermotionra'][0]['e_value']
                evpmds.append(bpstr)
            else:
                evpmds.append('?')
            rvs = [x for x in thisevent.get('velocity', []) if 'kind' not in x]
            if rvs:
                bpstr = str(np.round(float(rvs[0]['value']), 4))
                if 'upperlimit' in rvs[0]:
                    bpstr = '<' + bpstr
                if 'e_value' in rvs[0]:
                    bpstr += ' ± ' + rvs[0]['e_value']
                evrvs.append(bpstr)
            else:
                evrvs.append('?')
            gvs = [x for x in thisevent.get('velocity', []) if 'galactocentric' in x.get('kind', []) and 'total' in x.get('kind', [])]
            if gvs:
                bpstr = str(np.round(float(gvs[0]['value']), 4))
                if 'upperlimit' in gvs[0]:
                    bpstr = '<' + bpstr
                if 'e_value' in gvs[0]:
                    bpstr += ' ± ' + gvs[0]['e_value']
                evgvs.append(bpstr)
            else:
                evgvs.append('?')

rangepts = 100
raseps = 24
decseps = 18
rarange = [-pi + i * 2.0 * pi / rangepts for i in range(0, rangepts + 1)]
decrange = [-pi / 2.0 + i * pi / rangepts for i in range(0, rangepts + 1)]
ragrid = [-pi + i * 2.0 * pi / raseps for i in range(0, raseps + 1)]
decgrid = [-pi / 2.0 + i * pi / decseps for i in range(0, decseps + 1)]
mwgrid = [[[y.ra.radian - pi, y.dec.radian] for y in [coord(l=x*un.radian, b=0.0*un.radian, frame='galactic').icrs]][0] for x in np.linspace(-pi, pi, 200)]
mwcenter = [[y.ra.radian - pi, y.dec.radian] for y in [coord(l=0.0*un.radian, b=0.0*un.radian, frame='galactic').icrs]][0]
anticenter = [[y.ra.radian - pi, y.dec.radian] for y in [coord(l=pi*un.radian, b=0.0*un.radian, frame='galactic').icrs]][0]
m31center = [[y.ra.radian - pi, y.dec.radian] for y in [coord('121.1743 -21.5733', unit=un.deg, frame='galactic').icrs]][0]
lmccenter = [[y.ra.radian - pi, y.dec.radian] for y in [coord('280.4652 -32.8884', unit=un.deg, frame='galactic').icrs]][0]

for mi, mg in enumerate(mwgrid):
    if mi == 0:
        continue
    if abs(mg[0] - mwgrid[mi - 1][0]) > 1.0:
        mgi = mi
        break
mwgrid = mwgrid[mgi:] + mwgrid[:mgi]

tt = [
    ("Star" if args.catalog == 'hvs' else "Event", "@event"),
    ("Type", "@moduletype"),
    ("R.A. (deg)", "@ra{1.111}"),
    ("Dec. (deg)", "@dec{1.111}")
]
if args.catalog == 'hvs':
    tt.append(("R.A. proper motion (mas/yr)", "@propermotionra"))
    tt.append(("Dec. proper motion (mas/yr)", "@propermotiondec"))
    tt.append(("Radial velocity (km/s)", "@velocity"))
    tt.append(("Galactocentric velocity (km/s)", "@galactocentricvelocity"))
    tt.append(("Bound probability", "@boundprobability"))
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

p1.multi_line(raxs, rays, color='#bbbbbb')

decxs = []
decys = []
for dg in decgrid:
    decxs.append([2.0**1.5 * cos(dg) * sin(x / 2.0) /
                  sqrt(1.0 + cos(dg) * cos(x / 2.0)) for x in rarange])
    decys.append([sqrt(2.0) * sin(dg) / sqrt(1.0 + cos(dg) * cos(x / 2.0))
                  for x in rarange])

p1.multi_line(decxs, decys, color='#bbbbbb')

mwxs = []
mwys = []
for mg in mwgrid:
    mwxs.append(2.0**1.5 * cos(mg[1]) * sin(mg[0] / 2.0) / sqrt(1.0 + cos(mg[1]) * cos(mg[0] / 2.0)))
    mwys.append(sqrt(2.0) * sin(mg[1]) / sqrt(1.0 + cos(mg[1]) * cos(mg[0] / 2.0)))

mcx = 2.0**1.5 * cos(mwcenter[1]) * sin(mwcenter[0] / 2.0) / sqrt(1.0 + cos(mwcenter[1]) * cos(mwcenter[0] / 2.0))
mcy = sqrt(2.0) * sin(mwcenter[1]) / sqrt(1.0 + cos(mwcenter[1]) * cos(mwcenter[0] / 2.0))

acx = 2.0**1.5 * cos(anticenter[1]) * sin(anticenter[0] / 2.0) / sqrt(1.0 + cos(anticenter[1]) * cos(anticenter[0] / 2.0))
acy = sqrt(2.0) * sin(anticenter[1]) / sqrt(1.0 + cos(anticenter[1]) * cos(anticenter[0] / 2.0))

m31x = 2.0**1.5 * cos(m31center[1]) * sin(m31center[0] / 2.0) / sqrt(1.0 + cos(m31center[1]) * cos(m31center[0] / 2.0))
m31y = sqrt(2.0) * sin(m31center[1]) / sqrt(1.0 + cos(m31center[1]) * cos(m31center[0] / 2.0))

lmcx = 2.0**1.5 * cos(lmccenter[1]) * sin(lmccenter[0] / 2.0) / sqrt(1.0 + cos(lmccenter[1]) * cos(lmccenter[0] / 2.0))
lmcy = sqrt(2.0) * sin(lmccenter[1]) / sqrt(1.0 + cos(lmccenter[1]) * cos(lmccenter[0] / 2.0))

p1.line(mwxs, mwys, color='#555555', line_width=2)
p1.circle(mcx, mcy, color='#555555', size=10)
p1.x(acx, acy, color='#555555', size=10, line_width=4)
p1.circle(m31x, m31y, color='#999999', size=20, alpha=0.5)
p1.circle(lmcx, lmcy, color='#999999', size=20, alpha=0.5)
label_size = '20pt'
label = Label(x=m31x, y=m31y + 0.06, text='M31', text_color='#999999', text_align='center', text_font_size=label_size)
p1.add_layout(label)
label = Label(x=mcx + 0.06, y=mcy - 0.20, text='MW center', text_color='#555555', text_font_size=label_size)
p1.add_layout(label)
label = Label(x=acx - 0.06, y=acy, text='MW anticenter', text_color='#555555', text_font_size=label_size, text_align='right')
p1.add_layout(label)
label = Label(x=lmcx, y=lmcy + 0.06, text='LMC', text_color='#999999', text_align='center', text_font_size=label_size)
p1.add_layout(label)

if moduletype == 'spectraltype':
    claimedtypes = [[x, i] for i, x in enumerate(spectral_ordering) if x in evtypes] + [[untype, len(colors) - 1]]
    colors = [colors[ct[1]] for ct in claimedtypes]
    claimedtypes = [x[0] for x in claimedtypes]
else:
    claimedtypes = sorted(list(set(evtypes)))

glyphs = []
glsize = max(2.5, 7.0 - np.log10(len(evtypes)))
for ci, ct in enumerate(claimedtypes):
    ind = [i for i, t in enumerate(evtypes) if t == ct]

    cdsdict = dict(
        x=[evhxs[i] for i in ind],
        y=[evhys[i] for i in ind],
        ra=[evras[i] for i in ind],
        dec=[evdecs[i] for i in ind],
        event=[evnames[i] for i in ind],
        moduletype=[evtypes[i] for i in ind]
    )
    if args.catalog == 'hvs':
        cdsdict['boundprobability'] = [evbps[i] for i in ind]
        cdsdict['propermotionra'] = [evpmrs[i] for i in ind]
        cdsdict['propermotiondec'] = [evpmds[i] for i in ind]
        cdsdict['velocity'] = [evrvs[i] for i in ind]
        cdsdict['galactocentricvelocity'] = [evgvs[i] for i in ind]
    source = ColumnDataSource(data=cdsdict)
    if ct == 'Unknown':
        tcolor = 'black'
        falpha = 0.0
    elif ct == 'Unclassified':
        tcolor = 'lightslategrey'
        falpha = 1.0
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
