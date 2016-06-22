#!/usr/local/bin/python3.5

import json
import re
import sys
import os
import seaborn as sns
from tq import *
from digits import *
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.plotting import reset_output, Figure
from bokeh.resources import CDN
from bokeh.embed import file_html
from collections import OrderedDict
from palettable import colorbrewer

outdir = "../"

mincnt = 5

with open('../catalog.min.json', 'r') as f:
    filetext = f.read()
    meta = json.loads(filetext, object_pairs_hook=OrderedDict)
with open('type-synonyms.json', 'r') as f:
    typereps = json.loads(f.read(), object_pairs_hook=OrderedDict)
with open('non-sne-types.json', 'r') as f:
    nonsnetypes = json.loads(f.read(), object_pairs_hook=OrderedDict)
    nonsnetypes = [x.upper() for x in nonsnetypes]

sntypes = []

for event in meta:
    if 'claimedtype' in event and event['claimedtype']:
        for ct in event['claimedtype']:
            ctv = ct['value'].strip('?* ')
            for rep in typereps:
                if ctv in typereps[rep]:
                    ctv = rep
                    break
            if not ctv:
                continue
            if ctv not in sntypes and ctv.upper() not in nonsnetypes and ctv not in ['nIa'] \
                and not is_number(ctv) and '\\' not in ctv: #temporarily ignoring bad types from import
                sntypes.append(ctv) 

sntypes = sorted(sntypes)
snoffs = [[] for x in range(len(sntypes))]

tt = [  
        ("Type", "@ct")
     ]
hover = HoverTool(tooltips = tt, line_policy = 'interp')
p = Figure(x_range = [0., 100.], y_range = [0., 1.], title = 'Supernova Host Offsets',
    x_axis_label = 'Offset (kpc)', y_axis_label = 'CDF', plot_width = 980, plot_height = 500,
    title_text_font = 'futura', title_text_font_size = '14pt')
p.add_tools(hover)

for si, sntype in enumerate(sntypes):
    for event in meta:
        if 'hostoffkpc' in event and is_number(event['hostoffkpc']) \
            and 'claimedtype' in event and event['claimedtype'] and sntype in [x['value'] for x in event['claimedtype']]:
            snoffs[si].append(float(event['hostoffkpc']))
    snoffs[si] = sorted(snoffs[si])

colors = sns.color_palette("hls", n_colors = sum([1 if len(snoffs[i]) >= mincnt else 0 for i, x in enumerate(snoffs)])).as_hex()

cnt = 0
for si, sntype in enumerate(sntypes):
    if len(snoffs[si]) >= mincnt:
        data = dict(
            x = snoffs[si],
            y = [float(x)/float(len(snoffs[si])-1) for x in range(len(snoffs[si]))],
            ct = [sntype for x in range(len(snoffs[si]))]
        )
        p.line('x', 'y', source = ColumnDataSource(data), legend = sntype, color = colors[cnt], line_width = 1.5)
        cnt = cnt + 1

p.xaxis.axis_label_text_font = 'futura'
p.yaxis.axis_label_text_font = 'futura'
p.xaxis.major_label_text_font = 'futura'
p.yaxis.major_label_text_font = 'futura'
p.xaxis.axis_label_text_font_size = '11pt'
p.yaxis.axis_label_text_font_size = '11pt'
p.xaxis.major_label_text_font_size = '8pt'
p.yaxis.major_label_text_font_size = '8pt'
p.legend.label_text_font = 'futura'
p.legend.location = "bottom_right"
p.legend.label_text_font_size = '8pt'
p.legend.label_width = 20
p.legend.label_height = 10
p.legend.glyph_height = 10
p.legend.legend_spacing = 3

html = file_html(p, CDN, 'Offsets')

with open(outdir + "host-offsets.html", "w") as f:
    f.write(html)

# Necessary to clear Bokeh state
reset_output()
