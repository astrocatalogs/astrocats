from palettable import cubehelix, colorbrewer, wesanderson
from collections import OrderedDict
from random import shuffle, seed

__all__ = ['bandrepf', 'bandcolorf', 'radiocolorf', 'xraycolorf', 'bandaliasf', 'bandshortaliasf',
           'bandwavef', 'bandmetaf']

BAND_REPS = {
    'Ks': ['K_s'],
    'M2': ['uvm2', 'UVM2', 'UVm2', 'Um2', 'm2', 'um2'],
    'W1': ['uvw1', 'UVW1', 'UVw1', 'Uw1', 'w1', 'uw1'],
    'W2': ['uvw2', 'UVW2', 'UVw2', 'Uw2', 'w2', 'uw2']
}

# Some bands are uniquely tied to an instrument/telescope/system, add this info here.
BAND_META = {
    'M2':     {'telescope': 'Swift', 'instrument': 'UVOT'},
    'W1':     {'telescope': 'Swift', 'instrument': 'UVOT'},
    'W2':     {'telescope': 'Swift', 'instrument': 'UVOT'},
    'F110W':  {'telescope': 'Hubble', 'instrument': 'WFC3'},
    'F775W':  {'telescope': 'Hubble', 'instrument': 'WFC3'},
    'F850LP': {'telescope': 'Hubble', 'instrument': 'WFC3'}
}

BAND_CODES = [
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
    "C",
    "CR",
    "CV",
    "M2",
    "W1",
    "W2",
    "pg",
    "Mp",
    "w",
    "y",
    "Z",
    "F110W",
    "F775W",
    "F850LP",
    "VM",
    "RM",
    "Ks"
]

BAND_ALIASES = OrderedDict([
    ("u_SDSS", "u (SDSS)"),
    ("g_SDSS", "g (SDSS)"),
    ("r_SDSS", "r (SDSS)"),
    ("i_SDSS", "i (SDSS)"),
    ("z_SDSS", "z (SDSS)")
])

BAND_ALIASES_SHORT = OrderedDict([
    ("u_SDSS", "u"),
    ("g_SDSS", "g"),
    ("r_SDSS", "r"),
    ("i_SDSS", "i"),
    ("z_SDSS", "z"),
    ("G"     , "" )
])

BAND_WAVELENGTHS = {
    "u"      : 354.,
    "g"      : 475.,
    "r"      : 622.,
    "i"      : 763.,
    "z"      : 905.,
    "u'"     : 354.,
    "g'"     : 475.,
    "r'"     : 622.,
    "i'"     : 763.,
    "z'"     : 905.,
    "u_SDSS" : 354.3,
    "g_SDSS" : 477.0,
    "r_SDSS" : 623.1,
    "i_SDSS" : 762.5,
    "z_SDSS" : 913.4,
    "U"      : 365.,
    "B"      : 445.,
    "V"      : 551.,
    "R"      : 658.,
    "I"      : 806.,
    "Y"      : 1020.,
    "J"      : 1220.,
    "H"      : 1630.,
    "K"      : 2190.,
    "M2"     : 260.,
    "W1"     : 224.6,
    "W2"     : 192.8
}

RADIO_CODES = [
    "5.9"
]
XRAY_CODES = [
    "0.3 - 10",
    "0.5 - 8"
]

seed(101)
#bandcolors = ["#%06x" % round(float(x)/float(len(BAND_CODES))*0xFFFEFF) for x in range(len(BAND_CODES))]
bandcolors = cubehelix.cubehelix1_16.hex_colors[2:13] + cubehelix.cubehelix2_16.hex_colors[2:13] + cubehelix.cubehelix3_16.hex_colors[2:13]
shuffle(bandcolors)
bandcolors2 = cubehelix.perceptual_rainbow_16.hex_colors
shuffle(bandcolors2)
bandcolors = bandcolors + bandcolors2
bandcolordict = dict(list(zip(BAND_CODES,bandcolors)))

radiocolors = wesanderson.Zissou_5.hex_colors
shuffle(radiocolors)
radiocolordict = dict(list(zip(RADIO_CODES,radiocolors)))

xraycolors = colorbrewer.sequential.Oranges_9.hex_colors[2:]
shuffle(xraycolors)
xraycolordict = dict(list(zip(XRAY_CODES,xraycolors)))

def bandrepf(code):
    for rep in BAND_REPS:
        if code in BAND_REPS[rep]:
            return rep
    return code


def bandcolorf(code):
    newcode = bandrepf(code)
    if newcode in bandcolordict:
        return bandcolordict[newcode]
    return 'black'


def radiocolorf(code):
    if code in radiocolordict:
        return radiocolordict[code]
    return 'black'


def xraycolorf(code):
    if code in xraycolordict:
        return xraycolordict[code]
    return 'black'


def bandaliasf(code):
    newcode = bandrepf(code)
    if newcode in BAND_ALIASES:
        return BAND_ALIASES[newcode]
    return newcode


def bandshortaliasf(code):
    newcode = bandrepf(code)
    if newcode in BAND_ALIASES_SHORT:
        return BAND_ALIASES_SHORT[newcode]
    return newcode


def bandwavef(code):
    newcode = bandrepf(code)
    if newcode in BAND_WAVELENGTHS:
        return BAND_WAVELENGTHS[newcode]
    return 0.


def bandmetaf(band, field):
    if band in BAND_META:
        if field in BAND_META[band]:
            return BAND_META[band][field]
    return ''
