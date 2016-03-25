from palettable import cubehelix
from collections import OrderedDict
from random import shuffle, seed, randint, uniform

bandreps = {
    'M2':     ['uvm2', 'UVM2'],
    'W1':     ['uvw1', 'UVW1'],
    'W2':     ['uvw2', 'UVW2']
}

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
    "C",
    "CR",
    "CV",
    "M2",
    "W1",
    "W2",
    "pg",
    "Mp"
]

bandaliases = OrderedDict([
    ("u_SDSS", "u (SDSS)"),
    ("g_SDSS", "g (SDSS)"),
    ("r_SDSS", "r (SDSS)"),
    ("i_SDSS", "i (SDSS)"),
    ("z_SDSS", "z (SDSS)"),
    ("uvm2"  , "M2 (UVOT)"),
    ("uvw1"  , "W1 (UVOT)"),
    ("uvw2"  , "W2 (UVOT)"),
    ("UVM2"  , "M2 (UVOT)"),
    ("UVW1"  , "W1 (UVOT)"),
    ("UVW2"  , "W2 (UVOT)")
])

bandshortaliases = OrderedDict([
    ("u_SDSS", "u"),
    ("g_SDSS", "g"),
    ("r_SDSS", "r"),
    ("i_SDSS", "i"),
    ("z_SDSS", "z"),
    ("G"     , "" )
])

bandwavelengths = {
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
    "uvm2"   : 260.,
    "uvw1"   : 224.6,
    "uvw2"   : 192.8,
    "UVM2"   : 260.,
    "UVW1"   : 224.6,
    "UVW2"   : 192.8
}

seed(101)
#bandcolors = ["#%06x" % round(float(x)/float(len(bandcodes))*0xFFFEFF) for x in range(len(bandcodes))]
bandcolors = cubehelix.cubehelix1_16.hex_colors[2:13] + cubehelix.cubehelix2_16.hex_colors[2:13] + cubehelix.cubehelix3_16.hex_colors[2:13]
shuffle(bandcolors)

bandcolordict = dict(list(zip(bandcodes,bandcolors)))

def bandcolorf(code):
    if code in bandcolordict:
        return bandcolordict[code]
    return 'black'

def bandaliasf(code):
    if code in bandaliases:
        return bandaliases[code]
    return code

def bandshortaliasf(code):
    if code in bandshortaliases:
        return bandshortaliases[code]
    return code

def bandwavef(code):
    if code in bandwavelengths:
        return bandwavelengths[code]
    return 0.
