"""
"""

from collections import OrderedDict

from decimal import Decimal, localcontext
from random import seed, shuffle
from palettable import colorbrewer, cubehelix, wesanderson

from astrocats.catalog.utils import get_sig_digits
from astrocats.catalog.photometry import PHOTOMETRY


DEFAULT_UL_SIGMA = 5.0
DEFAULT_ZP = 30.0
D25 = Decimal('2.5')

BAND_REPS = {
    'Ks': ['K_s'],
    'UVM2': ['uvm2', 'UVM2', 'UVm2', 'Um2', 'um2'],
    'UVW1': ['uvw1', 'UVW1', 'UVw1', 'Uw1', 'uw1'],
    'UVW2': ['uvw2', 'UVW2', 'UVw2', 'Uw2', 'uw2']
}

# Some bands are uniquely tied to an instrument/telescope/system, add this
# info here.
BAND_META = {
    'UVM2': {
        PHOTOMETRY.TELESCOPE: 'Swift',
        PHOTOMETRY.INSTRUMENT: 'UVOT'
    },
    'UVW1': {
        PHOTOMETRY.TELESCOPE: 'Swift',
        PHOTOMETRY.INSTRUMENT: 'UVOT'
    },
    'UVW2': {
        PHOTOMETRY.TELESCOPE: 'Swift',
        PHOTOMETRY.INSTRUMENT: 'UVOT'
    },
    'F110W': {
        PHOTOMETRY.TELESCOPE: 'Hubble',
        PHOTOMETRY.INSTRUMENT: 'WFC3'
    },
    'F775W': {
        PHOTOMETRY.TELESCOPE: 'Hubble',
        PHOTOMETRY.INSTRUMENT: 'WFC3'
    },
    'F850LP': {
        PHOTOMETRY.TELESCOPE: 'Hubble',
        PHOTOMETRY.INSTRUMENT: 'WFC3'
    },
    'Kepler': {
        PHOTOMETRY.TELESCOPE: 'Kepler',
        PHOTOMETRY.INSTRUMENT: 'Kepler'
    },
    'G': {
        PHOTOMETRY.TELESCOPE: 'Gaia',
        PHOTOMETRY.INSTRUMENT: 'Astrometric'
    }
}

BAND_CODES = [
    "u", "g", "r", "i", "z", "u'", "g'", "r'", "i'", "z'", "u_SDSS", "g_SDSS",
    "r_SDSS", "i_SDSS", "z_SDSS", "U", "B", "V", "R", "I", "G", "Y", "J", "H",
    "K", "C", "CR", "CV", "M2", "UVW1", "UVW2", "pg", "Mp", "w", "y", "Z",
    "F110W", "F775W", "F850LP", "VM", "RM", "Ks", "Ic", "Rc"
]

BAND_ALIASES = OrderedDict([("u_SDSS", "u (SDSS)"), ("g_SDSS", "g (SDSS)"),
                            ("r_SDSS", "r (SDSS)"), ("i_SDSS", "i (SDSS)"),
                            ("z_SDSS", "z (SDSS)")])

BAND_ALIASES_SHORT = OrderedDict([("u_SDSS", "u"), ("g_SDSS", "g"),
                                  ("r_SDSS", "r"), ("i_SDSS", "i"),
                                  ("z_SDSS", "z"), ("G", "")])

BAND_WAVELENGTHS = {
    "u": 354.,
    "g": 475.,
    "r": 622.,
    "i": 763.,
    "z": 905.,
    "u'": 354.,
    "g'": 475.,
    "r'": 622.,
    "i'": 763.,
    "z'": 905.,
    "u_SDSS": 354.3,
    "g_SDSS": 477.0,
    "r_SDSS": 623.1,
    "i_SDSS": 762.5,
    "z_SDSS": 913.4,
    "U": 365.,
    "B": 445.,
    "V": 551.,
    "R": 658.,
    "I": 806.,
    "Y": 1020.,
    "J": 1220.,
    "H": 1630.,
    "K": 2190.,
    "UVM2": 260.,
    "UVW1": 224.6,
    "UVW2": 192.8,
    "Ic": 786.5,
    "Rc": 647.
}

RADIO_CODES = ["5.9"]
XRAY_CODES = ["0.3 - 10", "0.5 - 8"]

INSTRUMENT_REPS = {
    'Astrometric': 'Gaia-photometric'
}

seed(101)
# bandcolors = ["#%06x" % round(float(x)/float(len(BAND_CODES))*0xFFFEFF)
# for x in range(len(BAND_CODES))]
bandcolors = (cubehelix.cubehelix1_16.hex_colors[2:13] +
              cubehelix.cubehelix2_16.hex_colors[2:13] +
              cubehelix.cubehelix3_16.hex_colors[2:13])
shuffle(bandcolors)
bandcolors2 = cubehelix.perceptual_rainbow_16.hex_colors
shuffle(bandcolors2)
bandcolors = bandcolors + bandcolors2
bandcolordict = dict(list(zip(BAND_CODES, bandcolors)))

radiocolors = wesanderson.Zissou_5.hex_colors
shuffle(radiocolors)
radiocolordict = dict(list(zip(RADIO_CODES, radiocolors)))

xraycolors = colorbrewer.sequential.Oranges_9.hex_colors[2:]
shuffle(xraycolors)
xraycolordict = dict(list(zip(XRAY_CODES, xraycolors)))


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
    return None


def set_pd_mag_from_counts(photodict, c='', ec='', lec='', uec='',
                           zp=DEFAULT_ZP, sig=DEFAULT_UL_SIGMA):
    """Set photometry dictionary from a counts measurement."""
    with localcontext() as ctx:
        if lec == '' or uec == '':
            lec = ec
            uec = ec
        prec = max(
            get_sig_digits(str(c), strip_zeroes=False),
            get_sig_digits(str(lec), strip_zeroes=False),
            get_sig_digits(str(uec), strip_zeroes=False)) + 1
        ctx.prec = prec
        dlec = Decimal(str(lec))
        duec = Decimal(str(uec))
        if c != '':
            dc = Decimal(str(c))
        dzp = Decimal(str(zp))
        dsig = Decimal(str(sig))
        photodict[PHOTOMETRY.ZERO_POINT] = str(zp)
        if c == '' or float(c) < float(sig) * float(uec):
            photodict[PHOTOMETRY.UPPER_LIMIT] = True
            photodict[PHOTOMETRY.UPPER_LIMIT_SIGMA] = str(sig)
            photodict[PHOTOMETRY.MAGNITUDE] = str(dzp - (D25 * (dsig * duec
                                                                ).log10()))
            dnec = Decimal('10.0') ** (
                (dzp - Decimal(photodict[PHOTOMETRY.MAGNITUDE])) / D25)
            photodict[PHOTOMETRY.E_UPPER_MAGNITUDE] = str(D25 * (
                (dnec + duec).log10() - dnec.log10()))
        else:
            photodict[PHOTOMETRY.MAGNITUDE] = str(dzp - D25 * dc.log10())
            photodict[PHOTOMETRY.E_UPPER_MAGNITUDE] = str(D25 * (
                (dc + duec).log10() - dc.log10()))
            photodict[PHOTOMETRY.E_LOWER_MAGNITUDE] = str(D25 * (
                dc.log10() - (dc - dlec).log10()))
