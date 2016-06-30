"""Class for representing photometric data.
"""
from collections import OrderedDict
from random import seed, shuffle

from palettable import colorbrewer, cubehelix, wesanderson

# from .entry import KEYS
from .key import Key, KEY_TYPES, KeyCollection

# If `REQUIRE_KEY_IN_PHOTOMETRY` is 'True', then only parameters with names
#    included in `PHOTOMETRY` are allowed.  Others will raise an error.
#    If this parameter is 'False', then parameters corresponding to those in
#    `PHOTOMETRY` are still checked (for type etc), but additional parameters
#    are just tacked onto the `Photometry` object without any checks or errors.
REQUIRE_KEY_IN_PHOTOMETRY = True


class PHOTOMETRY(KeyCollection):
    TIME = Key('time', KEY_TYPES.NUMERIC)
    E_TIME = Key('e_time', KEY_TYPES.NUMERIC)
    MAGNITUDE = Key('magnitude', KEY_TYPES.NUMERIC)
    FLUX = Key('flux', KEY_TYPES.NUMERIC)
    FLUX_DENSITY = Key('fluxdensity', KEY_TYPES.NUMERIC)
    COUNTS = Key('counts', KEY_TYPES.NUMERIC)

    FREQUENCY = Key('frequency', KEY_TYPES.NUMERIC)
    NHMW = Key('nhmw', KEY_TYPES.NUMERIC)
    PHOTON_INDEX = Key('photonindex', KEY_TYPES.NUMERIC)
    UNABSORBED_FLUX = Key('unabsorbedflux', KEY_TYPES.NUMERIC)
    ENERGY = Key('energy', KEY_TYPES.NUMERIC)

    E_MAGNITUDE = Key('e_magnitude', KEY_TYPES.NUMERIC)
    E_FLUX = Key('e_flux', KEY_TYPES.NUMERIC)
    E_FLUX_DENSITY = Key('e_fluxdensity', KEY_TYPES.NUMERIC)
    E_COUNTS = Key('e_counts', KEY_TYPES.NUMERIC)
    E_UNABSORBED_FLUX = Key('e_unabsorbedflux', KEY_TYPES.NUMERIC)
    E_LOWER_MAGNITUDE = Key('e_lower_magnitude', KEY_TYPES.NUMERIC)
    E_UPPER_MAGNITUDE = Key('e_upper_magnitude', KEY_TYPES.NUMERIC)
    E_LOWER_TIME = Key('e_lower_time', KEY_TYPES.NUMERIC)
    E_UPPER_TIME = Key('e_upper_time', KEY_TYPES.NUMERIC)

    TELESCOPE = Key('telescope', KEY_TYPES.STRING)
    INSTRUMENT = Key('instrument', KEY_TYPES.STRING)
    BAND = Key('band', KEY_TYPES.STRING)
    SOURCE = Key('source', KEY_TYPES.STRING)
    SYSTEM = Key('system', KEY_TYPES.STRING)
    OBSERVATORY = Key('observatory', KEY_TYPES.STRING)
    OBSERVER = Key('observer', KEY_TYPES.STRING)
    SURVEY = Key('survey', KEY_TYPES.STRING)

    U_TIME = Key('u_time', KEY_TYPES.STRING)
    U_FLUX = Key('u_flux', KEY_TYPES.STRING)
    U_FLUX_DENSITY = Key('u_fluxdensity', KEY_TYPES.STRING)
    U_FREQUENCY = Key('u_frequency', KEY_TYPES.STRING)
    U_ENERGY = Key('u_energy', KEY_TYPES.STRING)

    SCORRECTED = Key('scorrected', KEY_TYPES.BOOL)
    KCORRECTED = Key('kcorrected', KEY_TYPES.BOOL)
    MCORRECTED = Key('mcorrected', KEY_TYPES.BOOL)
    UPPER_LIMIT = Key('upperlimit', KEY_TYPES.BOOL)
    HOST = Key('host', KEY_TYPES.BOOL)
    INCLUDES_HOST = Key('includeshost', KEY_TYPES.BOOL)


class Photometry(OrderedDict):
    """
    """

    def __init__(self, **kwargs):
        # Iterate over all `PHOTOMETRY` parameters, load each if given
        #    note that the stored `values` are the `Key` objects, referred to
        #    here with the name 'key'
        for key in PHOTOMETRY.vals():
            # If this key is given, process and store it
            if key in kwargs:
                if not key.check(kwargs[key]):
                    raise ValueError("Value for '{}' is invalid '{}'".format(
                        repr(key), kwargs[key]))

                # Handle Special Cases
                # --------------------
                # Only keep booleans if they are true
                if key.type == KEY_TYPES.BOOL and not kwargs[key]:
                    del kwargs[key]
                    continue

                # Check and store values
                # ----------------------
                # Remove key-value pair from `kwargs` dictionary
                value = kwargs.pop(key)
                # Make sure value is compatible with the 'Key' specification
                if key.check(value):
                    self[key] = value
                else:
                    raise ValueError("Value for '{}' is invalid '{}'".format(
                        repr(key), value))

        # If we require all parameters to be a key in `PHOTOMETRY`, then all
        #    elements should have been removed from `kwargs`
        if REQUIRE_KEY_IN_PHOTOMETRY and len(kwargs):
            raise ValueError(
                "All permitted keys stored, remaining: '{}'".format(kwargs))

        # Make sure that currently stored values are valid
        self._check()

        return

    def __repr__(self):
        pass

    def _check(self):
        REQ_KEY_TYPES = [
            [PHOTOMETRY.TIME, PHOTOMETRY.HOST],
            [PHOTOMETRY.MAGNITUDE, PHOTOMETRY.FLUX, PHOTOMETRY.FLUX_DENSITY,
             PHOTOMETRY.COUNTS]]

        for req_any in REQ_KEY_TYPES:
            if not any([req_key in self for req_key in req_any]):
                err_str = "Require one of: " + ",".join(
                    "'{}'".format(rk) for rk in req_any)
                raise ValueError(err_str)

        if PHOTOMETRY.SOURCE not in self:
            raise ValueError('Photometry must have source before being added!')

        return


BAND_REPS = {
    'Ks': ['K_s'],
    'M2': ['uvm2', 'UVM2', 'UVm2', 'Um2', 'm2', 'um2'],
    'W1': ['uvw1', 'UVW1', 'UVw1', 'Uw1', 'w1', 'uw1'],
    'W2': ['uvw2', 'UVW2', 'UVw2', 'Uw2', 'w2', 'uw2']
}

# Some bands are uniquely tied to an instrument/telescope/system, add this
# info here.
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
    ("G", "")
])

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
    "M2": 260.,
    "W1": 224.6,
    "W2": 192.8
}

RADIO_CODES = [
    "5.9"
]
XRAY_CODES = [
    "0.3 - 10",
    "0.5 - 8"
]

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
