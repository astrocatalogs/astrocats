"""Class for representing photometric data.
"""
from collections import OrderedDict
from random import seed, shuffle

from palettable import colorbrewer, cubehelix, wesanderson


# from .entry import KEYS
from .key import Key, KEY_TYPES

REQUIRE_KEY_IN_PHOTOMETRY = True


class PHOTOMETRY:
    TIME = Key('time', KEY_TYPES.NUMERIC)
    U_TIME = Key('u_time', KEY_TYPES.ANY)
    E_TIME = Key('e_time', KEY_TYPES.NUMERIC)
    TELESCOPE = Key('telescope', KEY_TYPES.ANY)
    INSTRUMENT = Key('instrument', KEY_TYPES.ANY)
    BAND = Key('band', KEY_TYPES.ANY)
    MAGNITUDE = Key('magnitude', KEY_TYPES.NUMERIC)
    E_MAGNITUDE = Key('e_magnitude', KEY_TYPES.NUMERIC)
    SOURCE = Key('source', KEY_TYPES.ANY)
    SYSTEM = Key('system', KEY_TYPES.ANY)
    SCORRECTED = Key('scorrected', KEY_TYPES.ANY)
    OBSERVATORY = Key('observatory', KEY_TYPES.ANY)
    OBSERVER = Key('observer', KEY_TYPES.ANY)
    SURVEY = Key('survey', KEY_TYPES.ANY)
    KCORRECTED = Key('kcorrected', KEY_TYPES.ANY)
    FLUX = Key('flux', KEY_TYPES.NUMERIC)
    FLUX_DENSITY = Key('fluxdensity', KEY_TYPES.NUMERIC)
    E_FLUX = Key('e_flux', KEY_TYPES.NUMERIC)
    E_FLUX_DENSITY = Key('e_fluxdensity', KEY_TYPES.NUMERIC)
    COUNTS = Key('counts', KEY_TYPES.NUMERIC)
    E_COUNTS = Key('e_counts', KEY_TYPES.NUMERIC)
    # U_FLUX = Key('u_flux', _type_)
    # U_FLUX_DENSITY = Key('u_fluxdensity', _type_)
    # FREQUENCY = Key('frequency', _type_)
    # U_FREQUENCY = Key('u_frequency', _type_)
    # NHMW = Key('nhmw', _type_)
    # PHOTON_INDEX = Key('photonindex', _type_)
    # UNABSORBED_FLUX = Key('unabsorbedflux', _type_)
    # E_UNABSORBED_FLUX = Key('e_unabsorbedflux', _type_)
    # ENERGY = Key('energy', _type_)
    # U_ENERGY = Key('u_energy', _type_)
    # E_LOWER_MAGNITUDE = Key('e_lower_magnitude', _type_)
    # E_UPPER_MAGNITUDE = Key('e_upper_magnitude', _type_)
    # E_LOWER_TIME = Key('e_lower_time', _type_)
    # E_UPPER_TIME = Key('e_upper_time', _type_)
    # MCORRECTED = Key('mcorrected', _type_)
    UPPERLIMIT = Key('upperlimit', KEY_TYPES.BOOL)
    HOST = Key('host', KEY_TYPES.BOOL)
    INCLUDESHOST = Key('includeshost', KEY_TYPES.BOOL)

    _keys = sorted([kk for kk in dir() if not kk.startswith('_')])


class Photometry(OrderedDict):
    """
    """

    def __init__(self, **kwargs):
        self._check_kwargs(kwargs)

        photo_keys = PHOTOMETRY._keys




        # Iterate over all passed keys, sanitize and store them
        for key, val in kwargs.items():
            # Make sure the key is one of the allowed keys from `PHOTOMETRY`
            if REQUIRE_KEY_IN_PHOTOMETRY and key not in photo_keys:
                raise KeyError("'{}' is not a valid `PHOTOMETRY` key".format(
                    key))

            # Make sure given value is compatible with the 'Key' specification
            # if

    def _check_kwargs(self, kwargs):
        REQ_KEY_TYPES.ANY = [
            [PHOTOMETRY.TIME, PHOTOMETRY.HOST],
            [PHOTOMETRY.MAGNITUDE, PHOTOMETRY.FLUX_DENSITY, PHOTOMETRY.FLUX,
             PHOTOMETRY.COUNTS, PHOTOMETRY.UNABSORBED_FLUX]]

        REQ_NUMBER = [PHOTOMETRY.TIME, PHOTOMETRY.MAGNITUDE, PHOTOMETRY.FLUX, PHOTOMETRY.FLUX_DENSITY, PHOTOMETRY.COUNTS,
                      PHOTOMETRY.E_MAGNITUDE, PHOTOMETRY.E_FLUX]

        for req_any in REQ_KEY_TYPES.ANY:
            if not any([req_key in kwargs for req_key in req_any]):
                err_str = "Require one of: " + ",".join(
                    "'{}'".format(rk) for rk in req_any)
                raise ValueError(err_str)



        if ((not host and not is_number(time)) or
            (not is_number(magnitude) and not is_number(flux) and not
             is_number(fluxdensity) and not is_number(counts))):
            warnings.warn('Time or brightness not numerical, not adding.')
            tprint('Name : "' + name + '", Time: "' + time + '", Band: "' +
                   band + '", AB magnitude: "' + magnitude + '"')
            return

        if (((e_magnitude and not is_number(e_magnitude)) or
             (e_flux and not is_number(e_flux)) or
             (e_fluxdensity and not is_number(e_fluxdensity)) or
             (e_counts and not is_number(e_counts)))):
            warnings.warn('Brightness error not numerical, not adding.')
            tprint('Name : "' + name + '", Time: "' + time +
                   '", Band: "' + band + '", AB error: "' + e_magnitude + '"')
            return

        if e_time and not is_number(e_time):
            warnings.warn('Time error not numerical, not adding.')
            tprint('Name : "' + name + '", Time: "' +
                   time + '", Time error: "' + e_time + '"')
            return

        if ((flux or fluxdensity) and ((not u_flux and not u_fluxdensity) or
                                       (not frequency and not band and not
                                        energy))):
            warnings.warn(
                "Unit and band/frequency must be set when adding photometry "
                "by flux or flux density, not adding.")
            tprint('Name : "' + name + '", Time: "' + time)
            return

        if not source:
            ValueError('Photometry must have source before being added!')

        if self.is_erroneous('photometry', source):
            return

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
