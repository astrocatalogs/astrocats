"""Class for representing photometric data."""
from collections import OrderedDict
from random import seed, shuffle

from astropy.time import Time as astrotime
from palettable import colorbrewer, cubehelix, wesanderson

from astrocats.catalog.catdict import CatDict, CatDictError
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection
from astrocats.catalog.utils import get_sig_digits, listify
from decimal import Decimal, localcontext

DEFAULT_UL_SIGMA = 3.0
DEFAULT_ZP = 30.0
D25 = Decimal('2.5')


class PHOTOMETRY(KeyCollection):
    """Keys for the `Photometry` class."""

    TIME = Key('time', KEY_TYPES.TIME, listable=True, priority=10)
    MAGNITUDE = Key('magnitude', KEY_TYPES.NUMERIC, priority=9)
    FLUX = Key('flux', KEY_TYPES.NUMERIC)
    FLUX_DENSITY = Key('fluxdensity', KEY_TYPES.NUMERIC)
    COUNT_RATE = Key('countrate', KEY_TYPES.NUMERIC)
    LUMINOSITY = Key('luminosity', KEY_TYPES.NUMERIC)
    ZERO_POINT = Key('zeropoint', KEY_TYPES.NUMERIC)
    UPPER_LIMIT_SIGMA = Key('upperlimitsigma', KEY_TYPES.NUMERIC)

    ENERGY = Key('energy', KEY_TYPES.NUMERIC, listable=True)
    FREQUENCY = Key('frequency', KEY_TYPES.NUMERIC, listable=True)
    WAVELENGTH = Key('wavelength', KEY_TYPES.NUMERIC, listable=True)
    NHMW = Key('nhmw', KEY_TYPES.NUMERIC)
    PHOTON_INDEX = Key('photonindex', KEY_TYPES.NUMERIC)
    UNABSORBED_FLUX = Key('unabsorbedflux', KEY_TYPES.NUMERIC)
    EXPOSURE_TIME = Key('exposuretime', KEY_TYPES.NUMERIC)
    OFF_AXIS_ANGLE = Key('offaxisangle', KEY_TYPES.NUMERIC)
    EXTRACTION_RADIUS = Key('extractionradius', KEY_TYPES.NUMERIC)

    E_COUNT_RATE = Key('e_countrate', KEY_TYPES.NUMERIC)
    E_FLUX = Key('e_flux', KEY_TYPES.NUMERIC)
    E_FLUX_DENSITY = Key('e_fluxdensity', KEY_TYPES.NUMERIC)
    E_LUMINOSITY = Key('e_luminosity', KEY_TYPES.NUMERIC)
    E_MAGNITUDE = Key('e_magnitude', KEY_TYPES.NUMERIC, priority=7)
    E_TIME = Key('e_time', KEY_TYPES.NUMERIC)
    E_UNABSORBED_FLUX = Key('e_unabsorbedflux', KEY_TYPES.NUMERIC)
    E_LOWER_COUNT_RATE = Key('e_lower_countrate', KEY_TYPES.NUMERIC)
    E_UPPER_COUNT_RATE = Key('e_upper_countrate', KEY_TYPES.NUMERIC)
    E_LOWER_MAGNITUDE = Key('e_lower_magnitude', KEY_TYPES.NUMERIC)
    E_UPPER_MAGNITUDE = Key('e_upper_magnitude', KEY_TYPES.NUMERIC)
    E_LOWER_FLUX = Key('e_lower_flux', KEY_TYPES.NUMERIC)
    E_UPPER_FLUX_DENSITY = Key('e_upper_fluxdensity', KEY_TYPES.NUMERIC)
    E_LOWER_FLUX_DENSITY = Key('e_lower_fluxdensity', KEY_TYPES.NUMERIC)
    E_UPPER_FLUX = Key('e_upper_flux', KEY_TYPES.NUMERIC)
    E_LOWER_LUMINOSITY = Key('e_lower_luminosity', KEY_TYPES.NUMERIC)
    E_UPPER_LUMINOSITY = Key('e_upper_luminosity', KEY_TYPES.NUMERIC)
    E_LOWER_TIME = Key('e_lower_time', KEY_TYPES.NUMERIC)
    E_UPPER_TIME = Key('e_upper_time', KEY_TYPES.NUMERIC)

    MODEL = Key('model', KEY_TYPES.STRING, compare=False)
    REALIZATION = Key('realization', KEY_TYPES.STRING, priority=15)
    SOURCE = Key('source', KEY_TYPES.STRING, compare=False)
    TELESCOPE = Key('telescope', KEY_TYPES.STRING, compare=False)
    INSTRUMENT = Key('instrument', KEY_TYPES.STRING, compare=False)
    MODE = Key('mode', KEY_TYPES.STRING, compare=False)
    BAND = Key('band', KEY_TYPES.STRING, priority=8)
    OBSERVATORY = Key('observatory', KEY_TYPES.STRING, compare=False)
    OBSERVER = Key('observer', KEY_TYPES.STRING, compare=False)
    SURVEY = Key('survey', KEY_TYPES.STRING, compare=False)
    BAND_SET = Key('bandset', KEY_TYPES.STRING)
    SYSTEM = Key('system', KEY_TYPES.STRING)

    DESCRIPTION = Key('description', KEY_TYPES.STRING, compare=False)

    U_COUNT_RATE = Key('u_countrate', KEY_TYPES.STRING)
    U_TIME = Key('u_time', KEY_TYPES.STRING)
    U_FLUX = Key('u_flux', KEY_TYPES.STRING)
    U_FLUX_DENSITY = Key('u_fluxdensity', KEY_TYPES.STRING)
    U_FREQUENCY = Key('u_frequency', KEY_TYPES.STRING)
    U_WAVELENGTH = Key('u_wavelength', KEY_TYPES.STRING)
    U_ENERGY = Key('u_energy', KEY_TYPES.STRING)
    U_LUMINOSITY = Key('u_luminosity', KEY_TYPES.STRING)
    U_EXPOSURE_TIME = Key('u_exposuretime', KEY_TYPES.STRING)
    U_OFF_AXIS_ANGLE = Key('u_offaxisangle', KEY_TYPES.STRING)
    U_EXTRACTION_RADIUS = Key('u_extractionradius', KEY_TYPES.STRING)

    SCORRECTED = Key('scorrected', KEY_TYPES.BOOL)
    KCORRECTED = Key('kcorrected', KEY_TYPES.BOOL)
    MCORRECTED = Key('mcorrected', KEY_TYPES.BOOL)
    SYNTHETIC = Key('synthetic', KEY_TYPES.BOOL)
    SIMULATED = Key('simulated', KEY_TYPES.BOOL)
    UPPER_LIMIT = Key('upperlimit', KEY_TYPES.BOOL, priority=6)
    LOWER_LIMIT = Key('lowerlimit', KEY_TYPES.BOOL)
    HOST = Key('host', KEY_TYPES.BOOL)
    INCLUDES_HOST = Key('includeshost', KEY_TYPES.BOOL)
    REST_FRAME = Key('restframe', KEY_TYPES.BOOL)
    HOST_NH_CORR = Key('hostnhcorr', KEY_TYPES.BOOL)


class Photometry(CatDict):
    """Container for a single photometric point with associated metadata.

    `Source` citation required.
    Photometry can be given as [magnitude, flux, flux-density, counts,
    luminosity].
    """

    _ALLOW_UNKNOWN_KEYS = True
    _KEYS = PHOTOMETRY

    def __init__(self, parent, **kwargs):
        """Initialize."""
        self._REQ_KEY_SETS = [[PHOTOMETRY.SOURCE, PHOTOMETRY.MODEL],
                              [PHOTOMETRY.TIME, PHOTOMETRY.HOST], [
                                  PHOTOMETRY.MAGNITUDE, PHOTOMETRY.FLUX,
                                  PHOTOMETRY.FLUX_DENSITY, PHOTOMETRY.COUNT_RATE,
                                  PHOTOMETRY.LUMINOSITY]]
        # Note: `_check()` is called at end of `super().__init__`
        super(Photometry, self).__init__(parent, **kwargs)

        # If `BAND` is given, but any of `bandmetaf_keys` is not, try to infer
        if self._KEYS.BAND in self:
            sband = self[self._KEYS.BAND]
            bandmetaf_keys = [
                self._KEYS.INSTRUMENT, self._KEYS.TELESCOPE, self._KEYS.SYSTEM
            ]

            for bmf in bandmetaf_keys:
                if bmf not in self:
                    temp = bandmetaf(sband, bmf)
                    if temp is not None:
                        self[bmf] = temp

        # Convert dates to MJD
        timestrs = [str(x) for x in listify(self.get(self._KEYS.TIME, ''))]
        for ti, timestr in enumerate(timestrs):
            if (any(x in timestr for x in ['-', '/'])
                    and not timestr.startswith('-')):
                timestrs[ti] = timestr.replace('/', '-')
                try:
                    timestrs[ti] = str(
                        astrotime(timestrs[ti], format='isot').mjd)
                except Exception:
                    raise CatDictError('Unable to convert date to MJD.')
            elif timestr:  # Make sure time is string
                timestrs[ti] = timestr
        if len(timestrs) > 0 and timestrs[0] != '':
            self[self._KEYS.TIME] = timestrs if len(
                timestrs) > 1 else timestrs[0]

        # Time unit is necessary for maximum time determination
        if self._KEYS.U_TIME not in self and self._KEYS.TIME in self:
            self._log.info('`{}` not found in photometry, assuming '
                           ' MJD.'.format(self._KEYS.U_TIME))
            self[self._KEYS.U_TIME] = 'MJD'

        if (self._KEYS.U_COUNT_RATE not in self and
                self._KEYS.COUNT_RATE in self):
            self._log.info('`{}` not found in photometry, assuming '
                           ' s^-1.'.format(self._KEYS.U_COUNT_RATE))
            self[self._KEYS.U_COUNT_RATE] = 's^-1'

        return

    def _check(self):
        """Check that entry attributes are legal."""
        # Run the super method
        super(Photometry, self)._check()

        err_str = None
        has_flux = self._KEYS.FLUX in self
        has_flux_dens = self._KEYS.FLUX_DENSITY in self
        has_u_flux = self._KEYS.U_FLUX in self
        has_u_flux_dens = self._KEYS.U_FLUX_DENSITY in self

        has_freq = self._KEYS.FREQUENCY in self
        has_band = self._KEYS.BAND in self
        has_ener = self._KEYS.ENERGY in self
        has_u_freq = self._KEYS.U_FREQUENCY in self
        has_u_ener = self._KEYS.U_ENERGY in self

        if has_flux or has_flux_dens:
            if not any([has_freq, has_band, has_ener]):
                err_str = ("Has `{}` or `{}`".format(self._KEYS.FLUX,
                                                     self._KEYS.FLUX_DENSITY) +
                           " but None of `{}`, `{}`, `{}`".format(
                               self._KEYS.FREQUENCY, self._KEYS.BAND,
                               self._KEYS.ENERGY))
            elif has_flux and not has_u_flux:
                err_str = "`{}` provided without `{}`.".format(
                    self._KEYS.FLUX, self._KEYS.U_FLUX)
            elif has_flux_dens and not has_u_flux_dens:
                err_str = "`{}` provided without `{}`.".format(
                    self._KEYS.FLUX_DENSITY, self._KEYS.U_FLUX_DENSITY)
            elif has_freq and not has_u_freq:
                err_str = "`{}` provided without `{}`.".format(
                    self._KEYS.FREQUENCY, self._KEYS.U_FREQUENCY)
            elif has_ener and not has_u_ener:
                err_str = "`{}` provided without `{}`.".format(
                    self._KEYS.ENERGY, self._KEYS.U_ENERGY)

        if err_str is not None:
            raise ValueError(err_str)

        return

    def _clean_value_for_key(self, key, value):
        value = super(Photometry, self)._clean_value_for_key(key, value)

        # Do some basic homogenization
        if key == self._KEYS.BAND:
            return bandrepf(value)
        elif key == self._KEYS.INSTRUMENT:
            return instrumentrepf(value)

        return value

    def sort_func(self, key):
        """Specify order for attributes."""
        if key == self._KEYS.TIME:
            return 'aaa'
        if key == self._KEYS.MODEL:
            return 'zzy'
        if key == self._KEYS.SOURCE:
            return 'zzz'
        return key


BAND_REPS = {
    'Ks': ['K_s'],
    'M2': ['uvm2', 'UVM2', 'UVm2', 'Um2', 'm2', 'um2'],
    'W1': ['uvw1', 'UVW1', 'UVw1', 'Uw1', 'w1', 'uw1'],
    'W2': ['uvw2', 'UVW2', 'UVw2', 'Uw2', 'w2', 'uw2']
}

# Some bands are uniquely tied to an instrument/telescope/system, add this
# info here.
BAND_META = {
    'M2': {
        PHOTOMETRY.TELESCOPE: 'Swift',
        PHOTOMETRY.INSTRUMENT: 'UVOT'
    },
    'W1': {
        PHOTOMETRY.TELESCOPE: 'Swift',
        PHOTOMETRY.INSTRUMENT: 'UVOT'
    },
    'W2': {
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
    "K", "C", "CR", "CV", "M2", "W1", "W2", "pg", "Mp", "w", "y", "Z", "F110W",
    "F775W", "F850LP", "VM", "RM", "Ks", "Ic", "Rc"
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
    "M2": 260.,
    "W1": 224.6,
    "W2": 192.8,
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


def instrumentrepf(code):
    for rep in INSTRUMENT_REPS:
        if code in INSTRUMENT_REPS[rep]:
            return rep
    return code


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
    return None


def get_ul_mag(ec, zp=DEFAULT_ZP, sig=DEFAULT_UL_SIGMA):
    dec = Decimal(str(ec))
    dzp = Decimal(str(zp))
    dsig = Decimal(str(sig))
    mag = str(dzp - (D25 * (dsig * dec).log10()))
    emag = str(dec)
    return mag, emag


def set_pd_mag_from_counts(photodict,
                           c,
                           ec='',
                           lec='',
                           uec='',
                           zp=DEFAULT_ZP,
                           sig=DEFAULT_UL_SIGMA):
    with localcontext() as ctx:
        if lec == '' or uec == '':
            lec = ec
            uec = ec
        prec = max(
            get_sig_digits(str(c)),
            get_sig_digits(str(lec)), get_sig_digits(str(uec)))
        ctx.prec = prec
        dlec = Decimal(str(lec))
        duec = Decimal(str(uec))
        dc = Decimal(str(c))
        dzp = Decimal(str(zp))
        dsig = Decimal(str(sig))
        photodict[PHOTOMETRY.ZERO_POINT] = str(zp)
        if float(c) < DEFAULT_UL_SIGMA * float(uec):
            photodict[PHOTOMETRY.UPPER_LIMIT] = True
            photodict[PHOTOMETRY.UPPER_LIMIT_SIGMA] = str(sig)
            photodict[PHOTOMETRY.MAGNITUDE] = str(dzp - (D25 * (dsig * duec
                                                                ).log10()))
            dnec = Decimal('10.0')**(
                (dzp - Decimal(photodict[PHOTOMETRY.MAGNITUDE])) / D25)
            photodict[PHOTOMETRY.E_UPPER_MAGNITUDE] = str(D25 * (
                (dnec + duec).log10() - dnec.log10()))
        else:
            photodict[PHOTOMETRY.MAGNITUDE] = str(dzp - D25 * dc.log10())
            photodict[PHOTOMETRY.E_UPPER_MAGNITUDE] = str(D25 * (
                (dc + duec).log10() - dc.log10()))
            photodict[PHOTOMETRY.E_LOWER_MAGNITUDE] = str(D25 * (
                dc.log10() - (dc - dlec).log10()))
