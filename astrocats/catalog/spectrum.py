"""Class for representing spectra.
"""
from astrocats.catalog.catdict import CatDict
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection

# If `_ALLOW_UNKNOWN_KEYS` is 'True', then only parameters with names
#    included in `PHOTOMETRY` are allowed.  Others will raise an error.
#    If this parameter is 'False', then parameters corresponding to those in
#    `PHOTOMETRY` are still checked (for type etc), but additional parameters
#    are just tacked onto the `Photometry` object without any checks or errors.
REQUIRE_KEY_IN_SPECTRA = False


class SPECTRUM(KeyCollection):
    # Arrays
    DATA = Key('data', compare=False)
    EXCLUDE = Key('exclude', compare=False)
    # Numbers
    E_LOWER_TIME = Key('e_lower_time', KEY_TYPES.NUMERIC, compare=False)
    E_TIME = Key('e_time', KEY_TYPES.NUMERIC, compare=False)
    E_UPPER_TIME = Key('e_upper_time', KEY_TYPES.NUMERIC, compare=False)
    SNR = Key('snr', KEY_TYPES.NUMERIC, compare=False)
    TIME = Key('time', KEY_TYPES.NUMERIC, listable=True)

    FILENAME = Key('filename', KEY_TYPES.STRING)
    FLUXUNIT = Key('fluxunit', KEY_TYPES.STRING, compare=False)
    INSTRUMENT = Key('instrument', KEY_TYPES.STRING, compare=False)
    OBSERVATORY = Key('observatory', KEY_TYPES.STRING, compare=False)
    OBSERVER = Key('observer', KEY_TYPES.STRING, compare=False)
    SOURCE = Key('source', KEY_TYPES.STRING, compare=False)
    SURVEY = Key('survey', KEY_TYPES.STRING, compare=False)
    TELESCOPE = Key('telescope', KEY_TYPES.STRING, compare=False)
    U_TIME = Key('u_time', KEY_TYPES.STRING, compare=False)
    WAVEUNIT = Key('waveunit', KEY_TYPES.STRING, compare=False)
    # Booleans
    DEREDDENED = Key('dereddened', KEY_TYPES.BOOL, compare=False)
    DEREDSHIFTED = Key('deredshifted', KEY_TYPES.BOOL, compare=False)
    HOST = Key('host', KEY_TYPES.BOOL, compare=False)
    INCLUDES_HOST = Key('includeshost', KEY_TYPES.BOOL, compare=False)


class Spectrum(CatDict):
    """
    """

    _KEYS = SPECTRUM

    def __init__(self, **kwargs):
        super().__init__(kwargs)
        self.REQ_KEY_TYPES = [
            [SPECTRUM.SOURCE],
            [SPECTRUM.TIME, SPECTRUM.HOST]
        ]
