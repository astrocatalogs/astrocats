"""Class for representing spectra.
"""
from astrocats.catalog.catdict import CatDict
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection
from astrocats.catalog.utils import trim_str_arr


class SPECTRUM(KeyCollection):
    # Arrays
    DATA = Key('data', compare=False)
    ERRORS = Key('errors', compare=False)
    EXCLUDE = Key('exclude', compare=False)
    WAVELENGTHS = Key('wavelengths', compare=False)
    FLUXES = Key('fluxes', compare=False)

    # Numbers
    E_LOWER_TIME = Key('e_lower_time', KEY_TYPES.NUMERIC, compare=False)
    E_TIME = Key('e_time', KEY_TYPES.NUMERIC, compare=False)
    E_UPPER_TIME = Key('e_upper_time', KEY_TYPES.NUMERIC, compare=False)
    SNR = Key('snr', KEY_TYPES.NUMERIC, compare=False)
    TIME = Key('time', KEY_TYPES.NUMERIC, compare=False, listable=True)

    FILENAME = Key('filename', KEY_TYPES.STRING)
    FLUX_UNIT = Key('fluxunit', KEY_TYPES.STRING, compare=False)
    ERROR_UNIT = Key('fluxunit', KEY_TYPES.STRING, compare=False)
    INSTRUMENT = Key('instrument', KEY_TYPES.STRING, compare=False)
    OBSERVATORY = Key('observatory', KEY_TYPES.STRING, compare=False)
    OBSERVER = Key('observer', KEY_TYPES.STRING, compare=False)
    SOURCE = Key('source', KEY_TYPES.STRING, compare=False)
    SURVEY = Key('survey', KEY_TYPES.STRING, compare=False)
    TELESCOPE = Key('telescope', KEY_TYPES.STRING, compare=False)
    U_TIME = Key('u_time', KEY_TYPES.STRING, compare=False)
    WAVE_UNIT = Key('waveunit', KEY_TYPES.STRING, compare=False)
    # Booleans
    DEREDDENED = Key('dereddened', KEY_TYPES.BOOL, compare=False)
    DEREDSHIFTED = Key('deredshifted', KEY_TYPES.BOOL, compare=False)
    HOST = Key('host', KEY_TYPES.BOOL, compare=False)
    INCLUDES_HOST = Key('includeshost', KEY_TYPES.BOOL, compare=False)


class Spectrum(CatDict):
    """Class for storing a single spectrum associated with an astrophysical
    entity.
    """

    _KEYS = SPECTRUM

    def __init__(self, parent, **kwargs):
        self._REQ_KEY_SETS = [
            [SPECTRUM.SOURCE, SPECTRUM.FILENAME],
            [SPECTRUM.FLUX_UNIT, SPECTRUM.FILENAME],
            [SPECTRUM.WAVE_UNIT, SPECTRUM.FILENAME],
        ]

        # FIX: add this back in
        # [SPECTRUM.TIME, SPECTRUM.HOST]

        # Note: `_check()` is called at end of `super().__init__`
        super().__init__(parent, **kwargs)

        # If `data` is not given, construct it from wavelengths, fluxes
        # [errors] `errors` is optional, but if given, then `errorunit` is also
        # req'd
        if SPECTRUM.DATA not in self:
            try:
                wavelengths = self[SPECTRUM.WAVELENGTHS]
                fluxes = self[SPECTRUM.FLUXES]
            except KeyError:
                if SPECTRUM.FILENAME in self:
                    return
                else:
                    err_str = "Neither data nor (wavelengths and fluxes) given"
                    self._log.error(err_str)
                    raise

            errors = self.get(SPECTRUM.ERRORS, None)
            if (errors is not None and
                    max([float(err) for err in errors]) > 0.0):
                if SPECTRUM.ERROR_UNIT not in self:
                    raise ValueError(
                        "Without `{}`,".format(SPECTRUM.DATA) +
                        " but with `{}`,".format(SPECTRUM.ERRORS) +
                        " `{}` also required".format(SPECTRUM.ERROR_UNIT))
                data = [trim_str_arr(wavelengths), trim_str_arr(fluxes),
                        trim_str_arr(errors)]
            else:
                data = [trim_str_arr(wavelengths), trim_str_arr(fluxes)]

            self[SPECTRUM.DATA] = [list(i) for i in zip(*data)]
            if SPECTRUM.WAVELENGTHS in self:
                del self[SPECTRUM.WAVELENGTHS]
            if SPECTRUM.FLUXES in self:
                del self[SPECTRUM.FLUXES]
            if SPECTRUM.ERRORS in self:
                del self[SPECTRUM.ERRORS]

        return

    def _check(self):
        """

        """
        # Run the super method
        super()._check()

        err_str = None
        has_data = self._KEYS.DATA in self
        has_wave = self._KEYS.WAVELENGTHS in self
        has_flux = self._KEYS.FLUXES in self
        has_filename = self._KEYS.FILENAME in self

        if not has_data:
            if (not has_wave or not has_flux) and not has_filename:
                err_str = (
                    "If `{}` not given".format(self._KEYS.DATA) +
                    "; `{}` or `{}` needed".format(
                        self._KEYS.WAVELENGTHS, self._KEYS.FLUXES))

        if err_str is not None:
            raise ValueError(err_str)

        return

    def sort_func(self, key):
        if key == self._KEYS.TIME:
            return 'aaa'
        if key == self._KEYS.DATA:
            return 'zzy'
        if key == self._KEYS.SOURCE:
            return 'zzz'
        return key
