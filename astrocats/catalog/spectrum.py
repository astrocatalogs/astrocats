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
    REDSHIFT = Key('redshift', KEY_TYPES.NUMERIC, compare=False)
    AIRMASS = Key('airmass', KEY_TYPES.NUMERIC, compare=False)

    FILENAME = Key('filename', KEY_TYPES.STRING)
    U_FLUXES = Key('u_fluxes', KEY_TYPES.STRING, compare=False)
    U_ERRORS = Key('u_errors', KEY_TYPES.STRING, compare=False)
    INSTRUMENT = Key('instrument', KEY_TYPES.STRING, compare=False)
    OBSERVATORY = Key('observatory', KEY_TYPES.STRING, compare=False)
    OBSERVER = Key('observer', KEY_TYPES.STRING, compare=False)
    MODEL = Key('model', KEY_TYPES.STRING, compare=False)
    REALIZATION = Key('realization', KEY_TYPES.STRING, priority=15)
    SOURCE = Key('source', KEY_TYPES.STRING, compare=False)
    REDUCER = Key('reducer', KEY_TYPES.STRING, compare=False)
    REDUCTION = Key('reduction', KEY_TYPES.STRING, compare=False)
    SURVEY = Key('survey', KEY_TYPES.STRING, compare=False)
    TELESCOPE = Key('telescope', KEY_TYPES.STRING, compare=False)
    U_TIME = Key('u_time', KEY_TYPES.STRING, compare=False)
    U_WAVELENGTHS = Key('u_wavelengths', KEY_TYPES.STRING, compare=False)
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
            [SPECTRUM.U_FLUXES, SPECTRUM.FILENAME],
            [SPECTRUM.U_WAVELENGTHS, SPECTRUM.FILENAME],
        ]

        # FIX: add this back in
        # [SPECTRUM.TIME, SPECTRUM.HOST]

        # Note: `_check()` is called at end of `super().__init__`
        super(Spectrum, self).__init__(parent, **kwargs)

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
                if SPECTRUM.U_ERRORS not in self:
                    raise ValueError(
                        "Without `{}`,".format(SPECTRUM.DATA) +
                        " but with `{}`,".format(SPECTRUM.ERRORS) +
                        " `{}` also required".format(SPECTRUM.U_ERRORS))
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

        if self._KEYS.U_TIME not in self:
            self._log.info('`{}` not found in spectrum, assuming '
                           ' MJD.'.format(self._KEYS.U_TIME))
            self[self._KEYS.U_TIME] = 'MJD'

        return

    def _check(self):
        """

        """
        # Run the super method
        super(Spectrum, self)._check()

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

    def is_duplicate_of(self, other):
        if super(Spectrum, self).is_duplicate_of(other):
            return True
        row_matches = 0
        for ri, row in enumerate(self.get(self._KEYS.DATA, [])):
            lambda1, flux1 = tuple(row[0:2])
            if (self._KEYS.DATA not in other or
                    ri > len(other[self._KEYS.DATA])):
                break
            lambda2, flux2 = tuple(other[self._KEYS.DATA][ri][0:2])
            minlambdalen = min(len(lambda1), len(lambda2))
            minfluxlen = min(len(flux1), len(flux2))
            if (lambda1[:minlambdalen + 1] == lambda2[:minlambdalen + 1] and
                    flux1[:minfluxlen + 1] == flux2[:minfluxlen + 1] and
                    float(flux1[:minfluxlen + 1]) != 0.0):
                row_matches += 1
            # Five row matches should be enough to be sure spectrum is a dupe.
            if row_matches >= 5:
                return True
            # Matches need to happen in the first 10 rows.
            if ri >= 10:
                break
        return False

    def sort_func(self, key):
        if key == self._KEYS.TIME:
            return 'aaa'
        if key == self._KEYS.DATA:
            return 'zzy'
        if key == self._KEYS.SOURCE:
            return 'zzz'
        return key
