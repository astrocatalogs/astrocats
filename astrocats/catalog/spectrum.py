"""Class for representing spectra.
"""
from Collections import OrderedDict

from .key import KEY_TYPES, Key, KeyCollection

# If `REQUIRE_KEY_IN_PHOTOMETRY` is 'True', then only parameters with names
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


class Spectrum(OrderedDict):
    """
    """

    def __init__(self, **kwargs):
        # Iterate over all `SPECTRA` parameters, load each if given note that
        # the stored `values` are the `Key` objects, referred to here with the
        # name 'key'
        for key in SPECTRA.vals():
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
        # elements should have been removed from `kwargs`
        if REQUIRE_KEY_IN_SPECTRA and len(kwargs):
            raise ValueError(
                "All permitted keys stored, remaining: '{}'".format(kwargs))

        # Make sure that currently stored values are valid
        self._check()

        return

    def __repr__(self):
        pass

    def _check(self):
        REQ_KEY_TYPES = [
            [SPECTRA.SOURCE],
            [SPECTRA.TIME, SPECTRA.HOST]
        ]

        for req_any in REQ_KEY_TYPES:
            if not any([req_key in self for req_key in req_any]):
                err_str = "Require one of: " + ",".join(
                    "'{}'".format(rk) for rk in req_any)
                raise ValueError(err_str)

        return
