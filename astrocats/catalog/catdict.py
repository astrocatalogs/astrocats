"""
"""
from collections import OrderedDict

from .key import Key, KEY_TYPES, KeyCollection


class CatDict(OrderedDict):
    """
    """

    # The `KeyCollection` object associated with this dictionary
    #    For example, `PHOTOMETRY` for the `Photometry(CatDict)` subclass
    KEYS = KeyCollection

    # If `REQUIRE_KEY_IN_PHOTOMETRY` is 'True', then only parameters with names
    #    included in `PHOTOMETRY` are allowed.  Others will raise an error.
    #    If this parameter is 'False', then parameters corresponding to those in
    #    `PHOTOMETRY` are still checked (for type etc), but additional parameters
    #    are just tacked onto the `Photometry` object without any checks or errors.
    REQUIRE_KEY_IN_PHOTOMETRY = True

    def __init__(self, **kwargs):
        # Iterate over all `PHOTOMETRY` parameters, load each if given
        #    note that the stored `values` are the `Key` objects, referred to
        #    here with the name 'key'
        for key in PHOTOMETRY.vals():
            # If this key is given, process and store it
            if key in kwargs:
                # Make sure value is compatible with the 'Key' specification
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
                self[key] = self._clean_value(key, value)

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
            [PHOTOMETRY.SOURCE],
            [PHOTOMETRY.TIME, PHOTOMETRY.HOST],
            [PHOTOMETRY.MAGNITUDE, PHOTOMETRY.FLUX, PHOTOMETRY.FLUX_DENSITY,
             PHOTOMETRY.COUNTS]]

        for req_any in REQ_KEY_TYPES:
            if not any([req_key in self for req_key in req_any]):
                err_str = "Require one or more of: " + ",".join(
                    "'{}'".format(rk) for rk in req_any)
                raise ValueError(err_str)

        return

    def _clean_value(self, key, value):
        """
        """
        # Store whether given value started as a list or not
        single = True
        # Make everything a list for conversions below
        if isinstance(value, list):
            # But if lists arent allowed, and this is, raise error
            if not key.listable:
                raise ValueError(
                    "`value` '{}' for '{}' shouldnt be a list.".format(
                        value, repr(key)))

            single = False
        else:
            single = True
            value = [value]

        # Store booleans as booleans, make sure each element of list is bool
        if key.type == KEY_TYPES.BOOL:
            if not all(isinstance(val, bool) for val in value):
                raise ValueError(
                    "`value` '{}' for '{}' should be boolean".format(
                        value, repr(key)))
        # Strings and numeric types should be stored as strings
        elif key.type in [KEY_TYPES.STRING, KEY_TYPES.NUMERIC]:
            value = [str(val) for val in value]

        # Convert back to single value, if thats how it started
        if single:
            value = value[0]

        return value
