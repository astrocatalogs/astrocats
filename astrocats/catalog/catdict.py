"""
"""
from collections import OrderedDict

from .key import KEY_TYPES, KeyCollection


class CatDict(OrderedDict):
    """
    """

    # The `KeyCollection` object associated with this dictionary
    #    For example, `PHOTOMETRY` for the `Photometry(CatDict)` subclass
    _KEYS = KeyCollection

    # If `_ALLOW_UNKNOWN_KEYS` is 'False', then only parameters with names
    #    included in `_KEYS` are allowed.  Others will raise an error.
    #    If this parameter is 'True', then parameters corresponding to those in
    #    `_KEYS` are still checked (for type etc), but additional parameters
    #    are just tacked onto the `CatDict` object without any checks or errors.
    _ALLOW_UNKNOWN_KEYS = True

    REQ_KEY_TYPES = []

    def __init__(self, **kwargs):
        # Iterate over all `_KEYS` parameters, load each if given
        #    note that the stored 'values' are the `Key` objects, referred to
        #    here with the name 'key'
        for key in self._KEYS.vals():
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
                self[key] = self._clean_value_for_key(key, value)

        # If we require all parameters to be a key in `PHOTOMETRY`, then all
        #    elements should have been removed from `kwargs`
        if not self._ALLOW_UNKNOWN_KEYS and len(kwargs):
            raise ValueError(
                "All permitted keys stored, remaining: '{}'".format(kwargs))

        # Make sure that currently stored values are valid
        self._check()

        return

    def __repr__(self):
        pass

    def is_duplicate_of(self, other):
        # If these are not the same type, return False
        if type(other) is not type(self):
            return False

        # Go over all expected parameters and check equality of each
        for key in self._KEYS.vals():
            # Skip parameters which shouldnt be compared
            if not key.compare:
                continue

            # If only one object has this parameter, not the same
            if (key in self) != (key in other):
                return False
            # If self doesnt have this parameter (and thus neither does), skip
            if key not in self:
                continue

            # Now, both objects have the same parameter, compare them
            if self[key] != other[key]:
                return False

        return True

    def _check(self):
        for req_any in self.REQ_KEY_TYPES:
            if not any([req_key in self for req_key in req_any]):
                err_str = "Require one or more of: " + ",".join(
                    "'{}'".format(rk) for rk in req_any)
                raise ValueError(err_str)

        return

    def _clean_value_for_key(self, key, value):
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
