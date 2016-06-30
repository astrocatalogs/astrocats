"""Class for representing sources of data.
"""
from Collections import OrderedDict

from .key import KEY_TYPES, Key, KeyCollection

# If `REQUIRE_KEY_IN_SOURCE` is 'True', then only parameters with names
# included in `SOURCE` are allowed.  Others will raise an error. If this
# parameter is 'False', then parameters corresponding to those in `SOURCE` are
# still checked (for type etc), but additional parameters are just tacked onto
# the `Source` object without any checks or errors.
REQUIRE_KEY_IN_SOURCE = True


class SOURCE(KeyCollection):
    # Strings
    NAME = Key('name', KEY_TYPES.STRING)
    BIBCODE = Key('bibcode', KEY_TYPES.STRING)
    URL = Key('url', KEY_TYPES.STRING)
    ACKNOWLEDGMENT = Key('acknowledgment', KEY_TYPES.STRING)
    # Numbers
    ALIAS = Key('alias', KEY_TYPES.NUMERIC)
    # Booleans
    SECONDARY = Key('secondary', KEY_TYPES.BOOL)


class Source(OrderedDict):
    """
    """

    def __init__(self, **kwargs):
        # Iterate over all `PHOTOMETRY` parameters, load each if given
        #    note that the stored `values` are the `Key` objects, referred to
        #    here with the name 'key'
        for key in SOURCE.vals():
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
        if REQUIRE_KEY_IN_SOURCE and len(kwargs):
            raise ValueError(
                "All permitted keys stored, remaining: '{}'".format(kwargs))

        # Make sure that currently stored values are valid
        self._check()

        return

    def __repr__(self):
        pass

    def _check(self):
        REQ_KEY_TYPES = [
            [SOURCE.BIBCODE, SOURCE.URL, SOURCE.NAME]]

        for req_any in REQ_KEY_TYPES:
            if not any([req_key in self for req_key in req_any]):
                err_str = "Require one of: " + ",".join(
                    "'{}'".format(rk) for rk in req_any)
                raise ValueError(err_str)

        return
