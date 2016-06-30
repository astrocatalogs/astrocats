"""Class for representing sources of data.
"""
from astrocats.catalog.catdict import CatDict
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection

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


class Source(CatDict):
    """
    """

    _KEYS = SOURCE

    def __init__(self, **kwargs):
        super().__init__(kwargs)
        self.REQ_KEY_TYPES = [
            [SOURCE.BIBCODE, SOURCE.URL, SOURCE.NAME]
        ]
