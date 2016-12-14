"""Handle known errors in catalog data.
"""
from astrocats.catalog.catdict import CatDict
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection


class ERROR(KeyCollection):
    # Any
    VALUE = Key('value')
    EXTRA = Key('extra')
    # Numeric
    ERROR = Key('error', KEY_TYPES.NUMERIC)
    # Booleans
    UPPER_LIMIT = Key('upperlimit', KEY_TYPES.BOOL)
    # Strings
    UNIT = Key('unit', KEY_TYPES.STRING)
    KIND = Key('kind', KEY_TYPES.STRING)
    SOURCE = Key('source', KEY_TYPES.STRING)
    SOURCE_KIND = Key('sourcekind', KEY_TYPES.STRING)


class Error(CatDict):
    """Class to represent (known) errors in catalog source data.
    """

    _KEYS = ERROR

    def __init__(self, parent, **kwargs):
        self._REQ_KEY_SETS = [
            [ERROR.VALUE]
        ]
        super(Error, self).__init__(parent, **kwargs)
