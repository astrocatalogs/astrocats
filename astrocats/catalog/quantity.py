"""Class for representing spectra.
"""
from astrocats.catalog.catdict import CatDict, CatDictError
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection


class QUANTITY(KeyCollection):
    # Any
    VALUE = Key('value')
    # Numeric
    ERROR = Key('error', KEY_TYPES.NUMERIC)
    PROB = Key('probability', KEY_TYPES.NUMERIC)
    # Booleans
    UPPER_LIMIT = Key('upperlimit', KEY_TYPES.BOOL)
    # Strings
    UNIT = Key('unit', KEY_TYPES.STRING)
    KIND = Key('kind', KEY_TYPES.STRING)
    SOURCE = Key('source', KEY_TYPES.STRING, compare=False)


class Quantity(CatDict):
    """Class for storing a single item of data associated with an astrophysical
    entity.
    """

    _KEYS = QUANTITY

    def __init__(self, parent, **kwargs):
        self._REQ_KEY_SETS = [
            [QUANTITY.VALUE],
            [QUANTITY.SOURCE]
        ]

        super().__init__(parent, **kwargs)

        # Aliases not added if in DISTINCT_FROM
        if self._key == parent._KEYS.ALIAS:
            value = parent.clean_entry_name(self[QUANTITY.VALUE])
            for df in parent.get(parent._KEYS.DISTINCT_FROM, []):
                if value == df[QUANTITY.VALUE]:
                    raise CatDictError(
                        "Alias '{}' in '{}'\' '{}' list".format(
                            value, parent[parent._KEYS.NAME],
                            parent._KEYS.DISTINCT_FROM))

        # Check that value exists
        if (not self[QUANTITY.VALUE] or self[QUANTITY.VALUE] == '--' or
                self[QUANTITY.VALUE] == '-'):
            raise CatDictError(
                "Value '{}' is empty, not adding to '{}'".format(
                    self[QUANTITY.VALUE], parent[parent._KEYS.NAME]))

        parent._clean_quantity(self)

    def sort_func(self, key):
        if key == self._KEYS.VALUE:
            return 'aaa'
        if key == self._KEYS.SOURCE:
            return 'zzz'
        return key
