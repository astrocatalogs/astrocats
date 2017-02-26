"""Class for representing spectra.
"""
from astrocats.catalog.catdict import CatDict, CatDictError
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection
from astrocats.catalog.utils import is_number


class QUANTITY(KeyCollection):
    """`KeyCollection` subclass for the the `Quantity` class.

    Attributes
    ----------
    VALUE : ANY
        The actual value being stored.
    E_VALUE : NUMERIC
        The 'error' (uncertainty) associate with the measured value.
    E_LOWER_VALUE : NUMERIC
        The 'error' (uncertainty) in the lower (more negative) direction.
    E_UPPER_VALUE : NUMERIC
        The 'error' (uncertainty) in the upper (more positive) direction.
    PROB : NUMERIC
    UPPER_LIMIT : BOOL
        If this value corresponds to a measured upper-limit.
    DESCRIPTION : STRING
        Verbal description or notes on a quantity.
    U_VALUE : STRING
        Unit of measurement associated with this value.
    KIND : STRING
    SOURCE : STRING
        The alias number(s) of the 'source(s)' (reference) for this value.
        NOTE: while source-aliases are integers referring to the actual
        source-entries, this entry is stored as a *single* string, where
        numerous souces are stored as a single, comma-separated string.
        e.g. ``Entry[QUANTITY.SOURCE] = '1, 4, 5'``.

    """
    # Any
    VALUE = Key('value', priority=10)
    # Numeric
    E_VALUE = Key('e_value', KEY_TYPES.NUMERIC)
    E_LOWER_VALUE = Key('e_lower_value', KEY_TYPES.NUMERIC)
    E_UPPER_VALUE = Key('e_upper_value', KEY_TYPES.NUMERIC)
    PROB = Key('probability', KEY_TYPES.NUMERIC)
    # Booleans
    UPPER_LIMIT = Key('upperlimit', KEY_TYPES.BOOL)
    LOWER_LIMIT = Key('lowerlimit', KEY_TYPES.BOOL)
    DERIVED = Key('derived', KEY_TYPES.BOOL)
    # Strings
    DESCRIPTION = Key('description', KEY_TYPES.STRING, compare=False)
    U_VALUE = Key('u_value', KEY_TYPES.STRING)
    KIND = Key('kind', KEY_TYPES.STRING, listable=True)
    SOURCE = Key('source', KEY_TYPES.STRING, compare=False)
    MODEL = Key('model', KEY_TYPES.STRING, compare=False)
    REALIZATION = Key('realization', KEY_TYPES.STRING)


class Quantity(CatDict):
    """Class to store a single item of data associated with an `Entry`.
    """

    _KEYS = QUANTITY

    def __init__(self, parent, **kwargs):
        self._REQ_KEY_SETS = [[QUANTITY.VALUE], [QUANTITY.SOURCE]]

        super(Quantity, self).__init__(parent, **kwargs)

        # Aliases not added if in DISTINCT_FROM
        if self._key == parent._KEYS.ALIAS:
            value = parent.catalog.clean_entry_name(self[QUANTITY.VALUE])
            for df in parent.get(parent._KEYS.DISTINCT_FROM, []):
                aliases = parent.catalog.entries[value].get_aliases(
                ) if value in parent.catalog.entries else [value]
                if df[QUANTITY.VALUE] in aliases:
                    raise CatDictError("Alias '{}' in '{}'\' '{}' list".format(
                        value, parent[
                            parent._KEYS.NAME], parent._KEYS.DISTINCT_FROM))
            if value in parent.catalog.entries:
                for df in parent.catalog.entries[value].get(
                        parent._KEYS.DISTINCT_FROM, []):
                    if df[QUANTITY.VALUE] in parent.get_aliases():
                        raise CatDictError(
                            "Alias '{}' in '{}'\' '{}' list".format(
                                value, parent[parent._KEYS.NAME],
                                parent._KEYS.DISTINCT_FROM))

        # Check that value exists
        if (not self[QUANTITY.VALUE] or self[QUANTITY.VALUE] == '--' or
                self[QUANTITY.VALUE] == '-'):
            raise CatDictError(
                "Value '{}' is empty, not adding to '{}'".format(self[
                    QUANTITY.VALUE], parent[parent._KEYS.NAME]))

        if not parent._clean_quantity(self):
            raise CatDictError(
                "Value '{}' did not survive cleaning process, not adding to "
                " '{}'.".format(self[QUANTITY.VALUE], parent[
                    parent._KEYS.NAME]))

        # Check that quantity value matches type after cleaning
        if (isinstance(self._key, Key) and
                self._key.type == KEY_TYPES.NUMERIC and
                not is_number(self[QUANTITY.VALUE])):
            raise CatDictError(
                "Value '{}' is not numeric, not adding to '{}'".format(self[
                    QUANTITY.VALUE], parent[parent._KEYS.NAME]))

    def sort_func(self, key):
        if key == self._KEYS.VALUE:
            return 'aaa'
        if key == self._KEYS.SOURCE:
            return 'zzz'
        return key
