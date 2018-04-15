"""Class for representing correlations."""
from astrocats.catalog.catdict import CatDict, CatDictError
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection


class CORRELATION(KeyCollection):
    """`KeyCollection` subclass for the the `Correlation` class.

    Attributes
    ----------
    VALUE : ANY
        The actual value being stored.
    KIND : STRING

    """

    # Any
    VALUE = Key('value', priority=10)
    # Numeric
    # Booleans
    DERIVED = Key('derived', KEY_TYPES.BOOL)
    # Strings
    KIND = Key('kind', KEY_TYPES.STRING, listable=True)
    QUANTITY = Key('quantity', KEY_TYPES.STRING)


class Correlation(CatDict):
    """Class to store correlation of a `Quantity` with another `Quantity`."""

    _KEYS = CORRELATION

    def __init__(self, parent, **kwargs):
        """Initialize `Quantity` object."""
        self._REQ_KEY_SETS = [[CORRELATION.VALUE], [CORRELATION.QUANTITY]]

        super(Correlation, self).__init__(parent, **kwargs)

        # Check that value exists
        if (not self[CORRELATION.VALUE] or self[CORRELATION.VALUE] == '--' or
                self[CORRELATION.VALUE] == '-'):
            raise CatDictError(
                "Value '{}' is empty, not adding to '{}'".format(self[
                    CORRELATION.VALUE], parent[parent._KEYS.NAME]))
