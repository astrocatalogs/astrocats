"""Class for representing model realizations.
"""

from astrocats.catalog.catdict import CatDict
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection


class REALIZATION(KeyCollection):
    SCORE = Key('score', KEY_TYPES.NUMERIC)
    PARAMETERS = Key('parameters', KEY_TYPES.DICT)
    ALIAS = Key('alias', KEY_TYPES.STRING)


class Realization(CatDict):
    """Container for a single model realization.
    """

    _ALLOW_UNKNOWN_KEYS = True
    _KEYS = REALIZATION

    def __init__(self, parent, **kwargs):
        self._REQ_KEY_SETS = []
        # Note: `_check()` is called at end of `super().__init__`
        super(Realization, self).__init__(parent, **kwargs)

        return

    def _check(self):
        """

        """
        # Run the super method
        super(Realization, self)._check()

        err_str = None

        if err_str is not None:
            raise ValueError(err_str)

        return

    def _clean_value_for_key(self, key, value):
        value = super(Realization, self)._clean_value_for_key(key, value)

        return value

    def sort_func(self, key):
        return key
