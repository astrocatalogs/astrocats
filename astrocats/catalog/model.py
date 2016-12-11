"""Class for representing models.
"""

from astrocats.catalog.catdict import CatDict
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection


class MODEL(KeyCollection):
    SOURCE = Key('source', KEY_TYPES.STRING, compare=False)
    CODE = Key('code', KEY_TYPES.STRING)
    NAME = Key('name', KEY_TYPES.STRING)
    MODEL = Key('model', KEY_TYPES.DICT, compare=False)
    PARAMETERS = Key('parameters', KEY_TYPES.DICT, compare=False)
    VERSION = Key('version', KEY_TYPES.STRING)
    DATE = Key('date', KEY_TYPES.STRING)

    DESC = Key('description', KEY_TYPES.STRING, compare=False)


class Model(CatDict):
    """Container for a model with associated metadata.

    `Source` citation required.
    """

    _ALLOW_UNKNOWN_KEYS = True
    _KEYS = MODEL

    def __init__(self, parent, **kwargs):
        self._REQ_KEY_SETS = [
            [MODEL.SOURCE],
            [MODEL.ALIAS],
            [MODEL.NAME, MODEL.CODE]
        ]
        # Note: `_check()` is called at end of `super().__init__`
        super().__init__(parent, **kwargs)

        return

    def _check(self):
        """

        """
        # Run the super method
        super()._check()

        err_str = None

        # Future checks go here

        if err_str is not None:
            raise ValueError(err_str)

        return

    def _clean_value_for_key(self, key, value):
        value = super()._clean_value_for_key(key, value)

        return value

    def sort_func(self, key):
        if key == self._KEYS.TIME:
            return 'aaa'
        if key == self._KEYS.SOURCE:
            return 'zzz'
        return key
