"""Class for representing sources of data.
"""
from astrocats.catalog.catdict import CatDict
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection


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
            [SOURCE.BIBCODE, SOURCE.URL, SOURCE.NAME],
            [SOURCE.ALIAS]
        ]
        return

    def append_sources_from(self, other):
        """`CatDict.append_sources_from` should never be called in `Source`.
        """
        raise RuntimeError("`Source.append_sources_from` called.")
