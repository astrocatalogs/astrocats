"""Class for representing sources of data.
"""
from astrocats.catalog.catdict import CatDict
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection


class SOURCE(KeyCollection):
    # Strings
    NAME = Key('name', KEY_TYPES.STRING)
    BIBCODE = Key('bibcode', KEY_TYPES.STRING)
    URL = Key('url', KEY_TYPES.STRING, compare=False)
    ACKNOWLEDGMENT = Key('acknowledgment', KEY_TYPES.STRING, compare=False)
    REFERENCE = Key('reference', KEY_TYPES.STRING, compare=False)
    # Numbers
    ALIAS = Key('alias', KEY_TYPES.NUMERIC, compare=False)
    # Booleans
    SECONDARY = Key('secondary', KEY_TYPES.BOOL, compare=False)


class Source(CatDict):
    """
    """

    _KEYS = SOURCE

    def __init__(self, parent, **kwargs):
        self._REQ_KEY_SETS = [
            [SOURCE.ALIAS],
            [SOURCE.BIBCODE, SOURCE.URL, SOURCE.NAME]
        ]
        super().__init__(parent, **kwargs)
        return

    def sort_func(self, key):
        if key == self._KEYS.NAME:
            return 'aaa'
        if key == self._KEYS.BIBCODE:
            return 'aab'
        if key == self._KEYS.ALIAS:
            return 'zzz'
        return key

    def append_sources_from(self, other):
        """`CatDict.append_sources_from` should never be called in `Source`.
        """
        raise RuntimeError("`Source.append_sources_from` called.")

    def is_duplicate_of(self, other):
        """Check if this Source is a duplicate of another.

        Unlike the function in the super class, this method will return True
        if *either* name or bibcode is the same.
        """
        # If these are not the same type, return False
        if type(other) is not type(self):
            return False

        # Go over all expected parameters and check equality of each
        for key in self._KEYS.vals():
            # Skip parameters which shouldnt be compared
            if not key.compare:
                continue

            # If only one object has this parameter, not the same
            # if (key in self) != (key in other):
            #     continue
            # If self doesnt have this parameter (and thus neither does), skip
            if key not in self or key not in other:
                continue

            # Now, both objects have the same parameter, compare them
            if self[key] == other[key]:
                return True

        return False
