"""Class for representing sources of data.
"""
from astrocats.catalog.catdict import CatDict
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection


class SOURCE(KeyCollection):
    """`KeyCollection` for the `Source` class.

    Attributes
    ----------
    NAME : STRING
    BIBCODE : STRING
    URL : STRING
    ACKNOWLEDGMENT : STRING
    REFERENCE : STRING
    ALIAS : NUMERIC
        Numerical alias (shorthand) for this entry.  Saved as a string (or
        list of strings), despite being stored as an integer.
    SECONDARY : BOOL
        Whether the given source is one which collected data from another,
        'Primary'-source, from which it actually originated

    """
    # Strings
    NAME = Key('name', KEY_TYPES.STRING)
    BIBCODE = Key('bibcode', KEY_TYPES.STRING)
    ARXIVID = Key('arxivid', KEY_TYPES.STRING)
    URL = Key('url', KEY_TYPES.STRING, compare=False)
    ACKNOWLEDGMENT = Key('acknowledgment', KEY_TYPES.STRING, compare=False)
    REFERENCE = Key('reference', KEY_TYPES.STRING, compare=False)
    # Numbers
    ALIAS = Key('alias', KEY_TYPES.NUMERIC, compare=False)
    # Booleans
    SECONDARY = Key('secondary', KEY_TYPES.BOOL, compare=False)
    PRIVATE = Key('private', KEY_TYPES.BOOL, compare=False)


class Source(CatDict):
    """Representation for the source/attribution of a data element.
    """

    _KEYS = SOURCE

    def __init__(self, parent, **kwargs):
        self._REQ_KEY_SETS = [
            [SOURCE.ALIAS],
            [SOURCE.BIBCODE, SOURCE.ARXIVID, SOURCE.URL, SOURCE.NAME]
        ]
        super(Source, self).__init__(parent, **kwargs)
        return

    def sort_func(self, key):
        if key == self._KEYS.NAME:
            return 'aaa'
        if key == self._KEYS.BIBCODE:
            return 'aab'
        if key == self._KEYS.ARXIVID:
            return 'aac'
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
        for key in self._KEYS.compare_vals():
            # If only one object has this parameter, not the same
            # This is commented out for sources because two sources are
            # considered the same if they share a name *or* a bibcode
            # if (key in self) != (key in other):
            #     continue
            # If self doesnt have this parameter (and thus neither does), skip
            if key not in self or key not in other:
                continue

            # Now, both objects have the same parameter, compare them
            if self[key] == other[key]:
                return True

        return False

    @classmethod
    def bibcode_from_url(cls, url):
        """Given a URL, try to find the ADS bibcode.

        Currently: only `ads` URLs will work, e.g.

        Returns
        -------
        code : str or 'None'
            The Bibcode if found, otherwise 'None'

        """
        try:
            code = url.split('/abs/')
            code = code[1].strip()
            return code
        except:
            return None
