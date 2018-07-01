"""Class for representing sources of data."""
import sys

# from astrocats.catalog.catdict import CatDict
# from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection


_PAS_PATH = "/Users/lzkelley/Research/catalogs/astroschema"
if _PAS_PATH not in sys.path:
    sys.path.append(_PAS_PATH)

import pyastroschema as pas  # noqa

SOURCE = pas.source.Source.get_keychain()


class Source(pas.source.Source):
    """Representation for the source/attribution of a data element."""

    _KEYS = SOURCE

    def __init__(self, parent, key=None, **kwargs):
        super(Source, self).__init__(extendable=True, **kwargs)
        self._key = key
        self._parent = parent
        # self._log = parent.catalog.log
        return

    def sort_func(self, key):
        if key == self._KEYS.NAME:
            return 'aaa'
        if key == self._KEYS.BIBCODE:
            return 'aab'
        if key == self._KEYS.ARXIVID:
            return 'aac'
        if key == self._KEYS.DOI:
            return 'aad'
        if key == self._KEYS.ALIAS:
            return 'zzz'
        return key

    def append_sources_from(self, other):
        """`CatDict.append_sources_from` should never be called in `Source`.
        """
        raise RuntimeError("`Source.append_sources_from` called.")

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
        except Exception:
            return None
