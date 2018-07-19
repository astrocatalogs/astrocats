"""
"""
import os
import sys

from astrocats.catalog.catdict import CatDictError
from astrocats.catalog.key import KEY_TYPES, Key
from astrocats.catalog.utils import is_number, uniq_cdl

_PAS_PATH = "/Users/lzkelley/Research/catalogs/astroschema"
if _PAS_PATH not in sys.path:
    sys.path.append(_PAS_PATH)

import pyastroschema as pas  # noqa

SCHEMA_DIR = "/Users/lzkelley/Research/catalogs/astrocats/astrocats/schema/"


class My_Meta_Struct(pas.struct.Meta_Struct):

    def __init__(self, parent, key=None, **kwargs):
        # super(Meta_Struct, self).__init__(self._SCHEMA_NAME, *args, **kwargs)
        super(My_Meta_Struct, self).__init__(extendable=True, **kwargs)
        self._key = key
        self._parent = parent
        return

    def append_sources_from(self, other):
        """Merge the source alias lists of two CatDicts."""
        # Get aliases lists from this `CatDict` and other
        self_aliases = self[self._KEYS.SOURCE].split(',')
        other_aliases = other[self._KEYS.SOURCE].split(',')

        # Store alias to `self`
        self[self._KEYS.SOURCE] = uniq_cdl(self_aliases + other_aliases)

        return


class _Source(My_Meta_Struct):

    _SCHEMA_NAME = os.path.join(SCHEMA_DIR, "source.json")
    # _SCHEMA_NAME = "source"

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


class _Quantity(My_Meta_Struct):

    _SCHEMA_NAME = os.path.join(SCHEMA_DIR, "quantity.json")
    # _SCHEMA_NAME = "quantity"

    def __init__(self, parent, key=None, **kwargs):
        super(_Quantity, self).__init__(parent, key=key, **kwargs)

        # Aliases not added if in DISTINCT_FROM
        if self._key == parent._KEYS.ALIAS:
            value = parent.catalog.clean_entry_name(self[self._KEYS.VALUE])
            for df in parent.get(parent._KEYS.DISTINCT_FROM, []):
                if value in parent.catalog.entries:
                    aliases = parent.catalog.entries[value].get_aliases()
                else:
                    aliases = [value]

                if df[self._KEYS.VALUE] in aliases:
                    raise CatDictError("Alias '{}' in '{}'\' '{}' list".format(
                        value, parent[parent._KEYS.NAME], parent._KEYS.DISTINCT_FROM))

            if value in parent.catalog.entries:
                for df in parent.catalog.entries[value].get(
                        parent._KEYS.DISTINCT_FROM, []):
                    if df[self._KEYS.VALUE] in parent.get_aliases():
                        err = "Alias '{}' in '{}'\' '{}' list".format(
                            value, parent[parent._KEYS.NAME], parent._KEYS.DISTINCT_FROM)
                        raise CatDictError(err)

        # Check that value exists
        if (not self[self._KEYS.VALUE] or self[self._KEYS.VALUE] == '--' or
                self[self._KEYS.VALUE] == '-'):
            err = "Value '{}' is empty, not adding to '{}'".format(
                self[self._KEYS.VALUE], parent[parent._KEYS.NAME])
            raise CatDictError(err)

        if not parent._clean_quantity(self):
            err = "Value '{}' did not survive cleaning process, not adding to '{}'.".format(
                self[self._KEYS.VALUE], parent[parent._KEYS.NAME])
            raise CatDictError(err)

        # Check that quantity value matches type after cleaning
        if isinstance(self._key, Key):
            if (self._key.type == KEY_TYPES.NUMERIC) and (not is_number(self[self._KEYS.VALUE])):
                err = "Value '{}' is not numeric, not adding to '{}'".format(
                    self[self._KEYS.VALUE], parent[parent._KEYS.NAME])
                raise CatDictError(err)

        # print("re-Validating...")
        self.validate()
        return

    def sort_func(self, key):
        """Sorting logic for `Quantity` objects."""
        if key == self._KEYS.VALUE:
            return 'aaa'
        if key == self._KEYS.SOURCE:
            return 'zzz'
        return key
