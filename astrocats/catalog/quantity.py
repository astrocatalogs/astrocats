"""Class for representing quantities."""
from astrocats.catalog.catdict import CatDictError
from astrocats.catalog.key import KEY_TYPES, Key
from astrocats.catalog.utils import is_number

import pyastroschema as pas  # noqa

QUANTITY = pas.struct.Quantity.get_keychain()


class Quantity(pas.struct.Quantity):
    """Class to store a single item of data associated with an `Entry`."""

    _KEYS = QUANTITY

    def __init__(self, parent, key=None, **kwargs):
        super(Quantity, self).__init__(extendable=True, **kwargs)
        # super(Quantity, self).__init__(parent, **kwargs)
        self._key = key
        self._parent = parent

        # Aliases not added if in DISTINCT_FROM
        if self._key == parent._KEYS.ALIAS:
            value = parent.catalog.clean_entry_name(self[QUANTITY.VALUE])
            for df in parent.get(parent._KEYS.DISTINCT_FROM, []):
                if value in parent.catalog.entries:
                    aliases = parent.catalog.entries[value].get_aliases()
                else:
                    aliases = [value]

                if df[QUANTITY.VALUE] in aliases:
                    raise CatDictError("Alias '{}' in '{}'\' '{}' list".format(
                        value, parent[parent._KEYS.NAME], parent._KEYS.DISTINCT_FROM))

            if value in parent.catalog.entries:
                for df in parent.catalog.entries[value].get(
                        parent._KEYS.DISTINCT_FROM, []):
                    if df[QUANTITY.VALUE] in parent.get_aliases():
                        err = "Alias '{}' in '{}'\' '{}' list".format(
                            value, parent[parent._KEYS.NAME], parent._KEYS.DISTINCT_FROM)
                        raise CatDictError(err)

        # Check that value exists
        if (not self[QUANTITY.VALUE] or self[QUANTITY.VALUE] == '--' or
                self[QUANTITY.VALUE] == '-'):
            err = "Value '{}' is empty, not adding to '{}'".format(
                self[QUANTITY.VALUE], parent[parent._KEYS.NAME])
            raise CatDictError(err)

        if not parent._clean_quantity(self):
            err = "Value '{}' did not survive cleaning process, not adding to '{}'.".format(
                self[QUANTITY.VALUE], parent[parent._KEYS.NAME])
            raise CatDictError(err)

        # Check that quantity value matches type after cleaning
        if isinstance(self._key, Key):
            if (self._key.type == KEY_TYPES.NUMERIC) and (not is_number(self[QUANTITY.VALUE])):
                err = "Value '{}' is not numeric, not adding to '{}'".format(
                    self[QUANTITY.VALUE], parent[parent._KEYS.NAME])
                raise CatDictError(err)

    def sort_func(self, key):
        """Sorting logic for `Quantity` objects."""
        if key == self._KEYS.VALUE:
            return 'aaa'
        if key == self._KEYS.SOURCE:
            return 'zzz'
        return key
