"""
"""
import os
import sys

import astrocats
from astrocats.catalog import utils

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
        self[self._KEYS.SOURCE] = utils.uniq_cdl(self_aliases + other_aliases)

        return


class _Source(My_Meta_Struct):

    _SCHEMA_NAME = os.path.join(SCHEMA_DIR, "source.json")

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
                    err = "Alias '{}' in '{}'\' '{}' list".format(
                        value, parent[parent._KEYS.NAME], parent._KEYS.DISTINCT_FROM)
                    raise astrocats.catalog.catdict.CatDictError(err)

            if value in parent.catalog.entries:
                for df in parent.catalog.entries[value].get(
                        parent._KEYS.DISTINCT_FROM, []):
                    if df[self._KEYS.VALUE] in parent.get_aliases():
                        err = "Alias '{}' in '{}'\' '{}' list".format(
                            value, parent[parent._KEYS.NAME], parent._KEYS.DISTINCT_FROM)
                        raise astrocats.catalog.catdict.CatDictError(err)

        # Check that value exists
        if (not self[self._KEYS.VALUE] or self[self._KEYS.VALUE] == '--' or
                self[self._KEYS.VALUE] == '-'):
            err = "Value '{}' is empty, not adding to '{}'".format(
                self[self._KEYS.VALUE], parent[parent._KEYS.NAME])
            raise astrocats.catalog.catdict.CatDictError(err)

        if not parent._clean_quantity(self):
            err = "Value '{}' did not survive cleaning process, not adding to '{}'.".format(
                self[self._KEYS.VALUE], parent[parent._KEYS.NAME])
            raise astrocats.catalog.catdict.CatDictError(err)

        # Check that quantity value matches type after cleaning
        '''
        if isinstance(self._key, Key):
            if (self._key.type == KEY_TYPES.NUMERIC) and (not is_number(self[self._KEYS.VALUE])):
                err = "Value '{}' is not numeric, not adding to '{}'".format(
                    self[self._KEYS.VALUE], parent[parent._KEYS.NAME])
                raise astrocats.catalog.catdict.CatDictError(err)
        '''

        self.validate()
        return

    def sort_func(self, key):
        """Sorting logic for `Quantity` objects."""
        if key == self._KEYS.VALUE:
            return 'aaa'
        if key == self._KEYS.SOURCE:
            return 'zzz'
        return key


class _Photometry(My_Meta_Struct):
    """Container for a single photometric point with associated metadata.

    `Source` citation required.
    Photometry can be given as [magnitude, flux, flux-density, counts, luminosity].

    """

    _SCHEMA_NAME = os.path.join(SCHEMA_DIR, "photometry.json")

    def __init__(self, parent, key=None, **kwargs):
        super(_Photometry, self).__init__(parent, key=key, **kwargs)

        # If `BAND` is given, but any of `bandmetaf_keys` is not, try to infer
        if self._KEYS.BAND in self:
            sband = self[self._KEYS.BAND]
            bandmetaf_keys = [self._KEYS.INSTRUMENT, self._KEYS.TELESCOPE, self._KEYS.SYSTEM]

            for bmf in bandmetaf_keys:
                if bmf not in self:
                    temp = utils.bandmetaf(sband, bmf)
                    if temp is not None:
                        self[bmf] = temp

        # Convert dates to MJD
        timestrs = [str(x) for x in utils.listify(self.get(self._KEYS.TIME, ''))]
        for ti, timestr in enumerate(timestrs):
            if (any(x in timestr for x in ['-', '/']) and not timestr.startswith('-')):
                timestrs[ti] = timestr.replace('/', '-')
                try:
                    # timestrs[ti] = str(utils.astrotime(timestrs[ti], format='isot').mjd)
                    timestrs[ti] = str(utils.astrotime(timestrs[ti], input='isot', output='mjd'))
                except Exception:
                    raise astrocats.catalog.catdict.CatDictError('Unable to convert date to MJD.')
            elif timestr:  # Make sure time is string
                timestrs[ti] = timestr

        if len(timestrs) > 0 and timestrs[0] != '':
            self[self._KEYS.TIME] = timestrs if len(timestrs) > 1 else timestrs[0]

        # Time unit is necessary for maximum time determination
        if self._KEYS.U_TIME not in self and self._KEYS.TIME in self:
            msg = '`{}` not found in photometry, assuming MJD.'.format(self._KEYS.U_TIME)
            self._parent._log.info(msg)
            self[self._KEYS.U_TIME] = 'MJD'

        if (self._KEYS.U_COUNT_RATE not in self and
                self._KEYS.COUNT_RATE in self):
            msg = '`{}` not found in photometry, assuming s^-1.'.format(self._KEYS.U_COUNT_RATE)
            self._parent._log.info(msg)
            self[self._KEYS.U_COUNT_RATE] = 's^-1'

        return

    def _clean_value_for_key(self, key, value):
        value = super(_Photometry, self)._clean_value_for_key(key, value)

        # Do some basic homogenization
        if key == self._KEYS.BAND:
            return utils.bandrepf(value)
        elif key == self._KEYS.INSTRUMENT:
            return utils.instrumentrepf(value)

        return value

    def sort_func(self, key):
        """Specify order for attributes."""
        if key == self._KEYS.TIME:
            return 'aaa'
        if key == self._KEYS.MODEL:
            return 'zzy'
        if key == self._KEYS.SOURCE:
            return 'zzz'
        return key
