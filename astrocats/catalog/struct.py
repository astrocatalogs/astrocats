"""
"""
import os
import sys
from collections import OrderedDict

import astrocats
from astrocats.catalog import utils

_PAS_PATH = "/Users/lzkelley/Research/catalogs/astroschema"
if _PAS_PATH not in sys.path:
    sys.path.append(_PAS_PATH)

import pyastroschema as pas  # noqa

PATH_SCHEMA_INPUT = os.path.join(astrocats._PATH_SCHEMA, "input", "")
PATH_SCHEMA_OUTPUT = os.path.join(astrocats._PATH_SCHEMA, "output", "")

_QUANTITY_RENAMES = [
    ["error_value", "e_value"],
    ["error_lower", "e_lower_value"],
    ["error_upper", "e_upper_value"],
    ["units_value", "u_value"],
    ["units_error", "u_e_value"],
]

path_my_source_schema = os.path.join(PATH_SCHEMA_INPUT, "astrocats_source.json")
path_my_quantity_schema = os.path.join(PATH_SCHEMA_INPUT, "astrocats_quantity.json")
path_my_photometry_schema = os.path.join(PATH_SCHEMA_INPUT, "astrocats_photometry.json")
path_my_spectrum_schema = os.path.join(PATH_SCHEMA_INPUT, "astrocats_spectrum.json")

path_error_schema = os.path.join(PATH_SCHEMA_INPUT, "error.json")
path_realization_schema = os.path.join(PATH_SCHEMA_INPUT, "realization.json")
path_model_schema = os.path.join(PATH_SCHEMA_INPUT, "model.json")
path_correlation_schema = os.path.join(PATH_SCHEMA_INPUT, "correlation.json")


class CatDictError(Exception):
    """Special Error class for non-fatal errors raised in CatDict."""

    def __init__(self, *args, **kwargs):
        """Initialize `CatDictError`."""
        # If `warn` is True, then a warning should be issues.  Otherwise ignore
        # completely
        self.warn = True
        if 'warn' in kwargs:
            self.warn = kwargs.pop('warn')
        Exception.__init__(self, *args, **kwargs)


class Meta_Struct(pas.struct.Struct):

    def __init__(self, parent, key=None, **kwargs):
        super(Meta_Struct, self).__init__(**kwargs)
        if not hasattr(self, "_KEYS"):
            err = ("`Meta_Struct` subclasses must have `_KEYS` set "
                   "before initialization!  `self` = '{}', `parent` = '{}'").format(
                       self, parent)
            raise RuntimeError(err)

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

    def validate(self):
        """Override the parent class `validate` method to raise a `CatDictError`.
        """
        try:
            super(Meta_Struct, self).validate()
        except pas.ValidationError as v_err:
            err = "Caught `ValidationError` during validation!  '{}'".format(str(v_err))
            raise CatDictError(err, warn=True)

        return


# Source
# =================================


@pas.struct.set_struct_schema("source", extensions=[path_my_source_schema])
class Source(Meta_Struct):

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


SOURCE = Source._KEYCHAIN
Source._KEYS = SOURCE


# Quantity
# =================================

class _Quantity_Schema_Rename(pas.schema.SchemaDict):

    def __init__(self, schema={}):
        # print("_Quantity_Schema_Rename.__init__(): '{}'".format(schema))
        # Make sure loaded object is a dict
        schema = pas.utils.get_schema_odict(schema)
        # Convert to str
        schema = pas.utils.json_dump_str(schema)
        # replace strings
        for old, new in _QUANTITY_RENAMES:
            schema = schema.replace(old, new)
        # Convert back to dict
        schema = pas.utils.json_load_str(schema)
        super(_Quantity_Schema_Rename, self).__init__(schema=schema)
        return


@pas.struct.set_struct_schema(
    "quantity", extensions=[path_my_quantity_schema], schema_class=_Quantity_Schema_Rename)
class Quantity(Meta_Struct):

    def __init__(self, parent, key=None, **kwargs):
        super(Quantity, self).__init__(parent, key=key, **kwargs)

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
                    raise CatDictError(err)

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

        self.validate()
        return

    def sort_func(self, key):
        """Sorting logic for `Quantity` objects."""
        if key == self._KEYS.VALUE:
            return 'aaa'
        if key == self._KEYS.SOURCE:
            return 'zzz'
        return key


QUANTITY = Quantity._KEYCHAIN
Quantity._KEYS = QUANTITY


# Photometry
# =================================


@pas.struct.set_struct_schema("photometry", extensions=[path_my_photometry_schema])
class Photometry(Meta_Struct):

    def __init__(self, parent, key=None, **kwargs):
        super(Photometry, self).__init__(parent, key=key, **kwargs)

        '''
        # If `BAND` is given, but any of `bandmetaf_keys` is not, try to infer
        if self._KEYS.BAND in self:
            sband = self[self._KEYS.BAND]
            bandmetaf_keys = [self._KEYS.INSTRUMENT, self._KEYS.TELESCOPE, self._KEYS.SYSTEM]

            for bmf in bandmetaf_keys:
                if bmf not in self:
                    temp = utils.bandmetaf(sband, bmf)
                    if temp is not None:
                        self[bmf] = temp
        '''

        # Convert dates to MJD
        timestrs = [str(x) for x in utils.listify(self.get(self._KEYS.TIME, ''))]
        for ti, timestr in enumerate(timestrs):
            if (any(x in timestr for x in ['-', '/']) and not timestr.startswith('-')):
                timestrs[ti] = timestr.replace('/', '-')
                try:
                    # timestrs[ti] = str(utils.astrotime(timestrs[ti], format='isot').mjd)
                    timestrs[ti] = str(utils.astrotime(timestrs[ti], input='isot', output='mjd'))
                except Exception:
                    raise CatDictError('Unable to convert date to MJD.')
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

    '''
    def _clean_value_for_key(self, key, value):
        value = super(_Photometry, self)._clean_value_for_key(key, value)

        # Do some basic homogenization
        if key == self._KEYS.BAND:
            return utils.bandrepf(value)
        elif key == self._KEYS.INSTRUMENT:
            return utils.instrumentrepf(value)

        return value
    '''

    def sort_func(self, key):
        """Specify order for attributes."""
        if key == self._KEYS.TIME:
            return 'aaa'
        if key == self._KEYS.MODEL:
            return 'zzy'
        if key == self._KEYS.SOURCE:
            return 'zzz'
        return key


PHOTOMETRY = Photometry._KEYCHAIN
Photometry._KEYS = PHOTOMETRY

PHOTOMETRY.FLUX_DENSITY = PHOTOMETRY.FLUXDENSITY
PHOTOMETRY.U_FLUX_DENSITY = PHOTOMETRY.U_FLUXDENSITY
PHOTOMETRY.E_FLUX_DENSITY = PHOTOMETRY.E_FLUXDENSITY
PHOTOMETRY.COUNT_RATE = PHOTOMETRY.COUNTRATE
PHOTOMETRY.E_COUNT_RATE = PHOTOMETRY.E_COUNTRATE
PHOTOMETRY.U_COUNT_RATE = PHOTOMETRY.U_COUNTRATE
PHOTOMETRY.UPPER_LIMIT = PHOTOMETRY.UPPERLIMIT
PHOTOMETRY.LOWER_LIMIT = PHOTOMETRY.LOWERLIMIT
PHOTOMETRY.PHOTON_INDEX = PHOTOMETRY.PHOTONINDEX
PHOTOMETRY.UNABSORBED_FLUX = PHOTOMETRY.UNABSORBEDFLUX
PHOTOMETRY.E_UNABSORBED_FLUX = PHOTOMETRY.E_UNABSORBEDFLUX
PHOTOMETRY.E_UPPER_UNABSORBED_FLUX = PHOTOMETRY.E_UPPER_UNABSORBEDFLUX
PHOTOMETRY.E_LOWER_UNABSORBED_FLUX = PHOTOMETRY.E_LOWER_UNABSORBEDFLUX
PHOTOMETRY.BAND_SET = PHOTOMETRY.BANDSET


# Spectrum
# ==================================


@pas.struct.set_struct_schema("spectrum", extensions=[path_my_spectrum_schema])
class Spectrum(Meta_Struct):
    """Class for storing a single spectrum."""

    def __init__(self, parent, key=None, **kwargs):
        super(Spectrum, self).__init__(parent, key=key, **kwargs)

        # If `data` is not given, construct it from wavelengths, fluxes
        # [errors] `errors` is optional, but if given, then `errorunit` is also
        # req'd
        if self._KEYS.DATA not in self:
            try:
                wavelengths = self[self._KEYS.WAVELENGTHS]
                fluxes = self[self._KEYS.FLUXES]
            except KeyError:
                if self._KEYS.FILENAME in self:
                    return

            errors = self.get(self._KEYS.ERRORS, None)
            if (errors is not None and max([float(err) for err in errors]) > 0.0):
                data = [utils.trim_str_arr(wavelengths), utils.trim_str_arr(fluxes),
                        utils.trim_str_arr(errors)]
            else:
                data = [utils.trim_str_arr(wavelengths), utils.trim_str_arr(fluxes)]

            self[self._KEYS.DATA] = [list(i) for i in zip(*data)]
            if self._KEYS.WAVELENGTHS in self:
                del self[self._KEYS.WAVELENGTHS]
            if self._KEYS.FLUXES in self:
                del self[self._KEYS.FLUXES]
            if self._KEYS.ERRORS in self:
                del self[self._KEYS.ERRORS]

        if self._KEYS.U_TIME not in self:
            msg = '`{}` not found in spectrum, assuming MJD.'.format(self._KEYS.U_TIME)
            self._parent._log.info(msg)
            self[self._KEYS.U_TIME] = 'MJD'

        return

    def is_duplicate_of(self, other):
        """Check if spectrum is duplicate of another."""
        if super(Spectrum, self).is_duplicate_of(other):
            return True

        row_matches = 0
        for ri, row in enumerate(self.get(self._KEYS.DATA, [])):
            lambda1, flux1 = tuple(row[0:2])
            if (self._KEYS.DATA not in other or
                    ri > len(other[self._KEYS.DATA])):
                break
            lambda2, flux2 = tuple(other[self._KEYS.DATA][ri][0:2])
            minlambdalen = min(len(lambda1), len(lambda2))
            minfluxlen = min(len(flux1), len(flux2))
            if (lambda1[:minlambdalen + 1] == lambda2[:minlambdalen + 1] and
                    flux1[:minfluxlen + 1] == flux2[:minfluxlen + 1] and
                    float(flux1[:minfluxlen + 1]) != 0.0):
                row_matches += 1
            # Five row matches should be enough to be sure spectrum is a dupe.
            if row_matches >= 5:
                return True
            # Matches need to happen in the first 10 rows.
            if ri >= 10:
                break
        return False

    def sort_func(self, key):
        """Logic for sorting keys in a `Spectrum` relative to one another."""
        if key == self._KEYS.TIME:
            return 'aaa'
        if key == self._KEYS.DATA:
            return 'zzy'
        if key == self._KEYS.SOURCE:
            return 'zzz'
        return key


SPECTRUM = Spectrum._KEYCHAIN
Spectrum._KEYS = SPECTRUM


# Error
# ==================================


@pas.struct.set_struct_schema(path_error_schema)
class Error(Meta_Struct):
    pass


ERROR = Error._KEYCHAIN
Error._KEYS = ERROR


# Realization
# ==================================


@pas.struct.set_struct_schema(path_realization_schema)
class Realization(Meta_Struct):
    pass


REALIZATION = Realization._KEYCHAIN
Realization._KEYS = REALIZATION


# Model
# ==================================


@pas.struct.set_struct_schema(path_model_schema)
class Model(Meta_Struct):
    """Container for a model with associated metadata.

    `Source` citation required.
    """

    def _init_cat_dict(self, cat_dict_class, key_in_self, **kwargs):
        """Initialize a CatDict object, checking for errors.
        """
        # Catch errors associated with crappy, but not unexpected data
        try:
            new_entry = cat_dict_class(self, key=key_in_self, **kwargs)
        except CatDictError as err:
            if err.warn:
                self._parent._log.info("'{}' Not adding '{}': '{}'".format(
                    self[self._KEYS.NAME], key_in_self, str(err)))
            return None
        return new_entry

    def _add_cat_dict(self, cat_dict_class, key_in_self, check_for_dupes=True, **kwargs):
        """Add a CatDict to this Entry if initialization succeeds and it
        doesn't already exist within the Entry.
        """
        # Try to create a new instance of this subclass of `CatDict`
        new_entry = self._init_cat_dict(cat_dict_class, key_in_self, **kwargs)
        if new_entry is None:
            return False

        # Compare this new entry with all previous entries to make sure is new
        if cat_dict_class != Error:
            for item in self.get(key_in_self, []):
                if new_entry.is_duplicate_of(item):
                    item.append_sources_from(new_entry)
                    # Return the entry in case we want to use any additional
                    # tags to augment the old entry
                    return new_entry

        self.setdefault(key_in_self, []).append(new_entry)

        return True

    def add_realization(self, **kwargs):
        self._add_cat_dict(Realization, self._KEYS.REALIZATIONS, **kwargs)

    def sort_func(self, key):
        if key == self._KEYS.SOURCE:
            return 'zzy'
        if key == self._KEYS.ALIAS:
            return 'zzz'
        return key


MODEL = Model._KEYCHAIN
Model._KEYS = MODEL


# Correlation
# ==================================


@pas.struct.set_struct_schema(path_correlation_schema)
class Correlation(Meta_Struct):
    """Class to store correlation of a `Quantity` with another `Quantity`."""

    def __init__(self, parent, key=None, **kwargs):
        super(Correlation, self).__init__(parent, key=key, **kwargs)

        # Check that value exists
        bad = False
        value = self[self._KEYS.VALUE]
        if len(value) == 0:
            bad = True
        elif (value == '--') or (value == '-'):
            bad = True

        if bad:
            err = "Value '{}' is empty, not adding to '{}'".format(
                value, parent[parent._KEYS.NAME])
            raise CatDictError(err)

        return


CORRELATION = Correlation._KEYCHAIN
Correlation._KEYS = CORRELATION


# Output astrocats schema
# -------------------------------------

def output_schema(path_out, strct_obj, verbose=True):
    # Order of dictionary keys in output file
    sort_order = ['$schema', 'id', 'title', 'description', 'version',
                  'type', 'definitions', 'properties', 'required']

    def sort_item_func(item):
        """Method to generate key for sorting each item of a dictionary
        """
        kk = item[0]

        # Find the index of this key in the `sort_order` list
        idx = None
        for ii, sval in enumerate(sort_order):
            # Ignore text-case
            if kk.lower().startswith(sval.lower()):
                idx = ii
                break

        # Sort things first if found in `sort_order`
        if idx is not None:
            rv = "a_" + str(idx)
        # Then just sort alphabetically by the key
        else:
            rv = "b_" + str(kk)

        return rv

    def sort_dict_func(odict):
        """Method to sort an entire dictionary
        """
        out = OrderedDict(sorted(odict.items(), key=sort_item_func))
        return out

    _schema = strct_obj._SCHEMA
    _title = _schema['title'].lower()
    fname = os.path.join(path_out, _title + ".json")
    _schema.dump(fname, sort_func=sort_dict_func)
    if verbose:
        print("struct.py:output_schema() - saved structure '{}' schema to '{}'".format(
            _title, fname))
    return fname


object_list = [Correlation, Spectrum, Source, Realization, Model, Quantity, Error]
for obj in object_list:
    output_schema(PATH_SCHEMA_OUTPUT, obj)
