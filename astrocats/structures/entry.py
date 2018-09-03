"""Definitions related to the `Entry` class for catalog entries."""
import codecs
import gzip as gz
# import hashlib
import json
import logging
import os
import sys
from collections import OrderedDict
from copy import deepcopy
# from decimal import Decimal
import math

import six

import pyastroschema as pas

from astrocats import utils
from astrocats.structures import struct
from astrocats.structures.struct import (Error, Model, Photometry, Quantity, Source, Spectrum)
from astrocats.structures.struct import (ERROR, MODEL, PHOTOMETRY, QUANTITY, SOURCE, SPECTRUM)

DEP_WARN_ARG = False


@struct.set_struct_schema("astroschema_entry", extensions="astrocats_entry")  # noqa
class _Entry(struct.Meta_Struct):

    _DEPRECATED_ADD_FUNCS = [
        "add_data", "add_alias", "add_error", "add_photometry", "add_quantity",
        "add_source", "add_model", "add_spectrum", "add_listed", "add_self_source",
        "_handle_addition_failure", "_append_additional_tags", "_init_cat_dict",
        "_check_cat_dict_source", "_add_cat_dict"
    ]

    _DEPRECATED_UNUSED_FUNCS = [
        "get_hash", "get_entry_text"
    ]

    def __init__(self, catalog=None, name=None, stub=False):
        # NOTE: FIX: LZK: cannot validate until a valid `name` is set
        super(_Entry, self).__init__(catalog, key=name, validate=False)
        self.catalog = catalog
        self.filename = None
        self.dupe_of = []
        self._stub = stub
        if catalog is not None:
            self._log = catalog.log
        else:
            from astrocats.structures.catalog import Catalog
            self._log = logging.getLogger()
            self._log.error("WARNING: Entry created without a catalog... creating catalog!")
            self.catalog = Catalog(None, self._log)

        self[self._KEYS.NAME] = name

        # NOTE: FIX: LZK: cannot validate until a valid `name` is set
        # self.validate()
        return

    def __repr__(self):
        """Return JSON representation of self."""
        jsonstring = utils.dict_to_pretty_string({self._KEYS.NAME: self})
        return jsonstring

    def __deepcopy__(self, memo):
        """Define how an `Entry` should be deep copied."""
        new_entry = self.__class__(self.catalog)
        for key in self:
            if not key.startswith('__') and key != 'catalog':
                new_entry[key] = deepcopy(self[key], memo)
        return new_entry

    def __getattribute__(self, name):

        try:
            return super(_Entry, self).__getattribute__(name)
        except AttributeError:
            if name in self._DEPRECATED_ADD_FUNCS:
                log = self._log
                msg = "Subclass using the `Entry_Old_Adder` class to preserve functionality."
                log.raise_error("`Entry.{}()` is deprecated! {}".format(name, msg),
                                struct.DeprecationError)
            else:
                raise

    def _get_save_path(self, bury=False):
        """Return the path that this Entry should be saved to."""
        filename = utils.get_filename(self[self._KEYS.NAME])

        # Put objects that shouldn't belong in this catalog in the boneyard
        if bury:
            outdir = self.catalog.get_repo_boneyard()

        # Get normal repository save directory
        else:
            repo_folders = self.catalog.PATHS.get_repo_output_folders()
            # If no repo folders exist, raise an error -- cannot save
            if not len(repo_folders):
                err_str = (
                    "No output data repositories found. Cannot save.\n"
                    "Make sure that repo names are correctly configured "
                    "in the `input/repos.json` file, and either manually or "
                    "automatically (using `astrocats CATALOG git-clone`) "
                    "clone the appropriate data repositories.")
                self.catalog.log.error(err_str)
                raise RuntimeError(err_str)

            outdir = repo_folders[0]

        return outdir, filename

    def _ordered(self, odict):
        """Convert the object into a plain OrderedDict."""
        ndict = OrderedDict()

        if hasattr(odict, "sort_func"):
            key = odict.sort_func
        else:
            key = None

        nkeys = list(sorted(odict.keys(), key=key))
        for key in nkeys:
            if isinstance(odict[key], OrderedDict):
                odict[key] = self._ordered(odict[key])
            if isinstance(odict[key], list):
                if (not (odict[key] and
                         not isinstance(odict[key][0], OrderedDict))):
                    nlist = []
                    for item in odict[key]:
                        if isinstance(item, OrderedDict):
                            nlist.append(self._ordered(item))
                        else:
                            nlist.append(item)
                    odict[key] = nlist
            ndict[key] = odict[key]

        return ndict

    # <<< ==================  STRUCT CREATION/ADDITION  ==================

    # < ----- From File -----

    @classmethod
    def init_from_file(cls, catalog, name=None, path=None, try_gzip=False, **kwargs):
        """Construct a new `Entry` instance from an input file.

        The input file can be given explicitly by `path`, or a path will
        be constructed appropriately if possible.

        Arguments
        ---------
        catalog : `astrocats.structures.catalog.Catalog` instance
            The parent catalog object of which this entry belongs.
        name : str or 'None'
            The name of this entry, e.g. `SN1987A` for a `Supernova` entry.
            If no `path` is given, a path is constructed by trying to find
            a file in one of the 'output' repositories with this `name`.
            note: either `name` or `path` must be provided.
        path : str or 'None'
            The absolutely path of the input file.
            note: either `name` or `path` must be provided.

        """
        if not catalog:
            from astrocats.structures.catalog import Catalog
            log = logging.getLogger()
            catalog = Catalog(None, log)

        log = catalog.log
        # log.debug("init_from_file()")
        if name is None and path is None:
            err = ("Either entry `name` or `path` must be specified to load entry.")
            log.raise_error(err)

        # If the path is given, use that to load from
        load_path = None
        if path is not None:
            load_path = path
            name = ''
        # If the name is given, try to find a path for it
        else:
            repo_paths = catalog.PATHS.get_repo_output_folders()
            for rep in repo_paths:
                filename = utils.get_filename(name)
                newpath = os.path.join(rep, filename + '.json')
                if os.path.isfile(newpath):
                    load_path = newpath
                    break

        if load_path is None or not os.path.isfile(load_path):
            # FIX: is this warning worthy?
            return None

        # Create a new `Entry` instance
        new_entry = cls(catalog, name)

        # Check if .gz file
        if try_gzip and not load_path.endswith('.gz'):
            try_gzip = False

        # Fill it with data from json file
        new_entry._load_data_from_json(load_path, gzip=try_gzip, **kwargs)

        return new_entry

    def _load_data_from_json(self, fname,
                             clean=False, merge=True, pop_schema=True, ignore_keys=[],
                             compare_to_existing=True, gzip=False, filter_on={}, validate=True):
        log = self._log
        log.debug("_load_data_from_json(): {}\n\t{}".format(self.name(), fname))
        # Store the filename this was loaded from
        self.filename = fname

        if gzip:
            jfil = gz.open(fname, 'rb')
        else:
            jfil = codecs.open(fname, 'r')

        data = json.load(jfil, object_pairs_hook=OrderedDict)
        name = list(data.keys())
        if len(name) != 1:
            err = "json file '{}' has multiple top-level keys: {}".format(
                fname, name)
            log.raise_error(err)

        name = name[0]
        # Remove the outmost dict level
        data = data[name]

        # Delete ignored keys
        for key in ignore_keys:
            if key in data:
                del data[key]

        # if name == 'SN2005fz':
        #     DEBUG = True
        # else:
        #     DEBUG = False

        # Convert the OrderedDict data from json into class structure i.e.
        # `Sources` will be extracted and created from the dict Everything
        # that remains afterwards should be okay to just store to this
        # `Entry`
        self._convert_odict_to_classes(
            data, clean=clean, merge=merge, pop_schema=pop_schema,
            compare_to_existing=compare_to_existing, filter_on=filter_on)

        if len(data) > 0:
            err_str = "Remaining entries in `data` after `_convert_odict_to_classes`."
            err_str += "\n{}".format(utils.dict_to_pretty_string(data))
            log.raise_error(err_str)

        jfil.close()

        # If object doesnt have a name yet, but json does, store it
        self_name = self.get(self._KEYS.NAME, None)
        if (self_name is None) or (len(self_name) == 0):
            self[self._KEYS.NAME] = name
        # Warn if there is a name mismatch
        elif self_name.lower().strip() != name.lower().strip():
            log.warning("Object name '{}' does not match name in json: '{}'".format(
                self_name, name))

        if validate:
            self.validate()
        # if DEBUG:
        #     print("\nAFTER VALIDATE\n")
        #     print(utils.dict_to_pretty_string(self[self._KEYS.ERRORS]))
        #    raise

        return

    def _convert_odict_to_classes(self, data, clean=False, merge=True, pop_schema=True,
                                  compare_to_existing=True, filter_on={}, DEBUG=False):
        """Convert `OrderedDict` into `Entry` or its derivative classes."""
        log = self._log
        # log.debug("_convert_odict_to_classes(): {}".format(self.name()))

        # if DEBUG:
        #     print("\n\nDEBUG!!! '{}' '{}' '{}'".format(
        #         self.get('name', None), data['name'], data.get('alias', None)))
        #
        # if DEBUG:
        #     rv = ('errors' in data)
        #     print("errors: ", rv)
        #     if rv:
        #         print(utils.dict_to_pretty_string(data['errors']))

        # Setup filters. Currently only used for photometry.
        fkeys = list(filter_on.keys())

        # Handle 'name'
        name_key = self._KEYS.NAME
        if name_key in data:
            self[name_key] = data.pop(name_key)

        # Cleanup 'internal' repository stuff
        if clean:
            # Add data to `self` in ways accomodating 'internal' formats and
            # leeway.  Removes each added entry from `data` so the remaining
            # stuff can be handled normally
            data = self.clean_internal(data)

        # Handle 'sources'
        # ----------------
        src_key = self._KEYS.SOURCES
        if src_key in data:
            # Remove from `data`
            sources = data.pop(src_key)
            for src in sources:
                self.add_source(allow_alias=True, **src)

        # Handle `photometry`
        # -------------------
        photo_key = self._KEYS.PHOTOMETRY
        if photo_key in data:
            photoms = data.pop(photo_key)
            for photo in photoms:
                skip = False
                for fkey in fkeys:
                    if fkey in photo and photo[fkey] not in filter_on[fkey]:
                        skip = True
                        break
                if skip:
                    continue

                self.create_and_add_struct(
                    Photometry, self._KEYS.PHOTOMETRY,
                    duplicate_check=compare_to_existing, **photo)

        # Handle `spectra`
        # ---------------
        spec_key = self._KEYS.SPECTRA
        if spec_key in data:
            spectra = data.pop(spec_key)
            for spec in spectra:
                self.create_and_add_struct(Spectrum, self._KEYS.SPECTRA, **spec)

        # Handle `error`
        # --------------
        err_key = self._KEYS.ERRORS
        if err_key in data:
            errors = data.pop(err_key)
            for err in errors:
                error, new_error = self.create_and_add_struct(Error, self._KEYS.ERRORS, **err)

        # Handle `models`
        # ---------------
        model_key = self._KEYS.MODELS
        if model_key in data:
            model = data.pop(model_key)
            for mod in model:
                self.create_and_add_struct(
                    Model, self._KEYS.MODELS, duplicate_check=compare_to_existing, **err)

        # Handle everything else --- should be `Quantity`s
        # ------------------------------------------------
        if len(data):
            log.debug("{} remaining entries, assuming `Quantity`".format(len(data)))
            # Iterate over remaining keys
            # `keys` must be converted to a list because `data` may be modified during iteration
            for key in list(data.keys()):
                vals = data.pop(key)
                # All quantities should be in lists of that quantity
                #    E.g. `aliases` is a list of alias quantities
                if not isinstance(vals, list):
                    vals = [vals]
                for vv in vals:
                    self.create_and_add_struct(
                        Quantity, key, duplicate_check=compare_to_existing, **vv)

        if merge and self.dupe_of:
            self.merge_dupes()

        return

    # > ----- from file -----

    def create_struct(self, struct_class, key_in_self, **kwargs):
        log = self._log

        # Perform checks and initialization operations before constructing struct class
        # -----------------------------------------------------------------------------
        valid, kwargs = self._create_struct_bef(struct_class, key_in_self, **kwargs)
        if not valid:
            log.info("`_create_struct_bef()` returned False, skipping")
            return None

        try:
            # `Source` class does not have its own sources
            if isinstance(struct_class, struct.Source) or issubclass(struct_class, struct.Source):
                source = None
            else:
                try:
                    src_key = struct_class._KEYS.SOURCES
                except AttributeError:
                    src_key = struct_class._KEYS.SOURCE

                source = kwargs.get(src_key, None)

        # NOTE: 2018-08-30: This is temporary, remove when working
        except AttributeError:
            print("struct_class = ", struct_class)
            print("_keys = ", struct_class._KEYS)
            print(repr(struct_class))
            print(repr(struct_class._KEYS))
            raise

        if source is not None:
            # If this source/data is erroneous, skip it
            if self.is_erroneous(key_in_self, source):
                log.info("This source is erroneous, skipping")
                return None

            # If this source/data is private, skip it
            if ((not self.catalog.args.private) and self.is_private(key_in_self, source)):
                log.info("This source is private, skipping")
                return None

        # Constructing struct class
        # -----------------------------------------------------------------------------
        try:
            new_struct = struct_class(self, key=key_in_self, **kwargs)
        except struct.CleaningError as err:
            log.info("Caught cleaning error: '{}'".format(err))
            msg = "Current task = '{}', ".format(self._parent.current_task)
            msg += "name = '{}', ".format(self[self._KEYS.NAME])
            msg += "`struct_class` = '{}', ".format(struct_class)
            msg += "`key_in_self` = '{}'".format(key_in_self)
            log.info(msg)
            new_struct = None
        except struct.CatDictError as err:
            log_lvl = log.WARNING if err.warn else log.INFO
            msg = "'{}' Not adding '{}': '{}'".format(self[self._KEYS.NAME], key_in_self, str(err))
            log.log(log_lvl, msg)

            err_str = "Entry.create_struct() failed!  "
            err_str += "Current task = '{}', ".format(self._parent.current_task)
            err_str += "name = '{}', ".format(self[self._KEYS.NAME])
            err_str += "`struct_class` = '{}', ".format(struct_class)
            err_str += "`key_in_self` = '{}'".format(key_in_self)
            if self.catalog.RAISE_ERROR_ON_ADDITION_FAILURE:
                log.raise_error(err_str)
            else:
                log.info(err_str)
                new_struct = None
        except Exception as err:
            log = self._log
            log.error("ERROR in `Entry.create_struct()`!")
            log.error("Current task = '{}'".format(self._parent.current_task))
            log.error("name = '{}'".format(self[self._KEYS.NAME]))
            log.error("`struct_class` = '{}'".format(struct_class))
            log.error("`key_in_self` = '{}'".format(key_in_self))
            # log.error("`kwargs = '{}'".format(kwargs))
            # log.error("self = \n'{}'".format(repr(self)))
            log.error("\n")
            log.error(str(err))
            raise

        return new_struct

    def add_struct(self, key_in_self, new_struct,
                   duplicate_check=True, duplicate_merge=True, **kwargs):
        """
        Returns
        -------
        stored_struct :
        given_struct :

        """

        log = self._log

        # Perform checks and finalization operations before adding struct to self
        # -----------------------------------------------------------------------------
        valid, new_struct = self._add_struct_bef(new_struct, **kwargs)
        if not valid:
            log.info("`_add_struct_bef()` returned False, not adding")
            # return None
            return None, new_struct

        # Compare this new entry with all previous entries to make sure is new
        if duplicate_check:
            for ii, item in enumerate(self.get(key_in_self, [])):
                if new_struct.is_duplicate_of(item):
                    log.info("Duplicate found, key: {}, index: {}".format(key_in_self, ii))
                    if duplicate_merge:
                        self._merge_structs(key_in_self, ii, new_struct)

                    # Return the entry in case we want to use any additional info from it
                    return item, new_struct

        self.setdefault(key_in_self, []).append(new_struct)
        return new_struct, None

    def create_and_add_struct(self, struct_class, key_in_self,
                              duplicate_check=True, duplicate_merge=True, **kwargs):
        new_struct = self.create_struct(struct_class, key_in_self, **kwargs)
        if new_struct is None:
            stored = None
        else:
            stored, new_struct = self.add_struct(
                key_in_self, new_struct,
                duplicate_check=duplicate_check, duplicate_merge=duplicate_merge)

        return stored, new_struct

    def _create_struct_bef(self, *args, **kwargs):
        return True, kwargs

    def _add_struct_bef(self, new_struct, **kwargs):
        return True, new_struct

    def _merge_structs(self, key_in_self, index, new_struct):
        item = self[key_in_self][index]
        item.append_aliases_from(new_struct)
        if isinstance(new_struct, Quantity):
            # self._append_additional_tags(name, source, new_struct)
            self._merge_quantities(item, new_struct)

        # self._log.debug("Merged from: '{}'".format(new_struct))
        # self._log.debug("Merged into: '{}'".format(item))

        return

    def _merge_quantities(self, dst, src):
        # NOTE: FIX: LZK: temporary method here!
        # See note in `supernovae.supernova.Supernova._merge_quantities`.
        return

    def _get_alias_from_add_struct_return(self, alias, struct_return):
        log = self._log

        # Found duplicate `Source`, the original is returned, return that alias
        if isinstance(struct_return, struct.Meta_Struct):
            alias = struct_return[SOURCE.ALIAS]
        # Source was added successfully, and `alias` (from above) is still accurate
        elif struct_return is None:
            log.raise_error("`add_struct` failed!")
        else:
            err = "Unexpected return type '{}' ({}) from `add_struct`".format(
                type(struct_return), struct_return)
            log.raise_error(err)

        return alias

    # >>> ==================  STRUCT CREATION/ADDITION  ==================

    def merge_dupes(self):
        """Merge two entries that correspond to the same entry."""
        for dupe in self.dupe_of:
            if dupe in self.catalog.entries:
                if self.catalog.entries[dupe]._stub:
                    # merge = False to avoid infinite recursion
                    self.catalog.load_entry_from_name(dupe, delete=True, merge=False)
                self.catalog.copy_entry_to_entry(self.catalog.entries[dupe], self)
                del self.catalog.entries[dupe]
        self.dupe_of = []
        return

    def clean_internal(self, data=None):
        """Clean input from 'internal', human added data.

        This is used in the 'Entry.init_from_file' method.
        """
        return data

    def extra_aliases(self):
        """Return aliases considered when merging duplicates."""
        return []

    def get_aliases(self, includename=True):
        """Retrieve the aliases of this object as a list of strings.

        Arguments
        ---------
        includename : bool
            Include the 'name' parameter in the list of aliases.
        """
        # empty list if doesnt exist
        alias_quanta = self.get(self._KEYS.ALIAS, [])
        aliases = [aq[QUANTITY.VALUE] for aq in alias_quanta]
        if includename and self[self._KEYS.NAME] not in aliases:
            aliases = [self[self._KEYS.NAME]] + aliases
        return aliases

    def get_source_by_alias(self, alias):
        """Given an alias, find the corresponding source in this entry.

        If the given alias doesn't exist (e.g. there are no sources), then a
        `ValueError` is raised.

        Arguments
        ---------
        alias : str
            The str-integer (e.g. '8') of the target source.

        Returns
        -------
        source : `astrocats.structures.source.Source` object
            The source object corresponding to the passed alias.

        """
        for source in self.get(self._KEYS.SOURCES, []):
            if source[self._KEYS.ALIAS] == alias:
                return source
        raise ValueError("Source '{}': alias '{}' not found!".format(self[self._KEYS.NAME], alias))

    def get_stub(self):
        """Get a new `Entry` which contains the 'stub' of this one.

        The 'stub' is only the name and aliases.

        Usage:
        -----
        To convert a normal entry into a stub (for example), overwrite the
        entry in place, i.e.
        >>> entries[name] = entries[name].get_stub()

        Returns
        -------
        stub : `astrocats.structures.entry.Entry` subclass object
            The type of the returned object is this instance's type.

        """
        stub = type(self)(self.catalog, self[self._KEYS.NAME], stub=True)
        if self._KEYS.ALIAS in self:
            stub[self._KEYS.ALIAS] = self[self._KEYS.ALIAS]
        if self._KEYS.DISTINCT_FROM in self:
            stub[self._KEYS.DISTINCT_FROM] = self[self._KEYS.DISTINCT_FROM]
        if self._KEYS.RA in self:
            stub[self._KEYS.RA] = self[self._KEYS.RA]
        if self._KEYS.DEC in self:
            stub[self._KEYS.DEC] = self[self._KEYS.DEC]
        if self._KEYS.DISCOVER_DATE in self:
            stub[self._KEYS.DISCOVER_DATE] = self[self._KEYS.DISCOVER_DATE]
        if self._KEYS.SOURCES in self:
            stub[self._KEYS.SOURCES] = self[self._KEYS.SOURCES]
        return stub

    def is_erroneous(self, field, sources):
        """Check if attribute has been marked as being erroneous."""
        if self._KEYS.ERRORS in self:
            my_errors = self[self._KEYS.ERRORS]
            for alias in sources.split(','):
                source = self.get_source_by_alias(alias)
                bib_err_values = [
                    err[ERROR.VALUE] for err in my_errors
                    if (err[ERROR.KIND] == SOURCE.BIBCODE) and (err[ERROR.EXTRA] == field)
                ]
                if (SOURCE.BIBCODE in source and source[SOURCE.BIBCODE] in bib_err_values):
                    return True

                name_err_values = [
                    err[ERROR.VALUE] for err in my_errors
                    if (err[ERROR.KIND] == SOURCE.NAME) and (err[ERROR.EXTRA] == field)
                ]
                if (SOURCE.NAME in source and source[SOURCE.NAME] in name_err_values):
                    return True

        return False

    def is_private(self, key, sources):
        """Check if attribute is private."""
        # aliases are always public.
        if key == self._KEYS.ALIAS:
            return False
        return all([SOURCE.PRIVATE in self.get_source_by_alias(x) for x in sources.split(',')])

    def name(self):
        """Return own name."""
        return self.get(self._KEYS.NAME, None)

    def num_sources(self):
        """Return the current number of sources stored in this instance.

        Returns
        -------
        len : int
            The *integer* number of existing sources.
        """
        return len(self.get(self._KEYS.SOURCES, []))

    def priority_prefixes(self):
        """Return prefixes to given priority when merging duplicate entries."""
        return ()

    def sanitize(self):
        """Sanitize the data (sort it, etc.) before writing it to disk.

        Template method that can be overridden in each catalog's subclassed
        `Entry` object.
        """
        name = self[self._KEYS.NAME]

        def sort_func_spec(spec):
            sort = (float(spec.get(SPECTRUM.TIME, 0.0)), spec.get(SPECTRUM.FILENAME, ''))
            return sort

        def sort_func_phot(phot):
            types = (six.string_types, float, int)
            time = phot.get(PHOTOMETRY.TIME, 0.0)
            if isinstance(time, types):
                time = float(time)
            else:
                time = min([float(y) for y in time])

            band = phot.get(PHOTOMETRY.BAND, '')
            mag = float(phot[PHOTOMETRY.MAGNITUDE]) if PHOTOMETRY.MAGNITUDE in phot else -math.inf
            sort = (time, band, mag)
            return sort

        aliases = self.get_aliases(includename=False)
        if name not in aliases:
            # Assign the first source to alias, if not available assign astrocats.
            if self._KEYS.SOURCES in self:
                self.add_quantity(self._KEYS.ALIAS, name, '1')
                if self._KEYS.ALIAS not in self:
                    source = self.add_self_source()
                    self.add_quantity(self._KEYS.ALIAS, name, source)
            else:
                source = self.add_self_source()
                self.add_quantity(self._KEYS.ALIAS, name, source)

        if self._KEYS.ALIAS in self:
            self[self._KEYS.ALIAS].sort(
                key=lambda key: utils.alias_priority(name, key[QUANTITY.VALUE]))
        else:
            self._log.raise_error('There should be at least one alias for `{}`.'.format(name))

        if self._KEYS.PHOTOMETRY in self:
            try:
                self[self._KEYS.PHOTOMETRY].sort(key=sort_func_phot)
            except TypeError:
                print("\n", self[self._KEYS.PHOTOMETRY], "\n")
                print("\n", [xx.get(PHOTOMETRY.TIME, None) for xx in self[self._KEYS.PHOTOMETRY]], "\n")
                raise

        if self._KEYS.SPECTRA in self:
            spec_times = [SPECTRUM.TIME in x for x in self[self._KEYS.SPECTRA]]
            if any(spec_times):
                self[self._KEYS.SPECTRA].sort(key=sort_func_spec)

        if self._KEYS.SOURCES in self:
            # Remove orphan sources
            source_aliases = [x[SOURCE.ALIAS] for x in self[self._KEYS.SOURCES]]
            # Sources with the `PRIVATE` attribute are always retained
            source_list = [
                x[SOURCE.ALIAS] for x in self[self._KEYS.SOURCES]
                if SOURCE.PRIVATE in x
            ]
            _SOURCES_SKIP_KEYS = [self._KEYS.NAME, self._KEYS.SOURCES, self._KEYS.ERRORS]
            for key in self:
                # if self._KEYS.get_key_by_name(key).no_source:
                if key in _SOURCES_SKIP_KEYS:
                    continue
                for item in self[key]:
                    if item._KEYS.SOURCE in item:
                        source_list += item[item._KEYS.SOURCE].split(',')
            new_src_list = sorted(list(set(source_aliases).intersection(source_list)))
            new_sources = []
            for source in self[self._KEYS.SOURCES]:
                if source[SOURCE.ALIAS] in new_src_list:
                    new_sources.append(source)
                else:
                    self._log.info('Removing orphaned source from `{}`.'.format(name))

            if not new_sources:
                del self[self._KEYS.SOURCES]

            self[self._KEYS.SOURCES] = new_sources
            return

    def save(self, bury=False, final=False):
        """Write entry to JSON file in the proper location.

        Arguments
        ---------
        bury : bool

        final : bool
            If this is the 'final' save, perform additional sanitization and
            cleaning operations.

        """
        outdir, filename = self._get_save_path(bury=bury)

        if final:
            self.sanitize()

        # jsonstring = json.dumps(
        jsonstring = pas.utils.json_dump_str(
            {self[self._KEYS.NAME]: self._ordered(self)},
            indent='\t' if sys.version_info[0] >= 3 else 4,
            separators=(',', ':'),
            ensure_ascii=False
        )
        if not os.path.isdir(outdir):
            raise RuntimeError("Output directory '{}' for event '{}' does "
                               "not exist.".format(outdir, self[self._KEYS.NAME]))
        save_name = os.path.join(outdir, filename + '.json')
        with codecs.open(save_name, 'w', encoding='utf8') as sf:
            sf.write(jsonstring)

        if not os.path.exists(save_name):
            raise RuntimeError("File '{}' was not saved!".format(save_name))

        return save_name

    def set_preferred_name(self):
        """Set a preferred name for the entry."""
        return self[self._KEYS.NAME]

    def sort_func(self, key):
        """Used to sort keys when writing Entry to JSON format.

        Should be supplemented/overridden by inheriting classes.
        """
        # if key == self._KEYS.SCHEMA:
        #     return 'aaa'
        if key == self._KEYS.NAME:
            return 'aab'
        if key == self._KEYS.SOURCES:
            return 'aac'
        if key == self._KEYS.ALIAS:
            return 'aad'
        if key == self._KEYS.MODELS:
            return 'aae'
        if key == self._KEYS.PHOTOMETRY:
            return 'zzy'
        if key == self._KEYS.SPECTRA:
            return 'zzz'
        return key

    # ======================  DEPRECATIONS  ========================
    # ==============================================================

    '''
    def get_hash(self, keys=[]):
        """Return a unique hash associated with the listed keys."""
        if not len(keys):
            keys = list(self.keys())

        string_rep = ''
        oself = self._ordered(deepcopy(self))
        for key in keys:
            string_rep += json.dumps(oself.get(key, ''), sort_keys=True)

        return hashlib.sha512(string_rep.encode()).hexdigest()[:16]
    '''

    '''
    def _clean_quantity(self, quantity):
        """Clean quantity value before it is added to entry."""
        value = quantity.get(QUANTITY.VALUE, '').strip()
        error = quantity.get(QUANTITY.E_VALUE, '').strip()
        unit = quantity.get(QUANTITY.U_VALUE, '').strip()
        kind = quantity.get(QUANTITY.KIND, '')

        if isinstance(kind, list) and not isinstance(kind, six.string_types):
            kind = [x.strip() for x in kind]
        else:
            kind = kind.strip()

        if not value:
            return False

        if utils.is_number(value):
            value = '%g' % Decimal(value)
        if error:
            error = '%g' % Decimal(error)

        if value:
            quantity[QUANTITY.VALUE] = value
        if error:
            quantity[QUANTITY.E_VALUE] = error
        if unit:
            quantity[QUANTITY.U_VALUE] = unit
        if kind:
            quantity[QUANTITY.KIND] = kind

        return True
    '''

    def _clean_quantity(self, *args, **kwargs):
        log = self._log
        msg = ""
        log.raise_error("`Entry._clean_quantity()` is deprecated! " + msg, struct.DeprecationError)

    '''
    def get_entry_text(self, fname):
        """Retrieve the raw text from a file."""
        if fname.split('.')[-1] == 'gz':
            with gz.open(fname, 'rt') as f:
                filetext = f.read()
        else:
            with codecs.open(fname, 'r') as f:
                filetext = f.read()
        return filetext
    '''


class _Entry_Old_Adder(_Entry):

    def add_data(self, key_in_self, value=None, source=None, struct_class=Quantity,
                 check_for_dupes=True, dupes_merge_tags=True, **kwargs):
        """Add a `CatDict` data container (dict) to this `Entry`.

        CatDict only added if initialization succeeds and it
        doesn't already exist within the Entry.
        """
        log = self._log
        # log.debug("Entry.add_data()")
        FAIL = 0
        DUPE = -1
        SUCC = 1

        if value is not None:
            if QUANTITY.VALUE in kwargs:
                log.raise_error("`value` given as both arg and kwarg!\n{}".format(kwargs))
            kwargs[QUANTITY.VALUE] = value

        if source is not None:
            if QUANTITY.SOURCE in kwargs:
                log.raise_error("`source` given as both arg and kwarg!\n{}".format(kwargs))
            kwargs[QUANTITY.SOURCE] = source

        # Make sure that a source, if given, is valid
        retval = self._check_cat_dict_source(struct_class, key_in_self, **kwargs)
        if not retval:
            self._handle_addition_failure(
                "Entry._check_cat_dict_source()", struct_class, key_in_self, **kwargs)
            return None, FAIL

        # Try to create a new instance of this subclass of `CatDict`
        new_data = self._init_cat_dict(struct_class, key_in_self, **kwargs)
        if new_data is None:
            self._handle_addition_failure(
                "Entry._init_cat_dict()", struct_class, key_in_self, **kwargs)
            return None, FAIL

        # Compare this new entry with all previous entries to make sure is new
        if check_for_dupes:
            for item in self.get(key_in_self, []):
                # Do not add duplicates
                if new_data.is_duplicate_of(item):
                    item.append_aliases_from(new_data)
                    return new_data, DUPE

        # Add data
        self.setdefault(key_in_self, []).append(new_data)

        return new_data, SUCC

    def add_alias(self, alias, source, clean=True, check_for_dupes=True, **kwargs):
        """Add an alias, optionally 'cleaning' the alias string.

        Calls the parent `catalog` method `clean_entry_name` - to apply the
        same name-cleaning as is applied to entry names themselves.

        Returns
        -------
        alias : str
            The stored version of the alias (cleaned or not).

        """
        if clean:
            alias = self.catalog.clean_entry_name(alias)

        # self.add_quantity(self._KEYS.ALIAS, alias, source)
        kwargs.update({QUANTITY.VALUE: alias, QUANTITY.SOURCE: source})
        new_data, retval = self.add_data(
            self._KEYS.ALIAS, check_for_dupes=check_for_dupes, **kwargs)

        # Check if this adding this alias makes us a dupe, if so mark ourselves as a dupe.
        if check_for_dupes and (alias in self.catalog.aliases):
            poss_dupe = self.catalog.aliases[alias]
            if (poss_dupe != self[self._KEYS.NAME]) and (poss_dupe in self.catalog.entries):
                self.dupe_of.append(poss_dupe)

        # If we're not checking for duplicates, warn about overwriting
        elif alias in self.catalog.aliases:
            self._log.warning("Warning, overwriting alias '{}': '{}'".format(
                alias, self.catalog.aliases[alias]))

        # Add this alias to parent catalog's reverse dictionary linking aliases to names
        self.catalog.aliases[alias] = self[self._KEYS.NAME]

        if self.dupe_of:
            self.merge_dupes()

        return alias

    def add_error(self, value, **kwargs):
        """Add an `Error` instance to this entry."""
        kwargs.update({ERROR.VALUE: value})
        self._add_cat_dict(Error, self._KEYS.ERRORS, **kwargs)
        return

    def add_photometry(self, compare_to_existing=True, **kwargs):
        """Add a `Photometry` instance to this entry."""
        self._add_cat_dict(Photometry, self._KEYS.PHOTOMETRY,
                           compare_to_existing=compare_to_existing, **kwargs)
        return

    def add_quantity(self, quantities, value, source,
                     check_for_dupes=True, compare_to_existing=True, **kwargs):
        """Add an `Quantity` instance to this entry."""
        success = True
        for quantity in utils.listify(quantities):
            kwargs.update({QUANTITY.VALUE: value, QUANTITY.SOURCE: source})
            cat_dict = self._add_cat_dict(
                Quantity, quantity,
                compare_to_existing=compare_to_existing, check_for_dupes=check_for_dupes,
                **kwargs)
            if isinstance(cat_dict, struct.Meta_Struct):
                self._append_additional_tags(quantity, source, cat_dict)
                success = False
            elif cat_dict is False:
                success = False

        return success

    def add_source(self, allow_alias=False, **kwargs):
        """Add a `Source` instance to this entry."""
        log = self._log
        # log.debug("Entry.add_source()")
        if (not allow_alias) and (SOURCE.ALIAS in kwargs):
            err_str = "`{}` passed in kwargs, this shouldn't happen!".format(SOURCE.ALIAS)
            log.error(err_str)
            raise RuntimeError(err_str)

        # Set alias number to be +1 of current number of sources
        if SOURCE.ALIAS not in kwargs:
            kwargs[SOURCE.ALIAS] = str(self.num_sources() + 1)

        # log.warning("Entry.add_source() - Calling `Entry._init_cat_dict()`")
        source_obj = self._init_cat_dict(Source, self._KEYS.SOURCES, **kwargs)
        if source_obj is None:
            return None

        for item in self.get(self._KEYS.SOURCES, []):
            if source_obj.is_duplicate_of(item):
                return item[SOURCE.ALIAS]

        self.setdefault(self._KEYS.SOURCES, []).append(source_obj)
        return source_obj[SOURCE.ALIAS]

    def add_model(self, allow_alias=False, **kwargs):
        """Add a `Model` instance to this entry."""
        if not allow_alias and MODEL.ALIAS in kwargs:
            err_str = "`{}` passed in kwargs, this shouldn't happen!".format(
                SOURCE.ALIAS)
            self._log.error(err_str)
            raise RuntimeError(err_str)

        # Set alias number to be +1 of current number of models
        if MODEL.ALIAS not in kwargs:
            kwargs[MODEL.ALIAS] = str(self.num_models() + 1)
        model_obj = self._init_cat_dict(Model, self._KEYS.MODELS, **kwargs)
        if model_obj is None:
            return None

        for item in self.get(self._KEYS.MODELS, ''):
            if model_obj.is_duplicate_of(item):
                return item[item._KEYS.ALIAS]

        self.setdefault(self._KEYS.MODELS, []).append(model_obj)
        return model_obj[model_obj._KEYS.ALIAS]

    def add_spectrum(self, compare_to_existing=True, **kwargs):
        """Add a `Spectrum` instance to this entry."""
        spec_key = self._KEYS.SPECTRA
        # Make sure that a source is given, and is valid (nor erroneous)
        retval = self._check_cat_dict_source(Spectrum, spec_key, **kwargs)
        if not retval:
            return None

        # Try to create a new instance of `Spectrum`
        new_spectrum = self._init_cat_dict(Spectrum, spec_key, **kwargs)
        if new_spectrum is None:
            return None

        is_dupe = False
        for item in self.get(spec_key, []):
            # Only the `filename` should be compared for duplicates. If a
            # duplicate is found, that means the previous `exclude` array
            # should be saved to the new object, and the old deleted
            if new_spectrum.is_duplicate_of(item):
                if SPECTRUM.EXCLUDE in new_spectrum:
                    item[SPECTRUM.EXCLUDE] = new_spectrum[SPECTRUM.EXCLUDE]
                elif SPECTRUM.EXCLUDE in item:
                    item.update(new_spectrum)
                is_dupe = True
                break

        if not is_dupe:
            self.setdefault(spec_key, []).append(new_spectrum)
        return

    def add_listed(self, key, value, check=True):
        listed = self.setdefault(key, [])
        # Make sure `value` isn't already in the list
        if check and (value in listed):
            return

        listed.append(value)
        return

    def add_self_source(self):
        """Add a source that refers to the catalog itself.

        For now this points to the Open Supernova Catalog by default.
        """
        return self.add_source(
            bibcode=self.catalog.OSC_BIBCODE,
            name=self.catalog.OSC_NAME,
            url=self.catalog.OSC_URL,
            secondary=True)

    def _handle_addition_failure(self, fail_loc, cat_class, cat_key, **kwargs):
        """Based on `catalog.ADDITION_FAILURE_BEHAVIOR`, react to a failure appropriately.

        A `logging.DEBUG` level message is always given.

        `ADDITION_FAILURE_BEHAVIOR` == `ADD_FAIL_ACTION.WARN`
            Then a `logging.WARNING` message is raised.
        `ADDITION_FAILURE_BEHAVIOR` == `ADD_FAIL_ACTION.IGNORE`
            No addition action is taken.
        `ADDITION_FAILURE_BEHAVIOR` == `ADD_FAIL_ACTION.RAISE`
            Then an error is raised.
            This is the default behavior that also acts if an unknown value is given.

        """
        err_str = "'{}' failed!\n".format(fail_loc)
        err_str += "class: '{}', key: '{}'\nkwargs: '{}'".format(cat_class, cat_key, kwargs)

        fail_flag = self.catalog.ADDITION_FAILURE_BEHAVIOR
        # Log a message at 'debug' level regardless
        self._log.debug(err_str)
        self._log.debug("`ADDITION_FAILURE_BEHAVIOR` = '{}'".format(fail_flag))

        # if `WARN` then also log a warning
        if fail_flag == utils.ADD_FAIL_ACTION.WARN:
            self._log.warning(err_str)
        # Raise an error
        elif fail_flag == utils.ADD_FAIL_ACTION.IGNORE:
            pass
        # default behavior is to raise an error
        else:
            self._log.raise_error(err_str, RuntimeError)

        return

    def _append_additional_tags(self, quantity, source, cat_dict):
        """Append additional bits of data to an existing quantity.

        Called when a newly added quantity is found to be a duplicate.
        """
        pass

    def _init_cat_dict(self, struct_class, key_in_self, **kwargs):
        """Initialize a CatDict object, checking for errors."""
        log = self._log
        # log.debug("Entry._init_cat_dict()")

        # Remove empty string values
        bad_keys = [kk for kk, vv in kwargs.items() if isinstance(vv, str) and len(vv) == 0]
        for bk in bad_keys:
            kwargs.pop(bk)

        # Catch errors associated with crappy, but not unexpected data
        try:
            new_struct = struct_class(self, key=key_in_self, **kwargs)
            if new_struct is None:
                err = "`struct_class` = '{}' returned None on init!".format(struct_class)
                raise RuntimeError(err)
        except struct.CatDictError as err:
            log_lvl = log.WARNING if err.warn else log.INFO
            msg = "Task: '{}', function: '{}'".format(
                self.catalog.get_current_task_str(), self.catalog.current_task.function)
            log.log(log_lvl, msg)
            msg = "'{}' Not adding '{}': '{}'".format(self[self._KEYS.NAME], key_in_self, str(err))
            log.log(log_lvl, msg)
            return None
        except Exception as err:
            log = self._log
            log.error("ERROR in `Entry._init_cat_dict()`!")
            log.error("Current task = '{}'".format(self._parent.current_task))
            log.error("name = '{}'".format(self[self._KEYS.NAME]))
            log.error("\n\n\n")
            log.error("`struct_class` = '{}'".format(struct_class))
            log.error("`key_in_self` = '{}'".format(key_in_self))
            log.error("`kwargs = '{}'".format(kwargs))
            log.error("\n\n")
            log.error("self = \n'{}'".format(repr(self)))
            log.error("\n")
            raise

        # NOTE: LZK COMMENTED OUT for testing new `Quantity` astroschema
        """
        try:
            new_struct = struct_class(self, key=key_in_self, **kwargs)
        except struct.CatDictError as err:
            if err.warn:
                self._log.info("'{}' Not adding '{}': '{}'".format(self[
                    self._KEYS.NAME], key_in_self, str(err)))
            return None
        """
        # log.warning("return new_struct")
        return new_struct

    def _check_cat_dict_source(self, struct_class, key_in_self, **kwargs):
        """Check that the quantity isn't erroneous or private."""
        # Make sure that a source is given
        source = kwargs.get(struct_class._KEYS.SOURCE, None)
        if source is not None:
            # If this source/data is erroneous, skip it
            if self.is_erroneous(key_in_self, source):
                self._log.info("This source is erroneous, skipping")
                return False
            # If this source/data is private, skip it
            if ((not self.catalog.args.private) and self.is_private(key_in_self, source)):
                self._log.info("This source is private, skipping")
                return False

        return True

    def _add_cat_dict(self, struct_class, key_in_self,
                      check_for_dupes=True, compare_to_existing=True, **kwargs):
        """Add a `CatDict` to this `Entry`.

        CatDict only added if initialization succeeds and it
        doesn't already exist within the Entry.
        """
        log = self._log
        # Make sure that a source is given, and is valid (nor erroneous)
        failure = False
        if struct_class != Error:
            try:
                source = self._check_cat_dict_source(struct_class, key_in_self, **kwargs)
            except struct.CatDictError as err:
                if err.warn:
                    msg = "'{}' Not adding '{}': '{}'".format(
                        self[self._KEYS.NAME], key_in_self, str(err))
                    log.info(msg)
                # return False
                failure = True

            if source is None:
                log.warning("Source `None` in `Entry._add_cat_dict()`!")
                log.warning("key: '{}', kwargs: '{}'".format(key_in_self, kwargs))
                # return False
                failure = True

        # Try to create a new instance of this subclass of `CatDict`
        new_struct = self._init_cat_dict(struct_class, key_in_self, **kwargs)
        if new_struct is None:
            # if self.catalog.RAISE_ERROR_ON_ADDITION_FAILURE:
            #     err_str = "Entry._init_cat_dict() failed!"
            #     err_str += "class: '{}', key_in_self: '{}', kwargs: '{}'".format(
            #         struct_class, key_in_self, kwargs)
            #     self._log.raise_error(err_str, RuntimeError)
            # return False
            failure = True

        if failure:
            err_str = "Entry._init_cat_dict() failed!"
            err_str += "class: '{}', key_in_self: '{}', kwargs: '{}'".format(
                struct_class, key_in_self, kwargs)
            if self.catalog.RAISE_ERROR_ON_ADDITION_FAILURE:
                log.raise_error(err_str)
            else:
                log.info(err_str)
                return False

        # Compare this new entry with all previous entries to make sure is new
        #    If it is NOT new, return the entry
        if compare_to_existing and struct_class != Error:
            for item in self.get(key_in_self, []):
                if new_struct.is_duplicate_of(item):
                    item.append_aliases_from(new_struct)
                    # Return the entry in case we want to use any additional
                    # tags to augment the old entry
                    return new_struct

        self.setdefault(key_in_self, []).append(new_struct)
        return True

    def _load_data_from_json(self, fhand,
                             clean=False, merge=True, pop_schema=True, ignore_keys=[],
                             compare_to_existing=True, gzip=False, filter_on={}):
        # FIX: check for overwrite??"""
        self._log.debug("_load_data_from_json(): {}\n\t{}".format(self.name(), fhand))
        # Store the filename this was loaded from
        self.filename = fhand

        if gzip:
            jfil = gz.open(fhand, 'rb')
        else:
            jfil = codecs.open(fhand, 'r')

        data = json.load(jfil, object_pairs_hook=OrderedDict)
        name = list(data.keys())
        if len(name) != 1:
            err = "json file '{}' has multiple keys: {}".format(fhand, list(name))
            self._log.error(err)
            raise ValueError(err)
        name = name[0]
        # Remove the outmost dict level
        data = data[name]
        # self._log.debug("Name: {}".format(name))

        # Delete ignored keys
        for key in ignore_keys:
            if key in data:
                del data[key]

        # Convert the OrderedDict data from json into class structure i.e.
        # `Sources` will be extracted and created from the dict Everything
        # that remains afterwards should be okay to just store to this
        # `Entry`
        self._convert_odict_to_classes(
            data, clean=clean, merge=merge, pop_schema=pop_schema,
            compare_to_existing=compare_to_existing, filter_on=filter_on)
        if len(data):
            err_str = "Remaining entries in `data` after `_convert_odict_to_classes`."
            err_str += "\n{}".format(utils.dict_to_pretty_string(data))
            self._log.error(err_str)
            raise RuntimeError(err_str)

        jfil.close()

        # If object doesnt have a name yet, but json does, store it
        self_name = self[self._KEYS.NAME]
        if len(self_name) == 0:
            self[self._KEYS.NAME] = name
        # Warn if there is a name mismatch
        elif self_name.lower().strip() != name.lower().strip():
            self._log.warning("Object name '{}' does not match name in json: '{}'".format(
                self_name, name))

        self.validate()
        return

    def _convert_odict_to_classes(self, data, clean=False, merge=True, pop_schema=True,
                                  compare_to_existing=True, filter_on={}):
        """Convert `OrderedDict` into `Entry` or its derivative classes."""
        # self._log.debug("_convert_odict_to_classes(): {}".format(self.name()))
        # self._log.debug("This should be a temporary fix.  Dont be lazy.")

        # Setup filters. Currently only used for photometry.
        fkeys = list(filter_on.keys())

        # Handle 'name'
        name_key = self._KEYS.NAME
        if name_key in data:
            self[name_key] = data.pop(name_key)

        # Cleanup 'internal' repository stuff
        if clean:
            # Add data to `self` in ways accomodating 'internal' formats and
            # leeway.  Removes each added entry from `data` so the remaining
            # stuff can be handled normally
            data = self.clean_internal(data)

        # Handle 'sources'
        # ----------------
        src_key = self._KEYS.SOURCES
        if src_key in data:
            # Remove from `data`
            sources = data.pop(src_key)
            self._log.debug("Found {} '{}' entries".format(
                len(sources), src_key))
            self._log.debug("{}: {}".format(src_key, sources))

            for src in sources:
                self.add_source(allow_alias=True, **src)

        # Handle `photometry`
        # -------------------
        photo_key = self._KEYS.PHOTOMETRY
        if photo_key in data:
            photoms = data.pop(photo_key)
            self._log.debug("Found {} '{}' entries".format(
                len(photoms), photo_key))
            phcount = 0
            for photo in photoms:
                skip = False
                for fkey in fkeys:
                    if fkey in photo and photo[fkey] not in filter_on[fkey]:
                        skip = True
                if skip:
                    continue
                self._add_cat_dict(
                    Photometry,
                    self._KEYS.PHOTOMETRY,
                    compare_to_existing=compare_to_existing,
                    **photo)
                phcount += 1
            self._log.debug("Added {} '{}' entries".format(
                phcount, photo_key))

        # Handle `spectra`
        # ---------------
        spec_key = self._KEYS.SPECTRA
        if spec_key in data:
            # When we are cleaning internal data, we don't always want to
            # require all of the normal spectrum data elements.
            spectra = data.pop(spec_key)
            self._log.debug("Found {} '{}' entries".format(
                len(spectra), spec_key))
            for spec in spectra:
                self._add_cat_dict(
                    Spectrum,
                    self._KEYS.SPECTRA,
                    compare_to_existing=compare_to_existing,
                    **spec)

        # Handle `error`
        # --------------
        err_key = self._KEYS.ERRORS
        if err_key in data:
            errors = data.pop(err_key)
            self._log.debug("Found {} '{}' entries".format(
                len(errors), err_key))
            for err in errors:
                self._add_cat_dict(Error, self._KEYS.ERRORS, **err)

        # Handle `models`
        # ---------------
        model_key = self._KEYS.MODELS
        if model_key in data:
            # When we are cleaning internal data, we don't always want to
            # require all of the normal spectrum data elements.
            model = data.pop(model_key)
            self._log.debug("Found {} '{}' entries".format(
                len(model), model_key))
            for mod in model:
                self._add_cat_dict(
                    Model, self._KEYS.MODELS, compare_to_existing=compare_to_existing, **mod)

        # Handle everything else --- should be `Quantity`s
        # ------------------------------------------------
        if len(data):
            self._log.debug("{} remaining entries, assuming `Quantity`".format(len(data)))
            # Iterate over remaining keys
            for key in list(data.keys()):
                vals = data.pop(key)
                # All quantities should be in lists of that quantity
                #    E.g. `aliases` is a list of alias quantities
                if not isinstance(vals, list):
                    vals = [vals]
                self._log.debug("{}: {}".format(key, vals))
                for vv in vals:
                    self._add_cat_dict(
                        Quantity, key,
                        check_for_dupes=merge, compare_to_existing=compare_to_existing, **vv)

        if merge and self.dupe_of:
            self.merge_dupes()

        return


class _Entry_New_Adder(_Entry):

    def add_alias(self, alias, source, clean=True, check_for_dupes=True, **kwargs):
        """Add an alias, optionally 'cleaning' the alias string.

        Calls the parent `catalog` method `clean_entry_name` - to apply the
        same name-cleaning as is applied to entry names themselves.

        Returns
        -------
        alias : str
            The stored version of the alias (cleaned or not).

        """
        if clean:
            alias = self.catalog.clean_entry_name(alias)

        # self.add_quantity(self._KEYS.ALIAS, alias, source)
        kwargs.update({QUANTITY.VALUE: alias, QUANTITY.SOURCE: source})
        '''
        new_data, retval = self.add_data(
            self._KEYS.ALIAS, check_for_dupes=check_for_dupes, **kwargs)
        '''

        # Check if this adding this alias makes us a dupe, if so mark ourselves as a dupe.
        if check_for_dupes and (alias in self.catalog.aliases):
            poss_dupe = self.catalog.aliases[alias]
            if (poss_dupe != self[self._KEYS.NAME]) and (poss_dupe in self.catalog.entries):
                self.dupe_of.append(poss_dupe)

        # If we're not checking for duplicates, warn about overwriting
        elif alias in self.catalog.aliases:
            self._log.warning("Warning, overwriting alias '{}': '{}'".format(
                alias, self.catalog.aliases[alias]))

        # Add this alias to parent catalog's reverse dictionary linking aliases to names
        self.catalog.aliases[alias] = self[self._KEYS.NAME]

        if self.dupe_of:
            self.merge_dupes()

        return alias

    def add_error(self, value, **kwargs):
        """Add an `Error` instance to this entry."""
        kwargs.update({ERROR.VALUE: value})
        # self._add_cat_dict(Error, self._KEYS.ERRORS, **kwargs)
        self.create_and_add_struct(Error, self._KEYS.ERRORS, **kwargs)
        return

    def add_photometry(self, **kwargs):
        """Add a `Photometry` instance to this entry."""
        # self._add_cat_dict(Photometry, self._KEYS.PHOTOMETRY,
        #                    compare_to_existing=compare_to_existing, **kwargs)
        val = kwargs.pop("compare_to_existing", None)
        if val is not None:
            if DEP_WARN_ARG:
                utils.log_deprecated_argument(
                    self._log, __file__,
                    "add_photometry", "compare_to_existing", "duplicate_check")
                kwargs["duplicate_check"] = val

        self.create_and_add_struct(Photometry, self._KEYS.PHOTOMETRY, **kwargs)
        return

    def add_quantity(self, name, value, source, **kwargs):
        """Add an `Quantity` instance to this entry."""
        success = True
        kwargs.update({QUANTITY.VALUE: value, QUANTITY.SOURCE: source})
        '''
        new_quantity = self._add_cat_dict(
            Quantity, name,
            compare_to_existing=compare_to_existing, check_for_dupes=check_for_dupes,
            **kwargs)
        '''

        val = kwargs.pop("compare_to_existing", None)
        if val is not None:
            if DEP_WARN_ARG:
                utils.log_deprecated_argument(
                    self._log, __file__,
                    "add_photometry", "compare_to_existing", "duplicate_check")
            kwargs["duplicate_check"] = val

        val = kwargs.pop("check_for_dupes", None)
        if val is not None:
            if DEP_WARN_ARG:
                utils.log_deprecated_argument(
                    self._log, __file__, "add_photometry", "check_for_dupes")

        new_quantity = self.create_and_add_struct(Quantity, name, **kwargs)

        if isinstance(new_quantity, Quantity):
            success = False
        elif new_quantity is False:
            success = False

        return success

    def add_source(self, allow_alias=False, **kwargs):
        """Add a `Source` instance to this entry."""
        log = self._log
        # log.debug("Entry.add_source()")
        if (not allow_alias) and (SOURCE.ALIAS in kwargs):
            err_str = "`{}` passed in kwargs, this shouldn't happen!".format(SOURCE.ALIAS)
            log.raise_error(err_str)

        # Set alias number to be +1 of current number of sources
        alias = kwargs.setdefault(SOURCE.ALIAS, str(self.num_sources() + 1))

        '''
        # Set alias number to be +1 of current number of sources
        if SOURCE.ALIAS not in kwargs:
            kwargs[SOURCE.ALIAS] = str(self.num_sources() + 1)

        source_obj = self._init_cat_dict(Source, self._KEYS.SOURCES, **kwargs)
        if source_obj is None:
            return None

        for item in self.get(self._KEYS.SOURCES, []):
            if source_obj.is_duplicate_of(item):
                return item[SOURCE.ALIAS]

        self.setdefault(self._KEYS.SOURCES, []).append(source_obj)
        return source_obj[SOURCE.ALIAS]
        '''

        # source = self.create_and_add_struct(Source, self._KEYS.SOURCES, duplicate_merge=False
        new_source = self.create_struct(Source, self._KEYS.SOURCES, **kwargs)
        if new_source is None:
            log.raise_error("Could not create_struct for source!")

        source, new_source = self.add_struct(self._KEYS.SOURCES, new_source,
                                             duplicate_check=True, duplicate_merge=False)
        alias = self._get_alias_from_add_struct_return(alias, source)

        return alias

    def add_model(self, allow_alias=False, **kwargs):
        """Add a `Model` instance to this entry."""
        log = self._log
        if not allow_alias and MODEL.ALIAS in kwargs:
            err_str = "`{}` passed in kwargs, this shouldn't happen!".format(SOURCE.ALIAS)
            log.raise_error(err_str)

        # Set alias number to be +1 of current number of models
        '''
        if MODEL.ALIAS not in kwargs:
            kwargs[MODEL.ALIAS] = str(self.num_models() + 1)
        model_obj = self._init_cat_dict(Model, self._KEYS.MODELS, **kwargs)
        if model_obj is None:
            return None

        for item in self.get(self._KEYS.MODELS, ''):
            if model_obj.is_duplicate_of(item):
                return item[item._KEYS.ALIAS]

        self.setdefault(self._KEYS.MODELS, []).append(model_obj)
        return model_obj[model_obj._KEYS.ALIAS]
        '''

        alias = kwargs.setdefault(MODEL.ALIAS, str(self.num_sources() + 1))
        new_model = self.create_struct(Source, self._KEYS.MODELS, **kwargs)
        if new_model is None:
            log.raise_error("Could not create_struct for model!")

        model, new_model = self.add_struct(self._KEYS.MODELS, new_model,
                                           duplicate_check=True, duplicate_merge=False)

        alias = self._get_alias_from_add_struct_return(alias, model)
        return alias

    def add_spectrum(self, compare_to_existing=True, **kwargs):
        """Add a `Spectrum` instance to this entry."""

        '''
        # Make sure that a source is given, and is valid (nor erroneous)
        retval = self._check_cat_dict_source(Spectrum, spec_key, **kwargs)
        if not retval:
            return None

        # Try to create a new instance of `Spectrum`
        new_spectrum = self._init_cat_dict(Spectrum, spec_key, **kwargs)
        if new_spectrum is None:
            return None

        is_dupe = False
        for item in self.get(spec_key, []):
            # Only the `filename` should be compared for duplicates. If a
            # duplicate is found, that means the previous `exclude` array
            # should be saved to the new object, and the old deleted
            if new_spectrum.is_duplicate_of(item):
                if SPECTRUM.EXCLUDE in new_spectrum:
                    item[SPECTRUM.EXCLUDE] = new_spectrum[SPECTRUM.EXCLUDE]
                elif SPECTRUM.EXCLUDE in item:
                    item.update(new_spectrum)
                is_dupe = True
                break

        if not is_dupe:
            self.setdefault(spec_key, []).append(new_spectrum)
        '''

        val = kwargs.pop("compare_to_existing", None)
        if val is not None:
            if DEP_WARN_ARG:
                utils.log_deprecated_argument(
                    self._log, __file__, "add_spectrum", "compare_to_existing", "duplicate_check")
            kwargs["duplicate_check"] = val

        self.create_and_add_struct(Spectrum, self._KEYS.SPECTRA, **kwargs)

        return

    '''
    def add_listed(self, key, value, check=True):
        listed = self.setdefault(key, [])
        # Make sure `value` isn't already in the list
        if check and (value in listed):
            return

        listed.append(value)
        return
    '''

    def add_self_source(self):
        """Add a source that refers to the catalog itself.

        For now this points to the Open Supernova Catalog by default.
        """
        return self.add_source(
            bibcode=self.catalog.OSC_BIBCODE,
            name=self.catalog.OSC_NAME,
            url=self.catalog.OSC_URL,
            secondary=True)

    def _handle_addition_failure(self, fail_loc, cat_class, cat_key, **kwargs):
        """Based on `catalog.ADDITION_FAILURE_BEHAVIOR`, react to a failure appropriately.

        A `logging.DEBUG` level message is always given.

        `ADDITION_FAILURE_BEHAVIOR` == `ADD_FAIL_ACTION.WARN`
            Then a `logging.WARNING` message is raised.
        `ADDITION_FAILURE_BEHAVIOR` == `ADD_FAIL_ACTION.IGNORE`
            No addition action is taken.
        `ADDITION_FAILURE_BEHAVIOR` == `ADD_FAIL_ACTION.RAISE`
            Then an error is raised.
            This is the default behavior that also acts if an unknown value is given.

        """
        err_str = "'{}' failed!\n".format(fail_loc)
        err_str += "class: '{}', key: '{}'\nkwargs: '{}'".format(cat_class, cat_key, kwargs)

        fail_flag = self.catalog.ADDITION_FAILURE_BEHAVIOR
        # Log a message at 'debug' level regardless
        self._log.debug(err_str)
        self._log.debug("`ADDITION_FAILURE_BEHAVIOR` = '{}'".format(fail_flag))

        # if `WARN` then also log a warning
        if fail_flag == utils.ADD_FAIL_ACTION.WARN:
            self._log.warning(err_str)
        # Raise an error
        elif fail_flag == utils.ADD_FAIL_ACTION.IGNORE:
            pass
        # default behavior is to raise an error
        else:
            self._log.raise_error(err_str, RuntimeError)

        return
