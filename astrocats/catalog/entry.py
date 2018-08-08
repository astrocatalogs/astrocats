"""Definitions related to the `Entry` class for catalog entries."""
import codecs
import gzip as gz
import hashlib
import json
import logging
import os
import sys
from collections import OrderedDict
from copy import deepcopy
from decimal import Decimal

import astrocats
from astrocats.catalog import utils, struct
from astrocats.catalog.struct import (Error, Model, Photometry, Quantity, Source, Spectrum)
from astrocats.catalog.struct import (ERROR, MODEL, QUANTITY, SOURCE, SPECTRUM)
# from past.builtins import basestring
from six import string_types

PATH_SCHEMA_INPUT = os.path.join(astrocats._PATH_SCHEMA, "input", "")
PATH_SCHEMA_OUTPUT = os.path.join(astrocats._PATH_SCHEMA, "output", "")

_PAS_PATH = "/Users/lzkelley/Research/catalogs/astroschema"
if _PAS_PATH not in sys.path:
    sys.path.append(_PAS_PATH)

import pyastroschema as pas


path_my_entry_schema = os.path.join(PATH_SCHEMA_INPUT, "astrocats_entry.json")
@pas.struct.set_struct_schema("entry", extensions=[path_my_entry_schema])  # noqa
class Entry(struct.Meta_Struct):

    def __init__(self, catalog=None, name=None, stub=False):
        # NOTE: FIX: LZK: cannot validate until a valid `name` is set
        super(Entry, self).__init__(catalog, key=name, validate=False)
        self.catalog = catalog
        self.filename = None
        self.dupe_of = []
        self._stub = stub
        if catalog is not None:
            self._log = catalog.log
        else:
            from astrocats.catalog.catalog import Catalog
            self._log = logging.getLogger()
            self._log.error("WARNING: Entry created without a catalog... creating catalog!")
            self.catalog = Catalog(None, self._log)

        self[self._KEYS.NAME] = name

        # NOTE: FIX: LZK: cannot validate until a valid `name` is set
        # self.validate()
        return

    def __repr__(self):
        """Return JSON representation of self."""
        jsonstring = utils.dict_to_pretty_string({ENTRY.NAME: self})
        return jsonstring

    def _append_additional_tags(self, quantity, source, cat_dict):
        """Append additional bits of data to an existing quantity.

        Called when a newly added quantity is found to be a duplicate.
        """
        pass

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

    def get_hash(self, keys=[]):
        """Return a unique hash associated with the listed keys."""
        if not len(keys):
            keys = list(self.keys())

        string_rep = ''
        oself = self._ordered(deepcopy(self))
        for key in keys:
            string_rep += json.dumps(oself.get(key, ''), sort_keys=True)

        return hashlib.sha512(string_rep.encode()).hexdigest()[:16]

    def _clean_quantity(self, quantity):
        """Clean quantity value before it is added to entry."""
        value = quantity.get(QUANTITY.VALUE, '').strip()
        error = quantity.get(QUANTITY.E_VALUE, '').strip()
        unit = quantity.get(QUANTITY.U_VALUE, '').strip()
        kind = quantity.get(QUANTITY.KIND, '')

        if isinstance(kind, list) and not isinstance(kind, string_types):
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

    def __deepcopy__(self, memo):
        """Define how an `Entry` should be deep copied."""
        new_entry = self.__class__(self.catalog)
        for key in self:
            if not key.startswith('__') and key != 'catalog':
                new_entry[key] = deepcopy(self[key], memo)
        return new_entry

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
        self._log.debug("Name: {}".format(name))

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
        self_name = self[ENTRY.NAME]
        if len(self_name) == 0:
            self[ENTRY.NAME] = name
        # Warn if there is a name mismatch
        elif self_name.lower().strip() != name.lower().strip():
            self._log.warning("Object name '{}' does not match name in json: '{}'".format(
                self_name, name))

        self.validate()
        return

    def _convert_odict_to_classes(self, data, clean=False, merge=True, pop_schema=True,
                                  compare_to_existing=True, filter_on={}):
        """Convert `OrderedDict` into `Entry` or its derivative classes."""
        self._log.debug("_convert_odict_to_classes(): {}".format(self.name()))
        self._log.debug("This should be a temporary fix.  Dont be lazy.")

        # Setup filters. Currently only used for photometry.
        fkeys = list(filter_on.keys())

        # Handle 'name'
        name_key = self._KEYS.NAME
        if name_key in data:
            self[name_key] = data.pop(name_key)

        # Handle 'schema'
        '''
        schema_key = self._KEYS.SCHEMA
        if schema_key in data:
            # Schema should be re-added every execution (done elsewhere) so
            # just delete the old entry
            if pop_schema:
                data.pop(schema_key)
            else:
                self[schema_key] = data.pop(schema_key)
        '''

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
                    Model,
                    self._KEYS.MODELS,
                    compare_to_existing=compare_to_existing,
                    **mod)

        # Handle everything else --- should be `Quantity`s
        # ------------------------------------------------
        if len(data):
            self._log.debug("{} remaining entries, assuming `Quantity`".format(
                len(data)))
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
                        Quantity,
                        key,
                        check_for_dupes=merge,
                        compare_to_existing=compare_to_existing,
                        **vv)

        if merge and self.dupe_of:
            self.merge_dupes()

        return

    def _check_cat_dict_source(self, cat_dict_class, key_in_self, **kwargs):
        """Check that the quantity isn't erroneous or private."""
        # Make sure that a source is given
        source = kwargs.get(cat_dict_class._KEYS.SOURCE, None)
        if source is not None:
            # If this source/data is erroneous, skip it
            if self.is_erroneous(key_in_self, source):
                self._log.info("This source is erroneous, skipping")
                return False
            # If this source/data is private, skip it
            if (self.catalog.args is not None and not self.catalog.args.private and
                    self.is_private(key_in_self, source)):
                self._log.info("This source is private, skipping")
                return False

        return True

    def _init_cat_dict(self, cat_dict_class, key_in_self, **kwargs):
        """Initialize a CatDict object, checking for errors."""

        # Remove empty string values
        bad_keys = [kk for kk, vv in kwargs.items()
                    if isinstance(vv, str) and len(vv) == 0]
        for bk in bad_keys:
            kwargs.pop(bk)

        # Catch errors associated with crappy, but not unexpected data
        try:
            new_entry = cat_dict_class(self, key=key_in_self, **kwargs)
        except struct.CatDictError as err:
            if err.warn:
                self._log.info("'{}' Not adding '{}': '{}'".format(
                    self[self._KEYS.NAME], key_in_self, str(err)))
            return None
        except Exception as err:
            log = self._log
            log.error("ERROR in `Entry._init_cat_dict()`!")
            log.error("Current task = '{}'".format(self._parent.current_task))
            log.error("name = '{}'".format(self[self._KEYS.NAME]))
            log.error("\n\n\n")
            log.error("`cat_dict_class` = '{}'".format(cat_dict_class))
            log.error("`key_in_self` = '{}'".format(key_in_self))
            log.error("`kwargs = '{}'".format(kwargs))
            log.error("\n\n")
            log.error("self = \n'{}'".format(repr(self)))
            log.error("\n")
            raise

        # NOTE: LZK COMMENTED OUT for testing new `Quantity` astroschema
        '''
        try:
            new_entry = cat_dict_class(self, key=key_in_self, **kwargs)
        except struct.CatDictError as err:
            if err.warn:
                self._log.info("'{}' Not adding '{}': '{}'".format(self[
                    self._KEYS.NAME], key_in_self, str(err)))
            return None
        '''

        return new_entry

    def _add_cat_dict(self, cat_dict_class, key_in_self,
                      check_for_dupes=True, compare_to_existing=True, **kwargs):
        """Add a `CatDict` to this `Entry`.

        CatDict only added if initialization succeeds and it
        doesn't already exist within the Entry.
        """
        # Make sure that a source is given, and is valid (nor erroneous)
        if cat_dict_class != Error:
            try:
                source = self._check_cat_dict_source(cat_dict_class, key_in_self, **kwargs)
            except struct.CatDictError as err:
                if err.warn:
                    msg = "'{}' Not adding '{}': '{}'".format(
                        self[self._KEYS.NAME], key_in_self, str(err))
                    self._log.info(msg)
                return False

            if source is None:
                return False

        # Try to create a new instance of this subclass of `CatDict`
        new_entry = self._init_cat_dict(cat_dict_class, key_in_self, **kwargs)
        if new_entry is None:
            if self.catalog.RAISE_ERROR_ON_ADDITION_FAILURE:
                err_str = "Entry._init_cat_dict() failed!"
                err_str += "class: '{}', key_in_self: '{}', kwargs: '{}'".format(
                    cat_dict_class, key_in_self, kwargs)
                utils.log_raise(self._log, err_str, RuntimeError)
            return False

        # Compare this new entry with all previous entries to make sure is new
        #    If it is NOT new, return the entry
        if compare_to_existing and cat_dict_class != Error:
            for item in self.get(key_in_self, []):
                if new_entry.is_duplicate_of(item):
                    item.append_sources_from(new_entry)
                    # Return the entry in case we want to use any additional
                    # tags to augment the old entry
                    return new_entry

        self.setdefault(key_in_self, []).append(new_entry)

        # NOTE: below code moved to `add_alias`
        '''
        # If this is an alias, add it to the parent catalog's reverse
        # dictionary linking aliases to names for fast lookup.
        if key_in_self == self._KEYS.ALIAS:
            # Check if this adding this alias makes us a dupe, if so mark
            # ourselves as a dupe.
            if (check_for_dupes and 'aliases' in dir(self.catalog) and
                    new_entry[QUANTITY.VALUE] in self.catalog.aliases):
                possible_dupe = self.catalog.aliases[new_entry[QUANTITY.VALUE]]
                if ((possible_dupe != self[self._KEYS.NAME]) and
                        (possible_dupe in self.catalog.entries)):
                    self.dupe_of.append(possible_dupe)

            if 'aliases' in dir(self.catalog):
                self.catalog.aliases[new_entry[QUANTITY.VALUE]] = self[self._KEYS.NAME]

        self.setdefault(key_in_self, []).append(new_entry)

        if (key_in_self == self._KEYS.ALIAS and check_for_dupes and self.dupe_of):
            self.merge_dupes()
        '''

        return True

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
            utils.log_raise(self._log, err_str, RuntimeError)

        return

    @classmethod
    def init_from_file(cls, catalog, name=None, path=None, try_gzip=False, **kwargs):
        """Construct a new `Entry` instance from an input file.

        The input file can be given explicitly by `path`, or a path will
        be constructed appropriately if possible.

        Arguments
        ---------
        catalog : `astrocats.catalog.catalog.Catalog` instance
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
            from astrocats.catalog.catalog import Catalog
            log = logging.getLogger()
            catalog = Catalog(None, log)

        catalog.log.debug("init_from_file()")
        if name is None and path is None:
            err = ("Either entry `name` or `path` must be specified to load entry.")
            log.error(err)
            raise ValueError(err)

        # If the path is given, use that to load from
        load_path = ''
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

    '''
    def add(self, key, value, **kwargs):
        log = self._log
        log.debug("Entry.add()")

        return
    '''

    def add_data(self, key_in_self, value=None, source=None, cat_dict_class=Quantity,
                 check_for_dupes=True, dupes_merge_tags=True, **kwargs):
        """Add a `CatDict` data container (dict) to this `Entry`.

        CatDict only added if initialization succeeds and it
        doesn't already exist within the Entry.
        """
        log = self._log
        log.debug("Entry.add_data()")
        FAIL = 0
        DUPE = -1
        SUCC = 1

        if value is not None:
            if QUANTITY.VALUE in kwargs:
                utils.log_raise(log, "`value` given as both arg and kwarg!\n{}".format(kwargs))
            kwargs[QUANTITY.VALUE] = value

        if source is not None:
            if QUANTITY.SOURCE in kwargs:
                utils.log_raise(log, "`source` given as both arg and kwarg!\n{}".format(kwargs))
            kwargs[QUANTITY.SOURCE] = source

        # Make sure that a source, if given, is valid
        retval = self._check_cat_dict_source(cat_dict_class, key_in_self, **kwargs)
        if not retval:
            self._handle_addition_failure(
                "Entry._check_cat_dict_source()", cat_dict_class, key_in_self, **kwargs)
            return None, FAIL

        # Try to create a new instance of this subclass of `CatDict`
        new_data = self._init_cat_dict(cat_dict_class, key_in_self, **kwargs)
        if new_data is None:
            self._handle_addition_failure(
                "Entry._init_cat_dict()", cat_dict_class, key_in_self, **kwargs)
            return None, FAIL

        # Compare this new entry with all previous entries to make sure is new
        if check_for_dupes:
            for item in self.get(key_in_self, []):
                # Do not add duplicates
                if new_data.is_duplicate_of(item):
                    item.append_sources_from(new_data)
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
        if (not allow_alias) and (SOURCE.ALIAS in kwargs):
            err_str = "`{}` passed in kwargs, this shouldn't happen!".format(SOURCE.ALIAS)
            self._log.error(err_str)
            raise RuntimeError(err_str)

        # Set alias number to be +1 of current number of sources
        if SOURCE.ALIAS not in kwargs:
            kwargs[SOURCE.ALIAS] = str(self.num_sources() + 1)
        source_obj = self._init_cat_dict(Source, self._KEYS.SOURCES, **kwargs)
        if source_obj is None:
            return None

        for item in self.get(self._KEYS.SOURCES, []):
            if source_obj.is_duplicate_of(item):
                return item[item._KEYS.ALIAS]

        self.setdefault(self._KEYS.SOURCES, []).append(source_obj)
        return source_obj[source_obj._KEYS.ALIAS]

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

    def add_self_source(self):
        """Add a source that refers to the catalog itself.

        For now this points to the Open Supernova Catalog by default.
        """
        return self.add_source(
            bibcode=self.catalog.OSC_BIBCODE,
            name=self.catalog.OSC_NAME,
            url=self.catalog.OSC_URL,
            secondary=True)

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

    def get_entry_text(self, fname):
        """Retrieve the raw text from a file."""
        if fname.split('.')[-1] == 'gz':
            with gz.open(fname, 'rt') as f:
                filetext = f.read()
        else:
            with codecs.open(fname, 'r') as f:
                filetext = f.read()
        return filetext

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
        source : `astrocats.catalog.source.Source` object
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
        stub : `astrocats.catalog.entry.Entry` subclass object
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
        if key == ENTRY.ALIAS:
            return False
        return all([SOURCE.PRIVATE in self.get_source_by_alias(x)
                    for x in sources.split(',')])

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

    def num_models(self):
        """Return the current number of models stored in this instance.

        Returns
        -------
        len : int
            The *integer* number of existing models.
        """
        return len(self.get(self._KEYS.MODELS, []))

    def priority_prefixes(self):
        """Return prefixes to given priority when merging duplicate entries."""
        return ()

    def sanitize(self):
        """Sanitize the data (sort it, etc.) before writing it to disk.

        Template method that can be overridden in each catalog's subclassed
        `Entry` object.
        """
        name = self[self._KEYS.NAME]

        aliases = self.get_aliases(includename=False)
        if name not in aliases:
            # Assign the first source to alias, if not available assign us.
            if self._KEYS.SOURCES in self:
                self.add_quantity(self._KEYS.ALIAS, name, '1')
                if self._KEYS.ALIAS not in self:
                    source = self.add_self_source()
                    self.add_quantity(self._KEYS.ALIAS, name, source)
            else:
                source = self.add_self_source()
                self.add_quantity(self._KEYS.ALIAS, name, source)

        if self._KEYS.ALIAS in self:
            self[self._KEYS.ALIAS].sort(key=lambda key: utils.alias_priority(name, key[QUANTITY.VALUE]))
        else:
            self._log.error('There should be at least one alias for `{}`.'.format(name))

        if self.catalog.name != 'BlackholeCatalog':
            self._log.error(self.catalog.name)
            self._log.error("WARNING: the `PHOTOMETRY` and `SPECTRA` "
                            "portion of sanitize have been removed!")
        '''
        if self._KEYS.PHOTOMETRY in self:
            self[self._KEYS.PHOTOMETRY].sort(
                key=lambda x: ((float(x[PHOTOMETRY.TIME]) if
                                isinstance(x[PHOTOMETRY.TIME],
                                           (basestring, float, int))
                                else min([float(y) for y in
                                          x[PHOTOMETRY.TIME]])) if
                               PHOTOMETRY.TIME in x else 0.0,
                               x[PHOTOMETRY.BAND] if PHOTOMETRY.BAND in
                               x else '',
                               float(x[PHOTOMETRY.MAGNITUDE]) if
                               PHOTOMETRY.MAGNITUDE in x else ''))

        if (self._KEYS.SPECTRA in self and list(
                filter(None, [SPECTRUM.TIME in x for x in self[self._KEYS.SPECTRA]]))):
            self[self._KEYS.SPECTRA].sort(
                key=lambda x: (float(x[SPECTRUM.TIME]) if
                               SPECTRUM.TIME in x else 0.0,
                               x[SPECTRUM.FILENAME] if
                               SPECTRUM.FILENAME in x else '')
            )
        '''

        if self._KEYS.SOURCES in self:
            # _NO_SRCS_REQUIRED = [self._KEYS.NAME, self._KEYS.SCHEMA, self._KEYS.SOURCES,
            #                      self._KEYS.ERRORS]
            # Remove orphan sources
            source_aliases = [x[SOURCE.ALIAS] for x in self[self._KEYS.SOURCES]]
            # Sources with the `PRIVATE` attribute are always retained
            source_list = [
                x[SOURCE.ALIAS] for x in self[self._KEYS.SOURCES]
                if SOURCE.PRIVATE in x
            ]
            # _SOURCES_SKIP_KEYS = [
            #     self._KEYS.NAME, self._KEYS.SCHEMA, self._KEYS.SOURCES, self._KEYS.ERRORS]
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

        # FIX: use 'dump' not 'dumps'
        jsonstring = json.dumps({self[self._KEYS.NAME]: self._ordered(self)},
                                indent='\t' if sys.version_info[0] >= 3 else 4,
                                separators=(',', ':'),
                                ensure_ascii=False)
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


ENTRY = Entry._KEYCHAIN
Entry._KEYS = ENTRY

struct.output_schema(PATH_SCHEMA_OUTPUT, Entry)
