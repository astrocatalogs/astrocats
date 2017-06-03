"""Definitions related to the `Entry` class for catalog entries."""
import codecs
import hashlib
import json
import logging
import os
import sys
import gzip as gz
from collections import OrderedDict
from copy import deepcopy

from astrocats.catalog.catdict import CatDict, CatDictError
from astrocats.catalog.error import ERROR, Error
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection
from astrocats.catalog.model import MODEL, Model
from astrocats.catalog.photometry import PHOTOMETRY, Photometry
from astrocats.catalog.quantity import QUANTITY, Quantity
from astrocats.catalog.source import SOURCE, Source
from astrocats.catalog.spectrum import SPECTRUM, Spectrum
from astrocats.catalog.utils import (alias_priority, dict_to_pretty_string,
                                     is_integer, is_number, listify)
from decimal import Decimal
from past.builtins import basestring


class ENTRY(KeyCollection):
    """General `CatDict` keys which should be relevant for all catalogs."""

    # Constants for use in key definitions
    _DIST_PREF_KINDS = [
        'heliocentric', 'cmb', 'spectroscopic', 'photometric', 'host',
        'cluster'
    ]
    _HOST_DIST_PREF_KINDS = [
        'heliocentric', 'cmb', 'spectroscopic', 'photometric', 'host',
        'cluster'
    ]

    # List of keys
    ALIAS = Key('alias', KEY_TYPES.STRING)
    COMOVING_DIST = Key('comovingdist',
                        KEY_TYPES.NUMERIC,
                        kind_preference=_DIST_PREF_KINDS,
                        replace_better=True)
    DEC = Key('dec', KEY_TYPES.STRING)
    DISCOVER_DATE = Key('discoverdate', KEY_TYPES.STRING, replace_better=True)
    DISCOVERER = Key('discoverer', KEY_TYPES.STRING)
    DISTINCT_FROM = Key('distinctfrom', KEY_TYPES.STRING)
    EBV = Key('ebv', KEY_TYPES.NUMERIC, replace_better=True)
    AV_CIRCUM = Key('avcircum', KEY_TYPES.NUMERIC, replace_better=True)
    ERRORS = Key('errors', no_source=True)
    HOST = Key('host', KEY_TYPES.STRING)
    HOST_DEC = Key('hostdec', KEY_TYPES.STRING)
    HOST_OFFSET_ANG = Key('hostoffsetang', KEY_TYPES.NUMERIC)
    HOST_OFFSET_DIST = Key('hostoffsetdist', KEY_TYPES.NUMERIC)
    HOST_RA = Key('hostra', KEY_TYPES.STRING)
    HOST_REDSHIFT = Key('hostredshift',
                        KEY_TYPES.NUMERIC,
                        kind_preference=_HOST_DIST_PREF_KINDS,
                        replace_better=True)
    HOST_VELOCITY = Key('hostvelocity',
                        KEY_TYPES.NUMERIC,
                        kind_preference=_HOST_DIST_PREF_KINDS,
                        replace_better=True)
    HOST_LUM_DIST = Key('hostlumdist',
                        KEY_TYPES.NUMERIC,
                        kind_preference=_HOST_DIST_PREF_KINDS,
                        replace_better=True)
    HOST_COMOVING_DIST = Key('hostcomovingdist',
                             KEY_TYPES.NUMERIC,
                             kind_preference=_HOST_DIST_PREF_KINDS,
                             replace_better=True)
    LUM_DIST = Key('lumdist',
                   KEY_TYPES.NUMERIC,
                   kind_preference=_DIST_PREF_KINDS,
                   replace_better=True)
    MAX_ABS_MAG = Key('maxabsmag', KEY_TYPES.NUMERIC)
    MAX_APP_MAG = Key('maxappmag', KEY_TYPES.NUMERIC)
    MAX_BAND = Key('maxband', KEY_TYPES.STRING)
    MAX_DATE = Key('maxdate', KEY_TYPES.STRING, replace_better=True)
    MODELS = Key('models')
    NAME = Key('name', KEY_TYPES.STRING, no_source=True)
    PHOTOMETRY = Key('photometry')
    RA = Key('ra', KEY_TYPES.STRING)
    REDSHIFT = Key('redshift',
                   KEY_TYPES.NUMERIC,
                   kind_preference=_DIST_PREF_KINDS,
                   replace_better=True)
    SCHEMA = Key('schema', no_source=True)
    SOURCES = Key('sources', no_source=True)
    SPECTRA = Key('spectra')
    VELOCITY = Key('velocity',
                   KEY_TYPES.NUMERIC,
                   kind_preference=_DIST_PREF_KINDS,
                   replace_better=True)


class Entry(OrderedDict):
    """Class representing an individual element of each Catalog.

    For example, a single supernova in the supernova catalog, this object
    handles and manages the addition of data for this `Entry`, using different
    `CatDict` instances (e.g. `Photometry`).

    Notes
    -----
    -   Stubs: a stub is the most minimal entry, containing an entry's 'name'
        and possible aliases.  These instances are used to represent entries
        which are known to exist (e.g. have already been saved) for cross
        referencing and duplicate removal.
        +   The `Entry.get_stub` method returns the 'stub' corresponding to the
            Entry instance.  i.e. it returns a *new object* with only the name
            and aliases copied over.

    Attributes
    ----------
    catalog : `astrocats.catalog.catalog.Catalog` object
        Pointer to the parent catalog object of which this entry is a member.
    filename : str or 'None'
        If this entry is loaded from a file, its (full path and) filename.
    _log : `logging.Logger` object
        Pointer to the logger from the parent catalog.
    _stub : bool
        Whether this instance represents a 'stub' (see above).
    _KEYS : `astrocats.catalog.key.KeyCollection` object
        The associated object which contains the different dictionary keys
        used in this type (e.g. `Supernova`) entry.

    """

    _KEYS = ENTRY

    def __init__(self, catalog=None, name=None, stub=False):
        """Create a new `Entry` object with the given `name`.

        Arguments
        ---------
        catalog : `astrocats.catalog.catalog.Catalog` instance
            The parent catalog object of which this entry belongs.
        name : str
            The name of this entry, e.g. `SN1987A` for a `Supernova` entry.
        stub : bool
            Whether or not this instance represents a 'stub' (see above).

        """
        super(Entry, self).__init__()
        self.catalog = catalog
        self.filename = None
        self.dupe_of = []
        self._stub = stub
        if catalog:
            self._log = catalog.log
        else:
            from astrocats.catalog.catalog import Catalog
            self._log = logging.getLogger()
            self.catalog = Catalog(None, self._log)
        self[self._KEYS.NAME] = name
        return

    def __repr__(self):
        """Return JSON representation of self."""
        jsonstring = dict_to_pretty_string({ENTRY.NAME: self})
        return jsonstring

    def _append_additional_tags(self, quantity, source, cat_dict):
        """Append additional bits of data to an existing quantity.

        Called when a newly added quantity is found to be a duplicate.
        """
        pass

    def _get_save_path(self, bury=False):
        """Return the path that this Entry should be saved to."""
        filename = self.get_filename(self[self._KEYS.NAME])

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

        if isinstance(odict, CatDict) or isinstance(odict, Entry):
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
        kind = quantity.get(QUANTITY.KIND, '').strip()

        if not value:
            return False

        if is_number(value):
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
                new_entry[key] = deepcopy(self[key])
        return new_entry

    def _load_data_from_json(self,
                             fhand,
                             clean=False,
                             merge=True,
                             pop_schema=True,
                             ignore_keys=[],
                             compare_to_existing=True,
                             gzip=False):
        # FIX: check for overwrite??"""
        self._log.debug("_load_data_from_json(): {}\n\t{}".format(self.name(),
                                                                  fhand))
        # Store the filename this was loaded from
        self.filename = fhand

        if gzip:
            jfil = gz.open(fhand, 'rb')
        else:
            jfil = open(fhand, 'r')

        data = json.load(jfil, object_pairs_hook=OrderedDict)
        name = list(data.keys())
        if len(name) != 1:
            err = "json file '{}' has multiple keys: {}".format(fhand,
                                                                list(name))
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
            data,
            clean=clean,
            merge=merge,
            pop_schema=pop_schema,
            compare_to_existing=compare_to_existing)
        if len(data):
            err_str = ("Remaining entries in `data` after "
                       "`_convert_odict_to_classes`.")
            err_str += "\n{}".format(dict_to_pretty_string(data))
            self._log.error(err_str)
            raise RuntimeError(err_str)

        jfil.close()

        # If object doesnt have a name yet, but json does, store it
        self_name = self[ENTRY.NAME]
        if len(self_name) == 0:
            self[ENTRY.NAME] = name
        # Warn if there is a name mismatch
        elif self_name.lower().strip() != name.lower().strip():
            self._log.warning("Object name '{}' does not match name in json:"
                              "'{}'".format(self_name, name))

        self.check()
        return

    def _convert_odict_to_classes(self,
                                  data,
                                  clean=False,
                                  merge=True,
                                  pop_schema=True,
                                  compare_to_existing=True):
        """Convert `OrderedDict` into `Entry` or its derivative classes."""
        self._log.debug("_convert_odict_to_classes(): {}".format(self.name()))
        self._log.debug("This should be a temporary fix.  Dont be lazy.")

        # Handle 'name'
        name_key = self._KEYS.NAME
        if name_key in data:
            self[name_key] = data.pop(name_key)

        # Handle 'schema'
        schema_key = self._KEYS.SCHEMA
        if schema_key in data:
            # Schema should be re-added every execution (done elsewhere) so
            # just delete the old entry
            if pop_schema:
                data.pop(schema_key)
            else:
                self[schema_key] = data.pop(schema_key)

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
            for photo in photoms:
                self._add_cat_dict(
                    Photometry,
                    self._KEYS.PHOTOMETRY,
                    compare_to_existing=compare_to_existing,
                    **photo)

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
        """Check that a source exists and that a quantity isn't erroneous."""
        # Make sure that a source is given
        source = kwargs.get(cat_dict_class._KEYS.SOURCE, None)
        if source is None:
            raise CatDictError(
                "{}: `source` must be provided!".format(self[self._KEYS.NAME]),
                warn=True)
        # Check that source is a list of integers
        for x in source.split(','):
            if not is_integer(x):
                raise CatDictError(
                    "{}: `source` is comma-delimited list of "
                    " integers!".format(self[self._KEYS.NAME]),
                    warn=True)
        # If this source/data is erroneous, skip it
        if self.is_erroneous(key_in_self, source):
            self._log.info("This source is erroneous, skipping")
            return None
        # If this source/data is private, skip it
        if (self.catalog.args is not None and not self.catalog.args.private and
                self.is_private(key_in_self, source)):
            self._log.info("This source is private, skipping")
            return None
        return source

    def _init_cat_dict(self, cat_dict_class, key_in_self, **kwargs):
        """Initialize a CatDict object, checking for errors."""
        # Catch errors associated with crappy, but not unexpected data
        try:
            new_entry = cat_dict_class(self, key=key_in_self, **kwargs)
        except CatDictError as err:
            if err.warn:
                self._log.info("'{}' Not adding '{}': '{}'".format(self[
                    self._KEYS.NAME], key_in_self, str(err)))
            return None
        return new_entry

    def _add_cat_dict(self,
                      cat_dict_class,
                      key_in_self,
                      check_for_dupes=True,
                      compare_to_existing=True,
                      **kwargs):
        """Add a `CatDict` to this `Entry`.

        CatDict only added if initialization succeeds and it
        doesn't already exist within the Entry.
        """
        # Make sure that a source is given, and is valid (nor erroneous)
        if cat_dict_class != Error:
            try:
                source = self._check_cat_dict_source(cat_dict_class,
                                                     key_in_self, **kwargs)
            except CatDictError as err:
                if err.warn:
                    self._log.info("'{}' Not adding '{}': '{}'".format(self[
                        self._KEYS.NAME], key_in_self, str(err)))
                return False

            if source is None:
                return False

        # Try to create a new instance of this subclass of `CatDict`
        new_entry = self._init_cat_dict(cat_dict_class, key_in_self, **kwargs)
        if new_entry is None:
            return False

        # Compare this new entry with all previous entries to make sure is new
        if compare_to_existing and cat_dict_class != Error:
            for item in self.get(key_in_self, []):
                if new_entry.is_duplicate_of(item):
                    item.append_sources_from(new_entry)
                    # Return the entry in case we want to use any additional
                    # tags to augment the old entry
                    return new_entry

        # If this is an alias, add it to the parent catalog's reverse
        # dictionary linking aliases to names for fast lookup.
        if key_in_self == self._KEYS.ALIAS:
            # Check if this adding this alias makes us a dupe, if so mark
            # ourselves as a dupe.
            if (check_for_dupes and 'aliases' in dir(self.catalog) and
                    new_entry[QUANTITY.VALUE] in self.catalog.aliases):
                possible_dupe = self.catalog.aliases[new_entry[QUANTITY.VALUE]]
                # print(possible_dupe)
                if (possible_dupe != self[self._KEYS.NAME] and
                        possible_dupe in self.catalog.entries):
                    self.dupe_of.append(possible_dupe)
            if 'aliases' in dir(self.catalog):
                self.catalog.aliases[new_entry[QUANTITY.VALUE]] = self[
                    self._KEYS.NAME]

        self.setdefault(key_in_self, []).append(new_entry)

        if (key_in_self == self._KEYS.ALIAS and check_for_dupes and
                self.dupe_of):
            self.merge_dupes()

        return True

    @classmethod
    def get_filename(cls, name):
        """Convert from an `Entry` name into an appropriate filename."""
        fname = name.replace('/', '_')
        return fname

    @classmethod
    def init_from_file(cls,
                       catalog,
                       name=None,
                       path=None,
                       clean=False,
                       merge=True,
                       pop_schema=True,
                       ignore_keys=[],
                       compare_to_existing=True,
                       try_gzip=False):
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
        clean : bool
            Whether special sanitization processing should be done on the input
            data.  This is mostly for input files from the 'internal'
            repositories.

        """
        if not catalog:
            from astrocats.catalog.catalog import Catalog
            log = logging.getLogger()
            catalog = Catalog(None, log)

        catalog.log.debug("init_from_file()")
        if name is None and path is None:
            err = ("Either entry `name` or `path` must be specified to load "
                   "entry.")
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
                filename = cls.get_filename(name)
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
        new_entry._load_data_from_json(
            load_path,
            clean=clean,
            merge=merge,
            pop_schema=pop_schema,
            ignore_keys=ignore_keys,
            compare_to_existing=compare_to_existing,
            gzip=try_gzip)

        return new_entry

    def add_alias(self, alias, source, clean=True):
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
        self.add_quantity(self._KEYS.ALIAS, alias, source)
        return alias

    def add_error(self, value, **kwargs):
        """Add an `Error` instance to this entry."""
        kwargs.update({ERROR.VALUE: value})
        self._add_cat_dict(Error, self._KEYS.ERRORS, **kwargs)
        return

    def add_photometry(self, compare_to_existing=True, **kwargs):
        """Add a `Photometry` instance to this entry."""
        self._add_cat_dict(
            Photometry,
            self._KEYS.PHOTOMETRY,
            compare_to_existing=compare_to_existing,
            **kwargs)
        return

    def merge_dupes(self):
        """Merge two entries that correspond to the same entry."""
        for dupe in self.dupe_of:
            if dupe in self.catalog.entries:
                if self.catalog.entries[dupe]._stub:
                    # merge = False to avoid infinite recursion
                    self.catalog.load_entry_from_name(
                        dupe, delete=True, merge=False)
                self.catalog.copy_entry_to_entry(self.catalog.entries[dupe],
                                                 self)
                del self.catalog.entries[dupe]
        self.dupe_of = []

    def add_quantity(self,
                     quantities,
                     value,
                     source,
                     check_for_dupes=True,
                     compare_to_existing=True,
                     **kwargs):
        """Add an `Quantity` instance to this entry."""
        success = True
        for quantity in listify(quantities):
            kwargs.update({QUANTITY.VALUE: value, QUANTITY.SOURCE: source})
            cat_dict = self._add_cat_dict(
                Quantity,
                quantity,
                compare_to_existing=compare_to_existing,
                check_for_dupes=check_for_dupes,
                **kwargs)
            if isinstance(cat_dict, CatDict):
                self._append_additional_tags(quantity, source, cat_dict)
                success = False

        return success

    def add_self_source(self):
        """Add a source that refers to the catalog itself.

        For now this points to the Open Supernova Catalog by default.
        """
        return self.add_source(
            bibcode=self.catalog.OSC_BIBCODE,
            name=self.catalog.OSC_NAME,
            url=self.catalog.OSC_URL,
            secondary=True)

    def add_source(self, allow_alias=False, **kwargs):
        """Add a `Source` instance to this entry."""
        if not allow_alias and SOURCE.ALIAS in kwargs:
            err_str = "`{}` passed in kwargs, this shouldn't happen!".format(
                SOURCE.ALIAS)
            self._log.error(err_str)
            raise RuntimeError(err_str)

        # Set alias number to be +1 of current number of sources
        if SOURCE.ALIAS not in kwargs:
            kwargs[SOURCE.ALIAS] = str(self.num_sources() + 1)
        source_obj = self._init_cat_dict(Source, self._KEYS.SOURCES, **kwargs)
        if source_obj is None:
            return None

        for item in self.get(self._KEYS.SOURCES, ''):
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
        source = self._check_cat_dict_source(Spectrum, spec_key, **kwargs)
        if source is None:
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

    def check(self):
        """Check that the entry has the required fields."""
        # Make sure there is a schema key in dict
        if self._KEYS.SCHEMA not in self:
            self[self._KEYS.SCHEMA] = self.catalog.SCHEMA.URL
        # Make sure there is a name key in dict
        if (self._KEYS.NAME not in self or len(self[self._KEYS.NAME]) == 0):
            raise ValueError("Entry name is empty:\n\t{}".format(
                json.dumps(
                    self, indent=2)))
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

    def get_entry_text(self, fname):
        """Retrieve the raw text from a file."""
        if fname.split('.')[-1] == 'gz':
            with gz.open(fname, 'rt') as f:
                filetext = f.read()
        else:
            with open(fname, 'r') as f:
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
        raise ValueError("Source '{}': alias '{}' not found!".format(self[
            self._KEYS.NAME], alias))

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
        return stub

    def is_erroneous(self, field, sources):
        """Check if attribute has been marked as being erroneous."""
        if self._KEYS.ERRORS in self:
            my_errors = self[self._KEYS.ERRORS]
            for alias in sources.split(','):
                source = self.get_source_by_alias(alias)
                bib_err_values = [
                    err[ERROR.VALUE] for err in my_errors
                    if err[ERROR.KIND] == SOURCE.BIBCODE and
                    err[ERROR.EXTRA] == field
                ]
                if (SOURCE.BIBCODE in source and
                        source[SOURCE.BIBCODE] in bib_err_values):
                    return True

                name_err_values = [
                    err[ERROR.VALUE] for err in my_errors
                    if err[ERROR.KIND] == SOURCE.NAME and err[ERROR.EXTRA] ==
                    field
                ]
                if (SOURCE.NAME in source and
                        source[SOURCE.NAME] in name_err_values):
                    return True

        return False

    def is_private(self, key, sources):
        """Check if attribute is private."""
        # aliases are always public.
        if key == ENTRY.ALIAS:
            return False
        return all([
            SOURCE.PRIVATE in self.get_source_by_alias(x)
            for x in sources.split(',')
        ])

    def name(self):
        """Return own name."""
        try:
            return self[self._KEYS.NAME]
        except KeyError:
            return None

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
            self[self._KEYS.ALIAS].sort(
                key=lambda key: alias_priority(name, key[QUANTITY.VALUE]))
        else:
            self._log.error(
                'There should be at least one alias for `{}`.'.format(name))

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
                filter(None, [
                    SPECTRUM.TIME in x for x in self[self._KEYS.SPECTRA]
                ]))):
            self[self._KEYS.SPECTRA].sort(
                key=lambda x: (float(x[SPECTRUM.TIME]) if
                               SPECTRUM.TIME in x else 0.0,
                               x[SPECTRUM.FILENAME] if
                               SPECTRUM.FILENAME in x else '')
            )

        if self._KEYS.SOURCES in self:
            # Remove orphan sources
            source_aliases = [
                x[SOURCE.ALIAS] for x in self[self._KEYS.SOURCES]
            ]
            # Sources with the `PRIVATE` attribute are always retained
            source_list = [
                x[SOURCE.ALIAS] for x in self[self._KEYS.SOURCES]
                if SOURCE.PRIVATE in x
            ]
            for key in self:
                # if self._KEYS.get_key_by_name(key).no_source:
                if (key in [
                        self._KEYS.NAME, self._KEYS.SCHEMA, self._KEYS.SOURCES,
                        self._KEYS.ERRORS
                ]):
                    continue
                for item in self[key]:
                    source_list += item[item._KEYS.SOURCE].split(',')
            new_src_list = sorted(
                list(set(source_aliases).intersection(source_list)))
            new_sources = []
            for source in self[self._KEYS.SOURCES]:
                if source[SOURCE.ALIAS] in new_src_list:
                    new_sources.append(source)
                else:
                    self._log.info('Removing orphaned source from `{}`.'
                                   .format(name))

            if not new_sources:
                del self[self._KEYS.SOURCES]

            self[self._KEYS.SOURCES] = new_sources

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
        jsonstring = json.dumps(
            {
                self[self._KEYS.NAME]: self._ordered(self)
            },
            indent='\t' if sys.version_info[0] >= 3 else 4,
            separators=(',', ':'),
            ensure_ascii=False)
        if not os.path.isdir(outdir):
            raise RuntimeError("Output directory '{}' for event '{}' does "
                               "not exist.".format(outdir, self[
                                   self._KEYS.NAME]))
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
        if key == self._KEYS.SCHEMA:
            return 'aaa'
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
