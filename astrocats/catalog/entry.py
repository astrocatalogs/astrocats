"""
"""
import codecs
import json
import os
from collections import OrderedDict

from astrocats.catalog.catdict import CatDict, CatDictError
from astrocats.catalog.error import ERROR, Error
from astrocats.catalog.key import KeyCollection
from astrocats.catalog.photometry import Photometry
from astrocats.catalog.quantity import QUANTITY, Quantity
from astrocats.catalog.source import SOURCE, Source
from astrocats.catalog.spectrum import SPECTRUM, Spectrum
from astrocats.catalog.utils import (alias_priority, dict_to_pretty_string,
                                     get_event_filename)


class KEYS(KeyCollection):
    """General `CatDict` keys which should be relevant for all catalogs.
    """
    ALIAS = 'alias'
    BIBCODE = 'bibcode'
    COMOVING_DIST = 'comovingdist'
    DEC = 'dec'
    DISCOVER_DATE = 'discoverdate'
    DISCOVERER = 'discoverer'
    DISTINCT_FROM = 'distinctfrom'
    EBV = 'ebv'
    ERRORS = 'errors'
    HOST = 'host'
    HOST_DEC = 'hostdec'
    HOST_OFFSET_ANG = 'hostoffsetang'
    HOST_OFFSET_DIST = 'hostoffsetdist'
    HOST_RA = 'hostra'
    LUM_DIST = 'lumdist'
    MAX_ABS_MAG = 'maxabsmag'
    MAX_APP_MAG = 'maxappmag'
    MAX_BAND = 'maxband'
    MAX_DATE = 'maxdate'
    NAME = 'name'
    PHOTOMETRY = 'photometry'
    RA = 'ra'
    REDSHIFT = 'redshift'
    SCHEMA = 'schema'
    SOURCES = 'sources'
    SPECTRA = 'spectra'
    URL = 'url'
    VELOCITY = 'velocity'


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

    _KEYS = KEYS

    def __init__(self, catalog, name, stub=False):
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
        self.catalog = catalog
        self.filename = None
        self._log = catalog.log
        self._stub = stub
        self[self._KEYS.NAME] = name
        return

    @classmethod
    def init_from_file(cls, catalog, name=None, path=None, clean=False):
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
        catalog.log.debug("init_from_file()")
        if name is None and path is None:
            raise ValueError("Either entry `name` or `path` must be specified "
                             "to load entry.")
        if name is not None and path is not None:
            raise ValueError("Either entry `name` or `path` should be "
                             "specified, not both.")

        # If the path is given, use that to load from
        load_path = ''
        if path is not None:
            load_path = path
            name = ''
        # If the name is given, try to find a path for it
        else:
            repo_paths = catalog.PATHS.get_repo_output_folders()
            for rep in repo_paths:
                filename = get_event_filename(name)
                newpath = os.path.join(rep, filename + '.json')
                if os.path.isfile(newpath):
                    load_path = newpath
                    break

        if load_path is None or not os.path.isfile(load_path):
            # FIX: is this warning worthy?
            return None

        # Create a new `Entry` instance
        new_entry = cls(catalog, name)
        # Fill it with data from json file
        new_entry._load_data_from_json(load_path, clean=clean)

        return new_entry

    def __repr__(self):
        jsonstring = dict_to_pretty_string({self[KEYS.NAME]: self})
        return jsonstring

    def add_error(self, quantity, value, **kwargs):
        """Add an `Error` instance to this entry.
        """
        kwargs.update({ERROR.VALUE: value})
        self._add_cat_dict(Quantity, quantity, **kwargs)
        return

    def add_photometry(self, **kwargs):
        """Add a `Photometry` instance to this entry.
        """
        self._add_cat_dict(Photometry, self._KEYS.PHOTOMETRY, **kwargs)
        return

    def add_quantity(self, quantity, value, source, **kwargs):
        """Add an `Quantity` instance to this entry.
        """
        # Aliases not added if in DISTINCT_FROM
        if quantity == self._KEYS.ALIAS:
            value = self.clean_entry_name(value)
            for df in self.get(self._KEYS.DISTINCT_FROM, []):
                if value == df[QUANTITY.VALUE]:
                    return False

        kwargs.update({QUANTITY.VALUE: value, QUANTITY.SOURCE: source})
        cat_dict = self._add_cat_dict(Quantity, quantity, **kwargs)
        if isinstance(cat_dict, CatDict):
            self._append_additional_tags(quantity, source, cat_dict)
            return False
        elif cat_dict:
            return True

        return False

    def _append_additional_tags(self, quantity, source, cat_dict):
        pass

    def add_source(self, allow_alias=False, **kwargs):
        """Add a `Source` instance to this entry.
        """
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

        # Set 'alias' number to be one higher than existing length of sources
        source_obj[SOURCE.ALIAS] = str(
            len(self.get(self._KEYS.SOURCES, [])) + 1)

        self.setdefault(self._KEYS.SOURCES, []).append(source_obj)
        return source_obj[source_obj._KEYS.ALIAS]

    def add_self_source(self):
        return self.add_source(
            bibcode=self.catalog.OSC_BIBCODE,
            name=self.catalog.OSC_NAME,
            url=self.catalog.OSC_URL, secondary=True)

    def add_spectrum(self, waveunit='', fluxunit='', **kwargs):
        """Add an `Spectrum` instance to this entry.
        """
        kwargs.update({SPECTRUM.WAVE_UNIT: waveunit,
                       SPECTRUM.FLUX_UNIT: fluxunit})
        spec_key = self._KEYS.SPECTRA
        # Make sure that a source is given, and is valid (nor erroneous)
        source = self._check_cat_dict_source(Spectrum, spec_key, **kwargs)
        if source is None:
            return None

        # Try to create a new instance of `Spectrum`
        new_spectrum = self._init_cat_dict(
            Spectrum, spec_key, **kwargs)
        if new_spectrum is None:
            return None

        num_spec = len(self.get(spec_key, []))
        for si in range(num_spec):
            item = self[spec_key][si]
            # Only the `filename` should be compared for duplicates If a
            # duplicate is found, that means the previous `exclude` array
            # should be saved to the new object, and the old deleted
            if new_spectrum.is_duplicate_of(item):
                if SPECTRUM.EXCLUDE in item:
                    new_spectrum[SPECTRUM.EXCLUDE] = item[SPECTRUM.EXCLUDE]
                del self[spec_key][si]
                break

        self.setdefault(spec_key, []).append(new_spectrum)
        return

    def check(self):
        """Check that the entry has the required fields.
        """
        # Make sure there is a schema key in dict
        if self._KEYS.SCHEMA not in self.keys():
            self[self._KEYS.SCHEMA] = self.catalog.SCHEMA.URL
        # Make sure there is a name key in dict
        if (self._KEYS.NAME not in self.keys() or
                len(self[self._KEYS.NAME]) == 0):
            raise ValueError("Entry name is empty:\n\t{}".format(
                json.dumps(self, indent=2)))
        return

    def clean_entry_name(self, name):
        """Template method to clean/sanitize an entry name before setting it.

        Should be overridden appropriately in subclasses `Entry` objects.
        """
        return name

    def clean_internal(self, data=None):
        """Clean input from 'internal', human added data.

        This is used in the 'Entry.init_from_file' method.
        """
        return data

    def get_aliases(self, includename=True):
        """Retrieve the aliases of this object as a list of strings.

        Arguments
        ---------
        includename : bool
            Include the 'name' parameter in the list of aliases.
        """
        # empty list if doesnt exist
        alias_quanta = self.get(self._KEYS.ALIAS, [])
        aliases = [aq['value'] for aq in alias_quanta]
        if includename and self[self._KEYS.NAME] not in aliases:
            aliases = [self[self._KEYS.NAME]] + aliases
        return aliases

    def get_entry_text(fname):
        """Retrieve the raw text from a file.
        """
        import gzip
        if fname.split('.')[-1] == 'gz':
            with gzip.open(fname, 'rt') as f:
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
        raise ValueError(
            "Source '{}': alias '{}' not found!".format(
                self[self._KEYS.NAME], alias))

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
        if self._KEYS.ALIAS in self.keys():
            stub[self._KEYS.ALIAS] = self[self._KEYS.ALIAS]
        return stub

    def name(self):
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
            else:
                source = self.add_self_source()
                self.add_quantity(self._KEYS.ALIAS, name, source)

        self[self._KEYS.ALIAS] = list(
            sorted(self[self._KEYS.ALIAS],
                   key=lambda key: alias_priority(name, key[QUANTITY.VALUE])))

    def sort_func(self, key):
        if key == self._KEYS.SCHEMA:
            return 'aaa'
        if key == self._KEYS.NAME:
            return 'aab'
        if key == self._KEYS.SOURCES:
            return 'aac'
        if key == self._KEYS.ALIAS:
            return 'aad'
        if key == self._KEYS.PHOTOMETRY:
            return 'zzy'
        if key == self._KEYS.SPECTRA:
            return 'zzz'
        return key

    def _get_save_path(self, bury=False):
        raise RuntimeError("This method must be overridden!")

    def _ordered(self, odict):
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
                                indent='\t', separators=(',', ':'),
                                ensure_ascii=False)
        if not os.path.isdir(outdir):
            raise RuntimeError("Output directory '{}' for event '{}' does "
                               "not exist.".format(outdir,
                                                   self[self._KEYS.NAME]))
        save_name = os.path.join(outdir, filename + '.json')
        with codecs.open(save_name, 'w', encoding='utf8') as sf:
            sf.write(jsonstring)

        return save_name

    def _clean_quantity(self, quantity):
        """Clean quantity value before it is added to entry.
        """
        pass

    def _load_data_from_json(self, fhand, clean=False):
        """FIX: check for overwrite??
        """
        self._log.debug("_load_data_from_json(): {}\n\t{}".format(
            self.name(), fhand))
        # Store the filename this was loaded from
        self.filename = fhand
        with open(fhand, 'r') as jfil:
            data = json.load(jfil, object_pairs_hook=OrderedDict)
            name = list(data.keys())
            if len(name) != 1:
                raise ValueError("json file '{}' has multiple keys: {}".format(
                    fhand, list(name)))
            name = name[0]
            # Remove the outmost dict level
            data = data[name]
            self._log.debug("Name: {}".format(name))

            # Convert the OrderedDict data from json into class structure i.e.
            # `Sources` will be extracted and created from the dict Everything
            # that remains afterwards should be okay to just store to this
            # `Entry`
            self._convert_odict_to_classes(data, clean=clean)
            if len(data):
                err_str = ("Remaining entries in `data` after "
                           "`_convert_odict_to_classes`.")
                err_str += "\n{}".format(dict_to_pretty_string(data))
                self._log.error(err_str)
                raise RuntimeError(err_str)

        # If object doesnt have a name yet, but json does, store it
        self_name = self[KEYS.NAME]
        if len(self_name) == 0:
            self[KEYS.NAME] = name
        # Warn if there is a name mismatch
        elif self_name.lower().strip() != name.lower().strip():
            self._log.warning("Object name '{}' does not match name in json:"
                              "'{}'".format(self_name, name))

        self.check()
        return

    def _convert_odict_to_classes(self, data, clean=False):
        """
        """
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
            data.pop(schema_key)

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
            # data['sources'] = newsources
            # self.setdefault(src_key, []).extend(newsources)

        # Handle `photometry`
        # -------------------
        photo_key = self._KEYS.PHOTOMETRY
        if photo_key in data:
            photoms = data.pop(photo_key)
            self._log.debug("Found {} '{}' entries".format(
                len(photoms), photo_key))
            new_photoms = []
            for photo in photoms:
                new_photoms.append(Photometry(self, **photo))
            # data[photo_key] = new_photoms
            self.setdefault(photo_key, []).extend(new_photoms)

        # Handle `spectra`
        # ---------------
        spec_key = self._KEYS.SPECTRA
        if spec_key in data:
            # When we are cleaning internal data, we don't always want to
            # require all of the normal spectrum data elements.
            spectra = data.pop(spec_key)
            self._log.debug("Found {} '{}' entries".format(
                len(spectra), spec_key))
            new_specs = []
            for spec in spectra:
                new_specs.append(Spectrum(self, **spec))
            # data[spec_key] = new_specs
            self.setdefault(spec_key, []).extend(new_specs)

        # Handle `error`
        # --------------
        err_key = self._KEYS.ERRORS
        if err_key in data:
            errors = data.pop(err_key)
            self._log.debug("Found {} '{}' entries".format(
                len(errors), err_key))
            new_errors = []
            for err in errors:
                new_errors.append(Error(self, **err))
            # data[err_key] = new_errors
            self.setdefault(err_key, []).extend(new_errors)

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
                new_quantities = []
                for vv in vals:
                    new_quantities.append(Quantity(self, name=key, **vv))

                self.setdefault(key, []).extend(new_quantities)

        return

    def _check_cat_dict_source(self, cat_dict_class, key_in_self, **kwargs):
        # Make sure that a source is given
        source = kwargs.get(cat_dict_class._KEYS.SOURCE, None)
        if source is None:
            raise ValueError("{}: `source` must be provided!".format(
                self[self._KEYS.NAME]))
        # If this source/data is erroneous, skip it
        if self.is_erroneous(key_in_self, source):
            self._log.info("This source is erroneous, skipping")
            return None
        return source

    def _init_cat_dict(self, cat_dict_class, key_in_self, **kwargs):
        # Catch errors associated with crappy, but not unexpected data
        # log warning if instructed
        try:
            new_entry = cat_dict_class(self, key=key_in_self, **kwargs)
        except CatDictError as err:
            if err.warn:
                self._log.info("'{}' Not adding '{}': '{}'".format(
                    self[self._KEYS.NAME], key_in_self, str(err)))
            return None
        return new_entry

    def _add_cat_dict(self, cat_dict_class, key_in_self, **kwargs):
        # Make sure that a source is given, and is valid (nor erroneous)
        source = self._check_cat_dict_source(
            cat_dict_class, key_in_self, **kwargs)
        if source is None:
            return False

        # Try to create a new instance of this subclass of `CatDict`
        new_entry = self._init_cat_dict(cat_dict_class, key_in_self, **kwargs)
        if new_entry is None:
            return False

        # Compare this new entry with all previous entries to make sure is new
        for item in self.get(key_in_self, []):
            if new_entry.is_duplicate_of(item):
                item.append_sources_from(new_entry)
                # Return the entry in case we want to use any additional tags
                # to augment the old entry
                return new_entry

        self.setdefault(key_in_self, []).append(new_entry)
        return True
