"""
"""
import codecs
import json
import os
import warnings
from collections import OrderedDict

from .utils import get_event_filename


class KEYS:
    ALIAS = 'alias'
    BIBCODE = 'bibcode'
    COMOVING_DIST = 'comovingdist'
    DEC = 'dec'
    DISCOVER_DATE = 'discoverdate'
    DISCOVERER = 'discoverer'
    DISTINCT_FROM = 'distinctfrom'
    EBV = 'ebv'
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
    SPECTRUM = 'spectrum'
    URL = 'url'
    VELOCITY = 'velocity'


class Entry(OrderedDict):

    # Whether or not this entry is a 'stub'.  Assume False
    _stub = False
    filename = None
    catalog = None

    _KEYS = KEYS

    def __init__(self, catalog, name, stub=False):
        """Create a new `Entry` object with the given `name`.
        """
        # if not name:
        #     raise ValueError("New `Entry` objects must have a valid name!")
        self[KEYS.NAME] = name
        self._stub = stub
        self.filename = None
        self.catalog = catalog
        return

    @classmethod
    def init_from_file(cls, catalog, name=None, path=None, clean=False):
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
        new_entry._load_data_from_json(load_path)

        if clean:
            new_entry.clean_internal()

        return new_entry

    def __repr__(self):
        jsonstring = json.dumps({self[KEYS.NAME]: self},
                                indent='\t', separators=(',', ':'),
                                ensure_ascii=False)
        return jsonstring

    def _load_data_from_json(self, fhand):
        """FIX: check for overwrite??
        """
        with open(fhand, 'r') as jfil:
            data = json.load(jfil, object_pairs_hook=OrderedDict)
            name = list(data.keys())
            if len(name) != 1:
                raise ValueError("json file '{}' has multiple keys: {}".format(
                    fhand, list(name)))
            name = name[0]
            data = data[name]
            self.update(data)
        self.filename = fhand
        # If object doesnt have a name yet, but json does, store it
        self_name = self[KEYS.NAME]
        if len(self_name) == 0:
            self[KEYS.NAME] = name
        # Warn if there is a name mismatch
        elif self_name.lower().strip() != name.lower().strip():
            warnings.warn(("Object name '{}' does not match name in json:"
                           "'{}'").format(
                self_name, name))

        self.check()
        return

    def save(self, empty=False, bury=False, gz=False, final=False):
        outdir, filename = self._get_save_path(bury=bury)

        if final:
            self.sanitize()

        # FIX: use 'dump' not 'dumps'
        jsonstring = json.dumps({self[KEYS.NAME]: self},
                                indent='\t', separators=(',', ':'),
                                ensure_ascii=False)
        if not os.path.isdir(outdir):
            raise RuntimeError("Output directory '{}' for event '{}' does "
                               "not exist.".format(outdir, self[KEYS.NAME]))
        save_name = os.path.join(outdir, filename + '.json')
        with codecs.open(save_name, 'w', encoding='utf8') as sf:
            sf.write(jsonstring)

        return save_name

    def sanitize(self):
        return

    def get_aliases(self, includename=True):
        """Retrieve the aliases of this object as a list of strings.
        """
        # empty list if doesnt exist
        alias_quanta = self.get(KEYS.ALIAS, [])
        aliases = [aq['value'] for aq in alias_quanta]
        if includename and self[KEYS.NAME] not in aliases:
            aliases = [self[KEYS.NAME]] + aliases
        return aliases

    def get_stub(self):
        """Get a new `Entry` which contains the 'stub' of this one.

        The 'stub' is *only* the name and aliases.

        Usage:
        -----
        To convert a normal entry into a stub (for example), overwrite the
        entry in place, i.e.
        >>> entries[name] = entries[name].get_stub()

        """
        stub = type(self)(self.catalog, self[KEYS.NAME], stub=True)
        if KEYS.ALIAS in self.keys():
            stub[KEYS.ALIAS] = self[KEYS.ALIAS]
        return stub

    def clean_internal(self):
        """Clean input from 'internal', human added data.

        This is used in the 'Entry.init_from_file' method.
        """
        # Rebuild the sources if they exist
        bibcodes = []
        try:
            old_sources = self[KEYS.SOURCES]
            del self[KEYS.SOURCES]
            for ss, source in enumerate(old_sources):
                if KEYS.BIBCODE in source:
                    bibcodes.append(source[KEYS.BIBCODE])
                    self.add_source(bibcode=source[KEYS.BIBCODE])
                else:
                    self.add_source(
                        srcname=source[KEYS.NAME], url=source[KEYS.URL])
        except KeyError:
            pass

        return

    def get_event_text(eventfile):
        import gzip
        if eventfile.split('.')[-1] == 'gz':
            with gzip.open(eventfile, 'rt') as f:
                filetext = f.read()
        else:
            with open(eventfile, 'r') as f:
                filetext = f.read()
        return filetext
