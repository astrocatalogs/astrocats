"""Testnovae specific catalog class."""
import codecs
import json
import os
from collections import OrderedDict
from datetime import datetime
# from subprocess import check_output

from astrocats.catalog.catalog import Catalog
from astrocats.catalog.struct import QUANTITY
from astrocats.catalog.utils import read_json_arr, read_json_dict

from .testnova import TESTNOVA, Testnova
from .utils import name_clean


class TestnovaCatalog(Catalog):
    """Catalog class for `Testnova` objects."""

    class PATHS(Catalog.PATHS):
        """Paths to catalog inputs/outputs."""

        PATH_BASE = os.path.abspath(os.path.dirname(__file__))

        def __init__(self, catalog):
            """Initialize paths."""
            super(TestnovaCatalog.PATHS, self).__init__(catalog)
            # auxiliary datafiles
            self.TYPE_SYNONYMS = os.path.join(self.PATH_INPUT, 'type-synonyms.json')
            self.SOURCE_SYNONYMS = os.path.join(self.PATH_INPUT, 'source-synonyms.json')
            self.URL_REDIRECTS = os.path.join(self.PATH_INPUT, 'url-redirects.json')
            self.NON_SNE_TYPES = os.path.join(self.PATH_INPUT, 'non-sne-types.json')
            self.NON_SNE_PREFIXES = os.path.join(self.PATH_INPUT, 'non-sne-prefixes.json')
            self.BIBERRORS = os.path.join(self.PATH_INPUT, 'biberrors.json')
            self.ATELS = os.path.join(self.PATH_INPUT, 'atels.json')
            self.CBETS = os.path.join(self.PATH_INPUT, 'cbets.json')
            self.IAUCS = os.path.join(self.PATH_INPUT, 'iaucs.json')
            # cached datafiles
            self.BIBAUTHORS = os.path.join(self.PATH_CACHE, 'bibauthors.json')
            self.EXTINCT = os.path.join(self.PATH_CACHE, 'extinctions.json')

        def get_all_repo_folders(self, *args, **kwargs):
            """
            """
            all_repos = super().get_all_repo_folders(*args, **kwargs)
            # repos = [rep for rep in all_repos if 'testcat-output' in rep]
            # self.catalog.log.warning("Skipping all non-targeted repos...")
            # return repos
            return all_repos

    def __init__(self, args, log):
        """Initialize catalog."""
        # Initialize super `astrocats.catalog.catalog.Catalog` object
        super(TestnovaCatalog, self).__init__(args, log, git_clone=True)
        self.proto = Testnova
        self._load_aux_data()
        return

    def should_bury(self, name):
        """Determine whether an entry should be "buried".

        An entry would be buried if it does not belong to the class of object
        associated with the given catalog.
        """
        (bury_entry, save_entry) = super(TestnovaCatalog, self).should_bury(name)

        ct_val = None
        if name.startswith(tuple(self.nonsneprefixes_dict)):
            self.log.debug(
                "Killing '{}', non-SNe prefix.".format(name))
            save_entry = False
        else:
            if TESTNOVA.CLAIMED_TYPE in self.entries[name]:
                for ct in self.entries[name][TESTNOVA.CLAIMED_TYPE]:
                    up_val = ct[QUANTITY.VALUE].upper().replace('?', '')
                    up_types = [x.upper() for x in self.nonsnetypes]
                    if up_val not in up_types and up_val != 'CANDIDATE':
                        bury_entry = False
                        save_entry = True
                        break
                    if up_val in up_types:
                        bury_entry = True
                        ct_val = ct[QUANTITY.VALUE]
            else:
                if (TESTNOVA.DISCOVER_DATE in self.entries[name] and
                    any([x.get(QUANTITY.VALUE).startswith('AT')
                         for x in self.entries[name][TESTNOVA.ALIAS]]) and
                    not any([x.get(QUANTITY.VALUE).startswith('SN')
                             for x in self.entries[name][TESTNOVA.ALIAS]])):
                    try:
                        try:
                            dd = datetime.strptime(self.entries[name][
                                TESTNOVA.DISCOVER_DATE][0].get('value', ''),
                                '%Y/%m/%d')
                        except ValueError:
                            dd = datetime.strptime(self.entries[name][
                                TESTNOVA.DISCOVER_DATE][0].get('value', '') +
                                '/12/31',
                                '%Y')
                    except ValueError:
                        pass
                    else:
                        diff = datetime.today() - dd
                        # Because of the TNS, many non-SNe beyond 2016.
                        if dd.year >= 2016 and diff.days > 180:
                            save_entry = False

            if not save_entry:
                self.log.warning(
                    "Not saving '{}', {}.".format(name, ct_val))
            elif bury_entry:
                self.log.info(
                    "Burying '{}', {}.".format(name, ct_val))

        return (bury_entry, save_entry)

    def _load_aux_data(self):
        """Load auxiliary dictionaries for use in this catalog."""
        # Create/Load auxiliary dictionaries
        self.nedd_dict = OrderedDict()
        self.bibauthor_dict = read_json_dict(self.PATHS.BIBAUTHORS)
        self.biberror_dict = read_json_dict(self.PATHS.BIBERRORS)
        self.extinctions_dict = read_json_dict(self.PATHS.EXTINCT)
        self.iaucs_dict = read_json_dict(self.PATHS.IAUCS)
        self.cbets_dict = read_json_dict(self.PATHS.CBETS)
        self.atels_dict = read_json_dict(self.PATHS.ATELS)
        self.source_syns = read_json_dict(self.PATHS.SOURCE_SYNONYMS)
        self.url_redirs = read_json_dict(self.PATHS.URL_REDIRECTS)
        self.type_syns = read_json_dict(self.PATHS.TYPE_SYNONYMS)
        # Create/Load auxiliary arrays
        self.nonsneprefixes_dict = read_json_arr(
            self.PATHS.NON_SNE_PREFIXES)
        self.nonsnetypes = read_json_arr(self.PATHS.NON_SNE_TYPES)
        return

    def save_caches(self):
        """Save caches to JSON files."""
        jsonstring = json.dumps(self.bibauthor_dict, indent='\t',
                                separators=(',', ':'), ensure_ascii=False)
        with codecs.open(self.PATHS.BIBAUTHORS, 'w', encoding='utf8') as f:
            f.write(jsonstring)
        jsonstring = json.dumps(self.extinctions_dict, indent='\t',
                                separators=(',', ':'), ensure_ascii=False)
        with codecs.open(self.PATHS.EXTINCT, 'w', encoding='utf8') as f:
            f.write(jsonstring)

    def clean_entry_name(self, name):
        """Clean entry's name."""
        return name_clean(name)
