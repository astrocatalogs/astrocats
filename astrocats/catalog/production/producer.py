"""
"""
import os
import json
from collections import OrderedDict

from . import utils as production_utils


class Producer_Base:

    def __init__(self, catalog, args):
        self.log = catalog.log
        self.log.debug("Producer_Base.__init__()")
        return

    def update(self, fname, event_name, event_data):
        self.log.debug("Producer_Base.update()")
        return

    def finish(self):
        self.log.debug("Producer_Base.finish()")
        return

    def _save_json(self, fname, data, expanded=False, **dump_kwargs):
        self.log.debug("Producer_Base._save_json()")
        dump_kwargs.setdefault('separators', (',', ':'))
        if expanded:
            dump_kwargs.setdefault('indent', '\t')

        self.log.debug("`dump_kwargs`: '{}'".format(dump_kwargs))

        jsonstring = json.dumps(data, **dump_kwargs)
        self.log.warning("Writing to json '{}'".format(fname))
        with open(fname, 'w') as ff:
            ff.write(jsonstring)

        fsize = os.path.getsize(fname)/1024/1024
        self.log.info("size: '{:.2f}' MB".format(fsize))
        return

    def _load_default_json(self, fname, default=OrderedDict()):
        """Attempt to load data from the given filename, if it doesn't exist, return the default.
        """
        self.log.debug("Producer_Base._load_default_json()")
        self.log.info("Trying to load data from '{}'".format(fname))

        try:
            with open(fname) as ff:
                file_text = ff.read()
            data = json.loads(file_text, object_pairs_hook=OrderedDict)
            self.log.debug("Loaded {} entries".format(len(data)))
        except FileNotFoundError:
            data = default
            self.log.debug("No file found, initializing default.")

        return data


class MD5_Pro(Producer_Base):
    """
    """

    def __init__(self, catalog, args):
        log = catalog.log
        self.log = log
        log.debug("MD5_Pro.__init__()")

        self.md5_fname = catalog.PATHS.MD5_FILE
        log.debug("`md5_fname` = '{}'".format(self.md5_fname))

        self.md5_data = self._load_default_json(self.md5_fname)
        self.count = 0
        return

    def update(self, fname, event_name, event_data):
        self.log.debug("MD5_Pro.update()")
        checksum = production_utils.load_md5_file(fname)
        changed = False
        if (fname not in self.md5_data) or (self.md5_data[fname] != checksum):
            self.md5_data[fname] = checksum
            self.count += 1
            changed = True
        return changed

    def finish(self):
        self.log.debug("MD5_Pro.finish()")
        self._save_json(self.md5_fname, self.md5_data)
        return


class Bib_Pro(Producer_Base):
    """
    """
    import ads

    def __init__(self, catalog, args):
        log = catalog.log
        self.log = log
        log.debug("Bib_Pro.__init__()")

        self.biblio_fname = catalog.PATHS.BIBLIO_FILE
        self.biblio_min_fname = catalog.PATHS.BIBLIO_MIN_FILE
        self.authors_fname = catalog.PATHS.AUTHORS_FILE
        self.all_authors_fname = catalog.PATHS.ALL_AUTHORS_FILE
        log.debug("`authors_fname` = '{}'".format(self.authors_fname))
        log.debug("`all_authors_fname` = '{}'".format(self.all_authors_fname))

        self.bib_data = OrderedDict()
        self.bib_authors_data = self._load_default_json(self.authors_fname)
        self.bib_all_authors_data = self._load_default_json(self.all_authors_fname)
        self.count = 0

        self.load_all_authors_flag = args.authors
        log.warning("`load_all_authors_flag` = '{}'".format(self.load_all_authors_flag))

        log.info("Initializing ADS...")
        self._init_ads()
        return

    def update(self, fname, event_name, event_data):
        self.log.debug("Bib_Pro.update()")
        bib_data = self.bib_data
        bib_authors_data = self.bib_authors_data
        bib_all_authors_data = self.bib_all_authors_data

        if 'sources' not in event_data:
            return

        for source in event_data['sources']:
            if 'bibcode' not in source:
                continue

            bc = source['bibcode']
            if bc not in bib_data:

                authors = ''
                if bc in bib_authors_data:
                    authors = bib_authors_data[bc]

                all_authors = []
                if bc in bib_all_authors_data and len(bib_all_authors_data[bc]):
                    all_authors = bib_all_authors_data[bc]
                elif self.load_all_authors_flag:
                    try:
                        q = list(self.ads.SearchQuery(bibcode=bc))
                        if not len(q):
                            q = list(self.ads.SearchQuery(alternate_bibcode=bc))
                        all_authors = q
                    except Exception as err:
                        self.log.debug("ADS search for '{}' failed: {}".foramt(bc, str(err)))

                    if all_authors and all_authors[0].author:
                        all_authors = all_authors[0].author

                    bib_all_authors_data[bc] = all_authors

                bib_data[bc] = OrderedDict(
                    [('authors', authors), ('all_authors', all_authors),
                     ('bibcode', bc), ('events', []), ('eventdates', []),
                     ('types', []), ('photocount', 0), ('spectracount', 0),
                     ('metacount', 0)])

            bib_data[bc]['events'].append(event_data['name'])

            for key in list(event_data.keys()):
                bcalias = source['alias']
                lc = 0
                if key in [
                        'name', 'sources', 'schema', 'photometry',
                        'spectra', 'errors'
                ]:
                    continue
                for quantum in event_data[key]:
                    if bcalias in quantum['source'].split(','):
                        lc += 1
                bib_data[bc]['metacount'] += lc

        return

    def finish(self):
        self.log.debug("Bib_Pro.finish()")
        # Each reference includes a list of events it has contributed, also include a count
        for key, vals in self.bib_data.items():
            num_events = len(vals['events'])
            vals['num_events'] = num_events

        self.log.warning("'{}' is not being constructed".format(self.authors_fname))
        # self._save_json(self.authors_fname, self.bib_authors_data)

        self._save_json(self.biblio_fname, self.bib_data, expanded=True)
        self._save_json(self.biblio_min_fname, self.bib_data, expanded=False)
        if self.load_all_authors_flag:
            self._save_json(self.all_authors_fname, self.bib_all_authors_data)
        return

    def _init_ads(self):
        path = 'ads.key'
        if os.path.isfile(path):
            with open(path, 'r') as ff:
                self.ads.config.token = ff.read().splitlines()[0]
        else:
            raise IOError(
                "Cannot find ads.key, please generate one at "
                "https://ui.adsabs.harvard.edu/#user/settings/token and place it in "
                "this file.")

        return
