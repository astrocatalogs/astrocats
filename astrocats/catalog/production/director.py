"""
"""
import os
from collections import OrderedDict
from copy import deepcopy

import numpy as np

from astrocats.catalog.utils import tq  # , dict_to_pretty_string
from . import utils as production_utils
from . import producer

catalog_fname_prefix = 'catalog'
names_fname_prefix = 'names'


class Director(producer.Producer_Base):

    SAVE_ENTRY_KEYS = [
        "check", "name", "alias", "discoverdate", "maxdate", "maxappmag",
        "maxabsmag", "host", "ra", "dec", "hostra", "hostdec", "hostoffsetang",
        "hostoffsetdist", "altitude", "azimuth", "airmass", "skybrightness",
        "instruments", "redshift", "velocity", "lumdist",
        "claimedtype", "ebv", "photolink", "spectralink", "radiolink", "xraylink",
        "references", "download", "responsive"
    ]

    DEL_QUANTITY_KEYS = ['source', 'u_value', 'e_value', 'e_upper_value',
                         'e_lower_value', 'derived']

    def __init__(self, catalog, args):
        log = catalog.log
        log.debug("Director.__init__()")
        self.log = log
        self.catalog = catalog

        self.output_data = OrderedDict()
        self.names_data = OrderedDict()

        self.web_table_min_fname = os.path.join(self.catalog.PATHS.PATH_OUTPUT, 'catalog.min.json')
        self.web_table_fname = os.path.join(self.catalog.PATHS.PATH_OUTPUT, 'catalog.json')
        self.names_fname = os.path.join(self.catalog.PATHS.PATH_OUTPUT, 'names.min.json')

        return

    def produce(self, args):
        catalog = self.catalog
        log = catalog.log
        log.debug("Director.produce()")

        log.warning("Running `produce` on catalog {} ('{}')".format(catalog.name, str(catalog)))

        event_filenames = catalog.PATHS.get_repo_output_file_list(
            normal=(not args.boneyard), bones=args.boneyard)
        event_filenames = sorted(event_filenames, key=lambda s: s.lower())

        # log.warning("Shuffling!")
        # np.random.shuffle(event_filenames)

        num_events = len(event_filenames)
        log.warning("{} Files, e.g. '{}'".format(num_events, np.random.choice(event_filenames)))

        md5_pro = producer.MD5_Pro(catalog, args)
        bib_pro = producer.Bib_Pro(catalog, args)

        # Iterate over all events
        # -----------------------
        for event_count, event_fname in enumerate(tq(event_filenames)):

            if args.travis and (event_count >= catalog.TRAVIS_QUERY_LIMIT):
                break

            event_name = production_utils.get_event_name_from_filename(event_fname)

            # Skip events that weren't specifically targetted
            if (args.event_list is not None) and (event_name not in args.event_list):
                log.debug("event_name '{}' not in event_list".format(event_name))
                continue
            else:
                log.info("Found targeted event '{}'".format(event_name))

            # Load this entry from the file
            #    `entry` is the event-name as recorded in the file, usually same as `event_name`.
            #    May differ if `entry` contains (e.g.) a slash, etc
            entry, event_data = production_utils.load_event_from_filename(event_fname, log)
            log.debug("entry = '{}' (from fname: '{}')".format(entry, event_name))

            # Generate checksum for each json input file, compare to previous (if they exist)
            #    to determine if a file needs to be reprocessed
            md5_pro.update(event_fname, entry, event_data)

            # Store all bibliography and authors information
            bib_pro.update(event_fname, entry, event_data)

            # Add this entry into the catalog after removing undesired elements
            self.update(event_fname, entry, event_data)

        # Do not save additional files if only targeting select events
        if args.event_list is not None:
            return

        # Write it all out at the end
        # ---------------------------
        md5_pro.finish()
        bib_pro.finish()
        self.finish()

        return

    def update(self, fname, event_name, event_data):
        self.log.debug("Director.update()")

        # Store names and aliases to a special dictionary
        if 'alias' in event_data:
            store_names = [x['value'] for x in event_data['alias']]
        else:
            store_names = [event_data['name']]
        self.names_data[event_data['name']] = store_names

        # For events with internal json files, add {'download' = 'e'} elements into dicts
        internal_event_fname = self.catalog.PATHS.get_filename_for_internal_event(event_name)
        self.log.debug("Checking for internal event '{}'".format(internal_event_fname))
        if os.path.isfile(internal_event_fname):
            event_data['download'] = 'e'
            self.log.debug("...exists")

        # Store the top 5 ranking references (sources)
        if 'sources' in event_data:
            ranked_sources = rank_sources(event_data)
            if ranked_sources is not None:
                event_data['references'] = ', '.join(ranked_sources[:5])

        # Store desired data
        self.output_data[event_name] = self._clean_event_dict(event_data)
        return

    def finish(self):
        self.log.debug("Director.finish()")

        self._finish_web_table_output()
        self._finish_names_output()
        return

    def _finish_web_table_output(self):
        # Convert to array since that's what datatables expects
        output = list(self.output_data.values())

        # Produce 'min'imal version of web-table
        self._save_json(self.web_table_min_fname, output, expanded=False)

        # Produce expanded version of web-table
        self._save_json(self.web_table_fname, output, expanded=True)
        return

    def _finish_names_output(self):
        self._save_json(self.names_fname, self.names_data, expanded=False)
        return

    def _clean_event_dict(self, event_data):
        new_data = OrderedDict()
        # Go through all entries in this event
        for key, vals in event_data.items():
            if key not in self.SAVE_ENTRY_KEYS:
                continue
            new_data[key] = deepcopy(vals)
            if vals and isinstance(vals, list):
                # For each quantity, remove the undesired elements
                for row in new_data[key]:
                    for tag in self.DEL_QUANTITY_KEYS:
                        if tag in row:
                            del row[tag]

        return new_data


def rank_sources(event_data):
    # Store sources for this entry to count their usage
    #    only store primary sources which contain a 'bibcode' and 'name'
    ranked_sources = {}
    for source in event_data['sources']:
        if ('name' not in source) or ('bibcode' not in source) or ('secondary' in source):
            continue

        alias = source['alias']
        ranked_sources[alias] = {
            'bibcode': source['bibcode'],
            'count': 0
        }

    if len(ranked_sources) == 0:
        return None

    # Go through all data in this entry, and count the contributions from each source
    for key, vals in event_data.items():
        if isinstance(vals, list):
            for row in vals:
                if 'source' in row:
                    _split = row['source'].split(',')
                    for source in ranked_sources:
                        if source in _split:
                            if key == 'spectra':
                                ranked_sources[source]['count'] += 10
                            else:
                                ranked_sources[source]['count'] += 1

    ranked_sources = [src['bibcode'] for src in
                      sorted(ranked_sources.values(), key=lambda x: x['count'], reverse=True)]

    return ranked_sources
