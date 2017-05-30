"""
"""
# import os
import sys
from collections import OrderedDict
from copy import deepcopy

import numpy as np

from astrocats.catalog.utils import tq
from astrocats.catalog.source import SOURCE
from . import utils as production_utils
from . import producer, html_pro
from .. entry import ENTRY
from .. quantity import QUANTITY


class Director(producer.Producer_Base):

    # Subclasses should modify these variables to add custom keys to saving/deleting.
    #    Combined with defaults (below) into `self.SAVE_ENTRY_KEYS` and `self.DEL_QUANTITY_KEYS`
    _SAVE_ENTRY_KEYS = []
    _DEL_QUANTITY_KEYS = []

    DEF_SAVE_ENTRY_KEYS = [
        "name", "alias", "ra", "dec", "altitude", "azimuth",
        "instruments", "redshift", "lumdist",
        "references", "sources"
    ]

    DEF_DEL_QUANTITY_KEYS = [
        'source', 'u_value', 'e_value', 'e_upper_value',
        'e_lower_value', 'derived'
    ]

    def __init__(self, catalog, args):
        log = catalog.log
        log.debug("Director.__init__()")
        self.log = log
        self.catalog = catalog

        self.output_data = OrderedDict()
        self.names_data = OrderedDict()

        self.SAVE_ENTRY_KEYS = self.DEF_SAVE_ENTRY_KEYS + self._SAVE_ENTRY_KEYS
        self.DEL_QUANTITY_KEYS = self.DEF_DEL_QUANTITY_KEYS + self._DEL_QUANTITY_KEYS
        log.info("`SAVE_ENTRY_KEYS` = '{}'".format(", ".join(self.SAVE_ENTRY_KEYS)))
        log.info("`DEL_QUANTITY_KEYS` = '{}'".format(", ".join(self.DEL_QUANTITY_KEYS)))

        self.web_table_min_fname = self.catalog.PATHS.WEB_TABLE_MIN_FILE
        self.web_table_fname = self.catalog.PATHS.WEB_TABLE_FILE
        self.names_fname = self.catalog.PATHS.NAMES_FILE
        self.names_min_fname = self.catalog.PATHS.NAMES_MIN_FILE
        log.info("`web_table_min_fname` = '{}'".format(self.web_table_min_fname))
        log.info("`web_table_fname` = '{}'".format(self.web_table_fname))
        log.info("`names_fname` = '{}'".format(self.names_fname))

        return

    def direct(self, args, write_collected=True, write_individual=True):
        catalog = self.catalog
        log = catalog.log
        _str = ["", "", "="*100, "", ""]
        for ss in _str:
            log.debug(ss)
        log.debug("Director.direct()")

        log.warning("Running `direct` on catalog {} ('{}')".format(catalog.name, type(catalog)))

        event_filenames = catalog.PATHS.get_repo_output_file_list(
            normal=(not args.boneyard), bones=args.boneyard)
        event_filenames = sorted(event_filenames, key=lambda s: s.lower())

        # log.warning("Shuffling!")
        # np.random.shuffle(event_filenames)

        num_events = len(event_filenames)
        log.warning("{} Files, e.g. '{}'".format(num_events, np.random.choice(event_filenames)))

        producers = []
        # MD5 File Checksums
        md5_pro = producer.MD5_Pro(catalog, args)
        producers.append(md5_pro)
        # Bibliography creator
        bib_pro = producer.Bib_Pro(catalog, args)
        producers.append(bib_pro)

        # Collect host-images
        img_pro = producer.Host_Image_Pro(catalog, args)
        # producers.append(img_pro)

        # Initialize an HTML Producer (for web html tables)
        web_pro = html_pro.HTML_Pro(catalog, args)
        # NOTE: this needs to be handled specially; dont add to list
        # producers.append(web_pro)

        # Iterate over all events
        # -----------------------
        for event_count, event_fname in enumerate(tq(event_filenames)):

            if args.travis and (event_count >= catalog.TRAVIS_QUERY_LIMIT):
                self.log.warning("Reached travis limit ({})".format(event_count))
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

            '''
            # Generate checksum for each json input file, compare to previous (if they exist)
            #    to determine if a file needs to be reprocessed
            md5_pro.update(event_fname, entry, event_data)

            # Store all bibliography and authors information
            bib_pro.update(event_fname, entry, event_data)
            '''

            # Collect host-images
            host_image_html = img_pro.update(event_fname, entry, event_data)

            # Generate HTML file for this event
            web_pro.update(event_fname, entry, event_data, host_image_html=host_image_html)

            for pro in producers:
                pro.update(event_fname, entry, event_data)

            # Add this entry into the catalog after removing undesired elements
            self.update(event_fname, entry, event_data)
            if args.test:
                sys.exit(232)

        # Do not save additional files if only targeting select events
        if args.event_list is not None:
            return

        # Write it all out at the end
        # ---------------------------
        '''
        md5_pro.finish()
        bib_pro.finish()
        if args.hosts:
            img_pro.finish()
        '''
        for pro in producers:
            pro.finish()

        self.finish()

        return

    def update(self, fname, event_name, event_data):
        self.log.debug("Director.update()")

        # Store names and aliases to a special dictionary
        if ENTRY.ALIAS in event_data:
            store_names = [x[QUANTITY.VALUE] for x in event_data[ENTRY.ALIAS]]
        else:
            store_names = [event_data[ENTRY.NAME]]
        self.names_data[event_data[ENTRY.NAME]] = store_names

        # For events with internal json files, add {'download' = 'e'} elements into dicts
        self.log.debug("Checking for internal event '{}'".format(event_name))
        if self.catalog.PATHS.is_internal_event(event_name):
            event_data['download'] = 'e'
            self.log.debug("...exists")

        # Store the top 5 ranking references (sources), and add an ID to each source
        if ENTRY.SOURCES in event_data:
            self.log.debug("Ranking sources")
            ranked_sources = rank_sources(event_data)
            num_rs = 0 if ranked_sources is None else len(ranked_sources)
            self.log.debug("Retrieved {} ranked sources".format(num_rs))
            if ranked_sources is not None:
                event_data['references'] = ', '.join(ranked_sources[:5])

        # Store desired data (filtering keys, parameters, etc)
        self.output_data[event_name] = self._clean_event_dict(event_data)
        return

    def finish(self):
        self.log.debug("Director.finish()")
        self.log.debug("finishing web-table outputs")
        self._finish_web_table_output()
        self.log.debug("finishing names output")
        self._finish_names_output()
        return

    def _finish_web_table_output(self):
        # Convert to array since that's what datatables expects
        output = list(self.output_data.values())

        # Produce 'min'imal version of web-table
        self.log.info("Writing web-table (min)")
        self._save_json(self.web_table_min_fname, output, expanded=False)

        # Produce expanded version of web-table
        self.log.info("Writing web-table (full)")
        self._save_json(self.web_table_fname, output, expanded=True)
        return

    def _finish_names_output(self):
        self.log.info("Writing names table (min)")
        self._save_json(self.names_min_fname, self.names_data, expanded=False)

        self.log.info("Writing names table (full)")
        self._save_json(self.names_fname, self.names_data, expanded=True)
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


def rank_sources(event_data, skip_secondary=False):
    # Store sources for this entry to count their usage
    #    only store primary sources
    ranked_sources = {}
    for source in event_data[ENTRY.SOURCES]:
        if (SOURCE.SECONDARY in source) and (skip_secondary):
            continue

        alias = source[SOURCE.ALIAS]
        ranked_sources[alias] = {
            'count': 0
        }
        # Store some sort of identification for this source (preferrence: first to last)
        #    The `Source` specification requires one of these to exist
        _source_id = [SOURCE.BIBCODE, SOURCE.ARXIVID, SOURCE.URL, SOURCE.NAME]
        for sid in _source_id:
            if sid in source:
                ranked_sources[alias][sid] = source[sid]
                # Also store this to a universal key (regardless of what type)
                ranked_sources[alias]['ID'] = source[sid]
                break

    if len(ranked_sources) == 0:
        return None

    # Go through all data in this entry, and count the contributions from each source
    for key, vals in event_data.items():
        if isinstance(vals, list):
            for row in vals:
                if QUANTITY.SOURCE in row:
                    _split = row[QUANTITY.SOURCE].split(',')
                    for source in ranked_sources:
                        if source in _split:
                            if key == 'spectra':
                                ranked_sources[source]['count'] += 10
                            else:
                                ranked_sources[source]['count'] += 1

    ranked_sources = [src['ID'] for src in
                      sorted(ranked_sources.values(), key=lambda x: x['count'], reverse=True)]

    return ranked_sources
