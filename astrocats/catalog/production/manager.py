import json
import os
from collections import OrderedDict
from copy import deepcopy

import numpy as np

from astrocats.catalog.utils import tq
from . import utils as production_utils

catalog_fname_prefix = 'catalog'
names_fname_prefix = 'names'

columnkey = [
    "check", "name", "alias", "discoverdate", "maxdate", "maxappmag",
    "maxabsmag", "host", "ra", "dec", "hostra", "hostdec", "hostoffsetang",
    "hostoffsetdist", "altitude", "azimuth", "airmass", "skybrightness",
    "instruments", "redshift", "velocity", "lumdist",
    "claimedtype", "ebv", "photolink", "spectralink", "radiolink", "xraylink",
    "references", "download", "responsive"
]

prune_tags = ['source', 'u_value', 'e_value', 'e_upper_value',
              'e_lower_value', 'derived']


class Manager:

    def __init__(self, catalog):
        log = catalog.log
        log.debug("Manager.__init__()")
        self.log = log
        self.catalog = catalog

        return

    def produce(self, args):
        catalog = self.catalog
        log = catalog.log
        log.debug("Manager.produce()")

        log.warning("Running `produce` on catalog {} ('{}')".format(catalog.name, str(catalog)))

        output = OrderedDict()
        output_copy = OrderedDict()

        event_filenames = catalog.PATHS.get_repo_output_file_list(
            normal=(not args.boneyard), bones=args.boneyard)
        event_filenames = sorted(event_filenames, key=lambda s: s.lower())

        # log.warning("Shuffling!")
        # np.random.shuffle(event_filenames)

        num_events = len(event_filenames)
        log.warning("{} Files, e.g. '{}'".format(num_events, np.random.choice(event_filenames)))

        # Load existing MD5 checksums
        md5_fname = catalog.PATHS.get_md5_filename()
        log.info("Checking MD5 file '{}'".format(md5_fname))
        try:
            with open(md5_fname) as f:
                filetext = f.read()
            md5dict = json.loads(filetext)
            log.info("Loaded {} MD5 checksums".format(len(md5dict)))
        except FileNotFoundError:
            log.info("MD5 file does not exist.")
            md5dict = {}

        # Iterate over all events
        # -----------------------
        for event_count, event_fname in enumerate(tq(event_filenames)):

            event_name = production_utils.get_event_name_from_filename(event_fname)

            # Skip events that weren't specifically targetted
            if (args.event_list is not None) and (event_name not in args.event_list):
                log.debug("event_name {} not in event_list".format(event_name))
                continue

            if args.travis and (event_count >= catalog.TRAVIS_QUERY_LIMIT):
                break

            # Generate checksum for each json input file, compare to previous (if they exist)
            #    to determine if a file needs to be reprocessed
            checksum = production_utils.md5file(event_fname)
            if (event_fname not in md5dict) or (md5dict[event_fname] != checksum):
                md5dict[event_fname] = checksum

            # Add this entry into the catalog
            filetext = production_utils.get_event_text(event_fname)
            output.update(json.loads(filetext, object_pairs_hook=OrderedDict))
            # `entry` is the event-name as recorded in the file usually same as `event_name`
            #    may differ if `entry` contains (e.g.) a slash, etc
            entry = next(reversed(output))

            log.debug("entry = '{}' (from fname: '{}')".format(entry, event_name))

            internal_event_fname = catalog.PATHS.get_filename_for_internal_event(event_name)
            log.debug("Checking for internal event '{}'".format(internal_event_fname))
            if os.path.isfile(internal_event_fname):
                output[entry]['download'] = 'e'
                log.debug("...exists")

            # Save this stuff because next line will delete it.
            if 'sources' in output[entry]:
                ranked_sources = rank_sources(output[entry])
                if ranked_sources is not None:
                    output[entry]['references'] = ', '.join(ranked_sources[:5])

            # Delete unneeded data from output, add blank entries when data missing.
            output_copy[entry] = OrderedDict()
            for col in columnkey:
                if col in output[entry]:
                    output_copy[entry][col] = deepcopy(output[entry][col])

            del output[entry]

        if args.event_list is not None:
            return

        # Write it all out at the end
        # ---------------------------

        output = deepcopy(output_copy)

        # Write the MD5 checksums
        md5_fname = catalog.PATHS.get_md5_filename()
        log.warning("Writing checksums to '{}'".format(md5_fname))
        jsonstring = json.dumps(md5dict, indent='\t', separators=(',', ':'))
        with open(md5_fname, 'w') as ff:
            ff.write(jsonstring)

        # Prune extraneous fields not required for main catalog file
        log.info("Pruning output")
        output = prune_output(output, prune_tags)

        # Convert to array since that's what datatables expects
        output = list(output.values())

        cat_fname = os.path.join(catalog.PATHS.PATH_OUTPUT, catalog_fname_prefix + '.min.json')
        log.warning("Writing 'minimum' catalog to '{}'".format(cat_fname))
        json_str = json.dumps(output, separators=(',', ':'))
        with open(cat_fname, 'w') as ff:
            ff.write(json_str)

        cat_fname = os.path.join(catalog.PATHS.PATH_OUTPUT, catalog_fname_prefix + '.json')
        log.warning("Writing catalog to '{}'".format(cat_fname))
        json_str = json.dumps(output, indent='\t', separators=(',', ':'))
        with open(cat_fname, 'w') as ff:
            ff.write(json_str)

        names = OrderedDict()
        for ev in output:
            names[ev['name']] = [
                x['value'] for x in ev.get('alias', [{'value': ev['name']}])
            ]
        names_fname = os.path.join(catalog.PATHS.PATH_OUTPUT, names_fname_prefix + '.min.json')
        log.warning("Writing names to '{}'".format(names_fname))
        json_str = json.dumps(names, separators=(',', ':'))
        with open(names_fname, 'w') as ff:
            ff.write(json_str)

        return


def prune_output(output, prune_tags):
    output_copy = OrderedDict()
    for entry in output:
        output_copy[entry] = OrderedDict()
        for col in output[entry]:
            output_copy[entry][col] = deepcopy(output[entry][col])
            if output_copy[entry][col]:
                for row in output_copy[entry][col]:
                    for tag in prune_tags:
                        if tag in row:
                            del row[tag]

    return deepcopy(output_copy)


def rank_sources(entry_data):
    # Store sources for this entry to count their usage
    #    only store primary sources which contain a 'bibcode' and 'name'
    ranked_sources = {}
    for source in entry_data['sources']:
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
    for key, vals in entry_data.items():
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
