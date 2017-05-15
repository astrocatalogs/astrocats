import argparse
import json
import os
import re

from collections import OrderedDict
from copy import deepcopy
from decimal import Decimal

from astrocats.catalog.utils import (bandshortaliasf, bandwavef,
                                     pretty_num, production, tq)
from astrocats.scripts.repos import (get_all_rep_folders, repo_file_list)
from astropy.time import Time as astrotime

from .catalog.utils import logger
log = logger.get_logger(tofile='producer.log')
log.warning("astrocats.producer()")

parser = argparse.ArgumentParser(
    description='Generate a catalog JSON file and plot HTML files from AstroCats data.'
)
parser.add_argument(
    '--no-write-catalog',
    '-nwc',
    dest='writecatalog',
    help='Don\'t write catalog file',
    default=True,
    action='store_false')
parser.add_argument(
    '--no-write-html',
    '-nwh',
    dest='writehtml',
    help='Don\'t write html plot files',
    default=True,
    action='store_false')
parser.add_argument(
    '--no-collect-hosts',
    '-nch',
    dest='collecthosts',
    help='Don\'t collect host galaxy images',
    default=True,
    action='store_false')
parser.add_argument(
    '--force-html',
    '-fh',
    dest='forcehtml',
    help='Force write html plot files',
    default=False,
    action='store_true')
parser.add_argument(
    '--event-list',
    '-el',
    dest='eventlist',
    help='Process a list of events',
    default=[],
    type=str,
    nargs='+')
parser.add_argument(
    '--test',
    '-te',
    dest='test',
    help='Test this script',
    default=True,
    action='store_true')
parser.add_argument(
    '--travis',
    '-tr',
    dest='travis',
    help='Set some options when using Travis',
    default=False,
    action='store_true')
parser.add_argument(
    '--boneyard',
    '-by',
    dest='boneyard',
    help='Make "boneyard" catalog',
    default=False,
    action='store_true')
parser.add_argument(
    '--delete-orphans',
    '-do',
    dest='deleteorphans',
    help='Delete orphan JSON files',
    default=False,
    action='store_true')
parser.add_argument(
    '--catalog',
    '-c',
    dest='catalog',
    help='Select which catalog to generate',
    default='sne',
    type=str)
args = parser.parse_args()

log.warning("args = {}".format(args))

if args.catalog == 'tde':
    catalog_dir = 'tidaldisruptions'
    catalog_name_short = 'tde'
    catalog_url = 'tde.space'
    catalog_title = 'TDE'
elif args.catalog == 'sne':
    catalog_dir = 'supernovae'
    catalog_name_short = 'sne'
    catalog_url = 'sne.space'
    catalog_title = 'Supernova'
else:
    raise ValueError('Unknown catalog!')

repo_folders = get_all_rep_folders(catalog_dir)
log.warning("repo_folders = '{}'".format(repo_folders))

outdir = "astrocats/" + catalog_dir + "/output/"
cachedir = "cache/"
jsondir = "json/"
html_subdir = "html/"

travislimit = 100
radiosigma = 3.0

linkdir = "https://" + catalog_url + "/" + catalog_name_short + "/"

testsuffix = '.test' if args.test else ''

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

catalog = OrderedDict()
catalogcopy = OrderedDict()
sourcedict = {}

if os.path.isfile(outdir + cachedir + 'hostimgs.json'):
    with open(outdir + cachedir + 'hostimgs.json', 'r') as f:
        filetext = f.read()
    hostimgdict = json.loads(filetext)
else:
    hostimgdict = {}

files = repo_file_list(
    catalog_dir, repo_folders, normal=(not args.boneyard), bones=args.boneyard)
log.warning("{} Files, e.g. '{}'".format(len(files), files[0]))

# Load existing MD5 checksums
if os.path.isfile(outdir + cachedir + 'md5s.json'):
    with open(outdir + cachedir + 'md5s.json', 'r') as f:
        filetext = f.read()
    md5dict = json.loads(filetext)
else:
    md5dict = {}

# Iterate over all events
# -----------------------
for event_count, event_fname in enumerate(tq(sorted(files, key=lambda s: s.lower()))):
    event_name_from_fname = os.path.splitext(os.path.basename(event_fname))[0].replace(
        '.json', '')
    # Skip events that weren't specifically targetted
    if args.eventlist and event_name_from_fname not in args.eventlist:
        continue

    if args.travis and event_count >= travislimit:
        break

    # Generate checksum for each json input file, compare to previous (if they exist)
    # to determine if a file needs to be reprocessed
    entry_changed = False
    checksum = production.md5file(event_fname)
    if event_fname not in md5dict or md5dict[event_fname] != checksum:
        entry_changed = True
        md5dict[event_fname] = checksum

    # Add this entry into the catalog
    filetext = production.get_event_text(event_fname)
    catalog.update(json.loads(filetext, object_pairs_hook=OrderedDict))
    entry = next(reversed(catalog))

    event_name = entry
    # log.warning("event_name = '{}'".format(event_name))

    # Skip events that aren't targetted
    if args.eventlist and event_name not in args.eventlist:
        continue

    if os.path.isfile("astrocats/" + catalog_dir + "/input/" + catalog_name_short +
                      "-internal/" + event_name_from_fname + ".json"):
        catalog[entry]['download'] = 'e'
    if 'discoverdate' in catalog[entry]:
        for d, date in enumerate(catalog[entry]['discoverdate']):
            catalog[entry]['discoverdate'][d]['value'] = catalog[entry][
                'discoverdate'][d]['value'].split('.')[0]
    if 'maxdate' in catalog[entry]:
        for d, date in enumerate(catalog[entry]['maxdate']):
            catalog[entry]['maxdate'][d]['value'] = catalog[entry]['maxdate'][
                d]['value'].split('.')[0]

    # Get host (galaxy) magnitude information
    hostmag = ''
    hosterr = ''
    if 'photometry' in catalog[entry]:
        for photo in catalog[entry]['photometry']:
            if 'host' in photo and ('upperlimit' not in photo or not photo['upperlimit']):
                hostmag = float(photo['magnitude'])
                hosterr = float(photo['e_magnitude']) if 'e_magnitude' in photo else 0.0

        # Delete the host magnitudes so they are not plotted as points
        catalog[entry]['photometry'][:] = [
            x for x in catalog[entry]['photometry'] if 'host' not in x
        ]

    # Find what parameters exist in this entry
    photoavail = 'photometry' in catalog[entry] and any(
        ['magnitude' in x for x in catalog[entry]['photometry']])
    radioavail = 'photometry' in catalog[entry] and any([
        'fluxdensity' in x and 'magnitude' not in x
        for x in catalog[entry]['photometry']
    ])
    xrayavail = 'photometry' in catalog[entry] and any(
        [('countrate' in x or 'flux' in x) and 'magnitude' not in x
         for x in catalog[entry]['photometry']])
    spectraavail = 'spectra' in catalog[entry]

    realizchecks = ''

    # Must be two sigma above host magnitude, if host magnitude known, to add
    # to phot count.
    numphoto = len([
        x for x in catalog[entry]['photometry']
        if 'upperlimit' not in x and 'magnitude' in x and
        (not hostmag or 'includeshost' not in x or float(x['magnitude']) <= (
            hostmag - 2.0 * hosterr))
    ]) if photoavail else 0
    numradio = len([
        x for x in catalog[entry]['photometry']
        if 'upperlimit' not in x and 'fluxdensity' in x and 'magnitude' not in
        x and (not x['e_fluxdensity'] or float(x[
            'fluxdensity']) > radiosigma * float(x['e_fluxdensity'])
        ) and (not hostmag or 'includeshost' not in x or float(x[
            'magnitude']) <= (hostmag - 2.0 * hosterr))
    ]) if photoavail else 0
    numxray = len([
        x for x in catalog[entry]['photometry']
        if 'upperlimit' not in x and ('countrate' in x or 'flux' in x
                                      ) and 'magnitude' not in x and
        (not hostmag or 'includeshost' not in x or float(x['magnitude']) <= (
            hostmag - 2.0 * hosterr))
    ]) if photoavail else 0
    numspectra = len(catalog[entry]['spectra']) if spectraavail else 0

    # log.warning("numphoto = {}".format(numphoto))
    # log.warning("numradio = {}".format(numradio))
    # log.warning("numxray = {}".format(numxray))
    # log.warning("numspectra = {}".format(numspectra))

    redshiftfactor = (1.0 / (
        1.0 + float(catalog[entry]['redshift'][0]['value']))) if (
            'redshift' in catalog[entry]) else 1.0

    mjdmax = ''
    if 'maxdate' in catalog[entry]:
        datestr = catalog[entry]['maxdate'][0]['value']
        datesplit = datestr.split('/')
        if len(datesplit) < 2:
            datestr += "/01"
        if len(datesplit) < 3:
            datestr += "/01"
        try:
            mjdmax = astrotime(datestr.replace('/', '-')).mjd
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception:
            pass

    minphotoep = ''
    maxphotoep = ''
    if mjdmax:
        photoeps = [(Decimal(x['time']) - Decimal(mjdmax + 0.5)) *
                    Decimal(redshiftfactor)
                    for x in catalog[entry]['photometry']
                    if 'upperlimit' not in x and 'includeshost' not in x and
                    'magnitude' in x and 'time' in x] if photoavail else []
        if photoeps:
            minphotoep = pretty_num(float(min(photoeps)), sig=3)
            maxphotoep = pretty_num(float(max(photoeps)), sig=3)

    minspectraep = ''
    maxspectraep = ''
    if mjdmax:
        spectraeps = ([(Decimal(x['time']) - Decimal(mjdmax + 0.5)) *
                       Decimal(redshiftfactor)
                       for x in catalog[entry]['spectra'] if 'time' in x]
                      if spectraavail else [])
        if spectraeps:
            minspectraep = pretty_num(float(min(spectraeps)), sig=3)
            maxspectraep = pretty_num(float(max(spectraeps)), sig=3)

    catalog[entry]['numphoto'] = numphoto
    catalog[entry]['numradio'] = numradio
    catalog[entry]['numxray'] = numxray
    catalog[entry]['numspectra'] = numspectra

    distancemod = 0.0
    if 'maxabsmag' in catalog[entry] and 'maxappmag' in catalog[entry]:
        distancemod = float(production.get_first_value(catalog, entry, 'maxappmag')) - \
            float(production.get_first_value(catalog, entry, 'maxabsmag'))

    plotlink = catalog_name_short + "/" + event_name_from_fname + "/"
    if photoavail:
        catalog[entry]['photolink'] = (str(numphoto) + (
            (',' + minphotoep + ',' + maxphotoep) if
            (minphotoep and maxphotoep and minphotoep != maxphotoep
             ) else ((',' + minphotoep) if minphotoep and maxphotoep else '')))
    if radioavail:
        catalog[entry]['radiolink'] = str(numradio)
    if xrayavail:
        catalog[entry]['xraylink'] = str(numxray)
    if spectraavail:
        catalog[entry]['spectralink'] = (str(numspectra) + (
            (',' + minspectraep + ',' + maxspectraep)
            if (minspectraep and maxspectraep and minspectraep != maxspectraep
                ) else ((',' + minspectraep) if minspectraep and maxspectraep else '')))

    # range of photometry elements
    prange = list(range(len(catalog[entry][
        'photometry']))) if 'photometry' in catalog[entry] else []

    instrulist = sorted([
        _f
        for _f in list({
            catalog[entry]['photometry'][x]['instrument']
            if 'instrument' in catalog[entry]['photometry'][x] else None
            for x in prange
        }) if _f
    ])
    if len(instrulist) > 0:
        instruments = ''
        for i, instru in enumerate(instrulist):
            instruments += instru
            bandlist = sorted(
                [
                    _f
                    for _f in list({
                        bandshortaliasf(catalog[entry]['photometry'][x][
                            'band'] if 'band' in catalog[entry]['photometry'][
                                x] else '') if 'instrument' in catalog[entry]
                        ['photometry'][x] and catalog[entry]['photometry'][x][
                            'instrument'] == instru else ""
                        for x in prange
                    }) if _f
                ],
                key=lambda y: (bandwavef(y), y))
            if bandlist:
                instruments += ' (' + ", ".join(bandlist) + ')'
            if i < len(instrulist) - 1:
                instruments += ', '

        # Now add bands without attached instrument
        obandlist = sorted(
            [
                _f
                for _f in list({
                    bandshortaliasf(catalog[entry]['photometry'][x]['band']
                                    if 'band' in catalog[entry]['photometry'][
                                        x] else '') if 'instrument' not in
                    catalog[entry]['photometry'][x] else ""
                    for x in prange
                }) if _f
            ],
            key=lambda y: (bandwavef(y), y))
        if obandlist:
            instruments += ", " + ", ".join(obandlist)
        catalog[entry]['instruments'] = instruments
        # log.warning("'{}' added instruments: {}".format(entry, instruments))
    else:
        bandlist = sorted(
            [
                _f
                for _f in list({
                    bandshortaliasf(catalog[entry]['photometry'][x]['band']
                                    if 'band' in catalog[entry]['photometry'][
                                        x] else '')
                    for x in prange
                }) if _f
            ],
            key=lambda y: (bandwavef(y), y))
        if len(bandlist) > 0:
            catalog[entry]['instruments'] = ", ".join(bandlist)
            # log.warning("'{}' added instruments (bandlist): {}".format(
            #     entry, catalog[entry]['instruments']))

    # Save this stuff because next line will delete it.
    if args.writecatalog:
        if 'sources' in catalog[entry]:
            lsourcedict = {}
            for sourcerow in catalog[entry]['sources']:
                if 'name' not in sourcerow:
                    continue
                strippedname = re.sub('<[^<]+?>', '', sourcerow['name'].encode(
                    'ascii', 'xmlcharrefreplace').decode("utf-8"))
                alias = sourcerow['alias']
                if 'bibcode' in sourcerow and 'secondary' not in sourcerow:
                    lsourcedict[alias] = {
                        'bibcode': sourcerow['bibcode'],
                        'count': 0
                    }
                if strippedname in sourcedict:
                    sourcedict[strippedname] += 1
                else:
                    sourcedict[strippedname] = 1

            for key in catalog[entry].keys():
                if isinstance(catalog[entry][key], list):
                    for row in catalog[entry][key]:
                        if 'source' in row:
                            for lsource in lsourcedict:
                                if lsource in row['source'].split(','):
                                    if key == 'spectra':
                                        lsourcedict[lsource]['count'] += 10
                                    else:
                                        lsourcedict[lsource]['count'] += 1

            ssources = sorted(
                list(lsourcedict.values()), key=lambda x: x['count'], reverse=True)
            if ssources:
                catalog[entry]['references'] = ','.join([y['bibcode'] for y in ssources[:5]])
                # log.warning("'{}' added references: {} (from {})".format(
                #     entry, catalog[entry]['references'], _old_refs))

        # Delete unneeded data from catalog, add blank entries when data missing.
        catalogcopy[entry] = OrderedDict()
        for col in columnkey:
            if col in catalog[entry]:
                vals = catalog[entry][col]
                if isinstance(vals, list):
                    vals = sorted(vals)
                catalogcopy[entry][col] = deepcopy(catalog[entry][col])

    # print()
    # for k1, v1 in catalog[entry].items():
    #     out_str = "{}: {}".format(k1, v1)
    #     if k1 in catalogcopy[entry]:
    #         out_str += "\n\t{}: {}".format(k1, catalogcopy[entry][k1])
    #     print(out_str)

    del catalog[entry]
    # sys.exit(2342)

# Write it all out at the end
if args.writecatalog and not args.eventlist:
    catalog = deepcopy(catalogcopy)

    # Write the MD5 checksums
    jsonstring = json.dumps(md5dict, indent='\t', separators=(',', ':'))
    with open(outdir + cachedir + 'md5s.json' + testsuffix, 'w') as f:
        f.write(jsonstring)

    # Prune extraneous fields not required for main catalog file
    catalogcopy = OrderedDict()
    for entry in catalog:
        catalogcopy[entry] = OrderedDict()
        for col in catalog[entry]:
            catalogcopy[entry][col] = deepcopy(catalog[entry][col])
            if catalogcopy[entry][col]:
                for row in catalogcopy[entry][col]:
                    for tag in prune_tags:
                        if tag in row:
                            del row[tag]
    catalog = deepcopy(catalogcopy)

    # Convert to array since that's what datatables expects
    catalog = list(catalog.values())

    if args.boneyard:
        catprefix = 'bones'
    else:
        catprefix = 'catalog'

    jsonstring = json.dumps(catalog, separators=(',', ':'))
    with open(outdir + catprefix + '.min.json' + testsuffix, 'w') as f:
        f.write(jsonstring)

    jsonstring = json.dumps(catalog, indent='\t', separators=(',', ':'))
    with open(outdir + catprefix + '.json' + testsuffix, 'w') as f:
        f.write(jsonstring)

    names = OrderedDict()
    for ev in catalog:
        names[ev['name']] = [
            x['value'] for x in ev.get('alias', [{'value': ev['name']}])
        ]
    jsonstring = json.dumps(names, separators=(',', ':'))
    boneyard_fname = outdir + 'names' + ('-by' if args.boneyard else '') + '.min.json' + testsuffix
    with open(boneyard_fname, 'w') as f:
        f.write(jsonstring)
