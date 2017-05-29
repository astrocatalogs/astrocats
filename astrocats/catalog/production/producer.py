"""
"""
import os
import json
import gzip
import re
from collections import OrderedDict

from astrocats.catalog import utils
from . import utils as production_utils
from .. source import SOURCE
from .. entry import ENTRY
from .. quantity import QUANTITY


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

    def _save(self, fname, data):
        """Write to text file.
        """
        self.log.debug("Producer_Base._save()")
        self.log.warning("Writing to '{}'".format(fname))
        with open(fname, 'w') as ff:
            ff.write(data)

        fsize = os.path.getsize(fname)/1024/1024
        self.log.info("size: '{:.2f}' MB".format(fsize))
        return fname

    def _save_json(self, fname, data, expanded=False, **dump_kwargs):
        """Write data to JSON file.
        """
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
        return fname

    def _save_gzip(self, fname, data, append_suffix=True):
        """Write to a gzipped file.

        NOTE: `data` must be gzip writable, i.e. bytes.  To convert from `str`, use `str.encode()`.
        """
        self.log.debug("Producer_Base._save_gzip()")
        if fname.endswith('.gz'):
            fname = fname[:-3]
        fname_gz = fname + '.gz'
        self.log.warning("Writing to gzip '{}'".format(fname_gz))
        with gzip.open(fname_gz, 'w') as ff:
            ff.write(data)

        fsize = os.path.getsize(fname)/1024/1024
        self.log.info("size: '{:.2f}' MB".format(fsize))
        return fname

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

    def touch(self, fname, times=None):
        self.log.debug("Producer_Base.touch(): touching '{}'".format(fname))
        with open(fname, 'a'):
            os.utime(fname, times)
        return


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

        # ads packagee only needed for loading author data
        if self.load_all_authors_flag:
            log.info("Initializing ADS...")
            self.ads = utils.import_ads()

        return

    def update(self, fname, event_name, event_data):
        self.log.debug("Bib_Pro.update()")
        bib_data = self.bib_data
        bib_authors_data = self.bib_authors_data
        bib_all_authors_data = self.bib_all_authors_data

        if ENTRY.SOURCES not in event_data:
            return

        for source in event_data[ENTRY.SOURCES]:
            if SOURCE.BIBCODE not in source:
                continue

            bc = source[SOURCE.BIBCODE]
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
                     (SOURCE.BIBCODE, bc), ('events', []), ('eventdates', []),
                     ('types', []), ('photocount', 0), ('spectracount', 0),
                     ('metacount', 0)])

            bib_data[bc]['events'].append(event_data['name'])

            for key in list(event_data.keys()):
                bcalias = source[SOURCE.ALIAS]
                lc = 0
                if key in ['name', 'sources', 'schema', 'photometry',
                           'spectra', 'errors']:
                    continue
                for quantum in event_data[key]:
                    if QUANTITY.SOURCE not in quantum:
                        break
                    if bcalias in quantum[QUANTITY.SOURCE].split(','):
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


class HTML_Pro(Producer_Base):

    COLUMN_KEY = [
        "check", "name", "alias", "discoverdate", "maxdate", "maxappmag",
        "maxabsmag", "host", "ra", "dec", "hostra", "hostdec", "hostoffsetang",
        "hostoffsetdist", "altitude", "azimuth", "airmass", "skybrightness", "instruments",
        "redshift", "velocity", "lumdist",
        "claimedtype", "ebv", "photolink", "spectralink", "radiolink", "xraylink",
        "references", "download", "responsive"
    ]

    EVENT_IGNORE_KEY = ["download"]

    _EVENT_PAGE_HEADER = [
        "", "Name", "Aliases", "Discovery Date", "Maximum Date [band]",
        r"<em>m</em><sub>max</sub> [band]", r"<em>M</em><sub>max</sub> [band]",
        "Host Name", "R.A.", "Dec.", "Host R.A.", "Host Dec.", "Host Offset (\")",
        "Host Offset (kpc)", "Alt. (°)", "Azi. (°)", "Airmass", "V<sub>sky</sub>",
        "Instruments/Bands", r"<em>z</em>",
        r"<em>v</em><sub>&#9737;</sub> (km/s)", r"<em>d</em><sub>L</sub> (Mpc)",
        "Claimed Type", "E(B-V)", "Photometry", "Spectra", "Radio", "X-ray",
        "References", "Download", ""
    ]

    _HEADER = [
        "", "Name", "Aliases", "Disc. Date", "Max Date",
        r"<em>m</em><sub>max</sub>", r"<em>M</em><sub>max</sub>", "Host Name",
        "R.A.", "Dec.", "Host R.A.", "Host Dec.", "Host Offset (\")",
        "Host Offset (kpc)", "Alt. (°)", "Azi. (°)", "Airmass", "V<sub>sky</sub>", "Instruments/Bands", r"<em>z</em>",
        r"<em>v</em><sub>&#9737;</sub> (km/s)", r"<em>d</em><sub>L</sub> (Mpc)",
        "Type", "E(B-V)", "Phot.", "Spec.", "Radio", "X-ray", "References", "Data",
        ""
    ]

    _TITLES = [
        "", "Name (IAU name preferred)", "Aliases",
        "Discovey Date (year-month-day)", "Date of Maximum (year-month-day)",
        "Maximum apparent AB magnitude", "Maximum absolute AB magnitude",
        "Host Name", "{moduletitle} J2000 Right Ascension (h:m:s)",
        "{moduletitle} J2000 Declination (d:m:s)",
        "Host J2000 Right Ascension (h:m:s)", "Host J2000 Declination (d:m:s)",
        "Host Offset (Arcseconds)", "Host Offset (kpc)",
        "Altitude (Degrees)", "Azimuth (Degrees)", "Airmass", "Sky Brightness in V (Mags per arcsecond^2)",
        "List of Instruments and Bands", "Redshift",
        "Heliocentric velocity (km/s)", "Luminosity distance (Mpc)",
        "Claimed Type", "Milky Way Reddening", "Photometry", "pectra", "Radio",
        "X-rays", "Bibcodes of references with most data on event",
        "Download and edit data", ""
    ]

    def __init__(self, catalog, args):
        log = catalog.log
        self.log = log
        log.debug("HTML_Pro.__init__()")

        self.module_dir = 'blackholes'
        self.module_url = 'holes.space'
        self.module_name = 'bh'
        self.module_title = 'Blackhole'

        self.HEADER = OrderedDict(list(zip(self.COLUMN_KEY, self._HEADER)))
        self.EVENT_PAGE_HEADER = OrderedDict(list(zip(self.COLUMN_KEY, self._EVENT_PAGE_HEADER)))
        _titles = [tt.format(moduletitle=self.module_title) for tt in self._TITLES]
        self.TITLES = OrderedDict(list(zip(self.COLUMN_KEY, _titles)))

        # outdir = "astrocats/" + self.module_dir + "/output/"
        # htmldir = "html/"
        self.HTML_OUT_DIR = catalog.PATHS.PATH_HTML

        html = '<html><title></title><body></body></html>'

        html = html.replace(
            '<body>',
            '''<body class='event-body'><div style="padding-bottom:8px;"><strong>Disclaimer:</strong> All data collected by the OSC was originally generated by others, if you intend to use this data in a publication, we ask that you please cite the linked sources and/or contact the sources of the data directly. Data sources are revealed by hovering over the data with your cursor.</div>'''
        )

        self.html = html

        return

    def produce(self, fname, event_name, event_data):
        self.log.debug("HTML_Pro.produce()")
        module_url = self.module_url
        module_dir = self.module_dir
        module_name = self.module_name
        html = self.html

        html = re.sub(
            r'(\<\/title\>)', r'''\1\n
            <base target="_parent" />\n
            <link rel="stylesheet" href="https://''' + module_url +
            '''/wp-content/themes/astrocats-child-theme/event-iframe.css" type="text/css">\n
            <script type="text/javascript" src="https://''' + module_url +
            '''/wp-content/plugins/transient-table/transient-table.js" type="text/js"></script>\n
            <script type="text/javascript">\n
                if(top==self)\n
                this.location="''' + event_name + '''"\n
            </script>''', html)

        # repfolder = get_rep_folder(catalog[entry], repofolders)
        html = re.sub(
            r'(\<\/body\>)', '<div class="event-download">' + r'<a href="' +
            r'../json/' + event_name + r'.json" download>' +
            r'Download all data for ' + event_name + r'</a></div>\n\1', html)
        issueargs = '?title=' + ('[' + event_name + '] <Descriptive issue title>').encode('ascii', 'xmlcharrefreplace').decode("utf-8") + '&body=' + \
            ('Please describe the issue with ' + event_name + '\'s data here, be as descriptive as possible! ' +
             'If you believe the issue appears in other events as well, please identify which other events the issue possibly extends to.').encode('ascii', 'xmlcharrefreplace').decode("utf-8")
        html = re.sub(r'(\<\/body\>)', '<div class="event-issue">' +
                      r'<a href="https://github.com/astrocatalogs/' + module_dir
                      + '/issues/new' + issueargs + r'" target="_blank">' +
                      r'Report an issue with ' + event_name + r'</a></div>\n\1',
                      html)

        newhtml = r'<div class="event-tab-div"><h3 class="event-tab-title">Event metadata</h3><table class="event-table"><tr><th width=100px class="event-cell">Quantity</th><th class="event-cell">Value<sup>Sources</sup> [Kind]</th></tr>\n'
        # Check if this corresponds to an 'internal' file
        edit = self.catalog.PATHS.is_internal_event(event_name)

        for key in self.COLUMN_KEY:
            if key in event_data and key not in self.EVENT_IGNORE_KEY and len(event_data[key]) > 0:
                keyhtml = ''
                if isinstance(event_data[key], str):
                    subentry = re.sub('<[^<]+?>', '', event_data[key])
                    keyhtml = keyhtml + subentry
                else:
                    for r, row in enumerate(event_data[key]):
                        if 'value' in row and 'source' in row:
                            sources = [
                                str(x)
                                for x in sorted(
                                    [x.strip() for x in row['source'].split(',')],
                                    key=lambda x: float(x) if utils.is_number(x) else float("inf")
                                )
                            ]
                            sourcehtml = ''
                            for s, source in enumerate(sources):
                                sourcehtml = sourcehtml + \
                                    (', ' if s > 0 else '') + r'<a href="#source' + \
                                    source + r'" target="_self">' + source + r'</a>'
                            keyhtml = keyhtml + (r'<br>' if r > 0 else '')
                            keyhtml = keyhtml + "<div class='stt'>"
                            if 'derived' in row and row['derived']:
                                keyhtml = keyhtml + '<span class="derived">'
                            keyhtml = keyhtml + row['value']
                            if ((key == 'maxdate' or key == 'maxabsmag' or
                                 key == 'maxappmag') and
                                    'maxband' in event_data and
                                    event_data['maxband']):
                                keyhtml = keyhtml + r' [' + event_data['maxband'][0]['value'] + ']'
                            if 'e_value' in row:
                                keyhtml = keyhtml + r' ± ' + row['e_value']
                            if 'derived' in row and row['derived']:
                                keyhtml = keyhtml + '</span>'

                            # Mark erroneous button
                            sourceids = []
                            idtypes = []
                            for alias in row['source'].split(','):
                                for source in event_data['sources']:
                                    if source['alias'] == alias:
                                        if 'bibcode' in source:
                                            sourceids.append(source['bibcode'])
                                            idtypes.append('bibcode')
                                        elif 'arxivid' in source:
                                            sourceids.append(source['arxivid'])
                                            idtypes.append('arxivid')
                                        else:
                                            sourceids.append(source['name'])
                                            idtypes.append('name')
                            if not sourceids or not idtypes:
                                raise (ValueError(
                                    'Unable to find associated source by alias!'
                                ))
                            keyhtml = (
                                keyhtml +
                                "<span class='sttt'><button class='sme' type='button' onclick='markError(\""
                                + event_name + "\", \"" + key + "\", \"" +
                                ','.join(idtypes) + "\", \"" +
                                ','.join(sourceids) + "\", \"" + edit + "\", \"" + module_name +
                                "\")'>Flag as erroneous</button></span>")
                            keyhtml = keyhtml + r'</div><sup>' + sourcehtml + r'</sup>'
                        elif isinstance(row, str):
                            keyhtml = keyhtml + (r'<br>' if r > 0 else '') + row.strip()

                if keyhtml:
                    newhtml = newhtml + r'<tr><td class="event-cell">'
                    if key not in [
                            'photolink', 'spectralink', 'radiolink',
                            'xraylink', 'name'
                    ]:
                        newhtml = (
                            newhtml + '<div class="stt">' +
                            self.EVENT_PAGE_HEADER[key] +
                            "<span class='sttright'><button class='saq' type='button' onclick='addQuantity(\""
                            + event_name + "\", \"" + key + "\", \"" + edit + "\", \"" + module_name +
                            "\")'>Add new value</button></span></div>")
                    else:
                        newhtml = newhtml + self.EVENT_PAGE_HEADER[key]
                    newhtml = newhtml + r'</td><td width=250px class="event-cell">' + keyhtml

                newhtml = newhtml + r'</td></tr>\n'
        newhtml = newhtml + r'</table><p><em>Values that are colored <span class="derived">purple</span> were computed by the OSC using values provided by the specified sources.</em></p>\n'
        newhtml = newhtml + r'<p><em>*Absolute magnitudes take into account luminosity distance and redshift decrements but not SED shape, thus the K-corrections used to determine absolute magnitudes are approximate.</em></p></div>\n\1'
        html = re.sub(r'(\<\/body\>)', newhtml, html)

        if 'sources' in event_data and len(event_data['sources']):
            newhtml = r'<div class="event-tab-div"><h3 class="event-tab-title">Sources of data</h3><table class="event-table"><tr><th width=30px class="event-cell">ID</th><th class="event-cell">Source Info</th></tr><tr><th colspan="2" class="event-cell">Primary Sources</th></tr>\n'
            first_secondary = False
            for source in event_data['sources']:
                biburl = ''
                if 'bibcode' in source:
                    biburl = 'http://adsabs.harvard.edu/abs/' + \
                        source['bibcode']

                refurl = ''
                if 'url' in source:
                    refurl = source['url']

                sourcename = source['name'] if 'name' in source else source[
                    'bibcode'] if 'bibcode' in source else source['arxivid']
                if not first_secondary and source.get('secondary', False):
                    first_secondary = True
                    newhtml += r'<tr><th colspan="2" class="event-cell">Secondary Sources</th></tr>\n'

                sourcelines = []
                if ('bibcode' not in source or
                        sourcename != source['bibcode']):
                    sourcelines.append(
                        ((r'<a href="' + refurl + '">') if refurl else '') +
                        sourcename.encode('ascii', 'xmlcharrefreplace').decode(
                            "utf-8") + ((r'</a>\n') if refurl else ''))
                if 'reference' in source:
                    sourcelines.append(source['reference'])
                if 'bibcode' in source:
                    sourcelines.append(r'\n[' + (
                        ('<a href="' + biburl + '">') if 'reference' in source
                        else '') + source['bibcode'] + (
                            r'</a>' if 'reference' in source else '') + ']')
                if source.get('private', False):
                    sourcelines.append(
                        r'<span class="private-source">Data from source not yet public</span>\n'
                    )

                newhtml = (newhtml + r'<tr><td class="event-cell" id="source' +
                           source['alias'] + '">' + source['alias'] +
                           r'</td><td width=250px class="event-cell">' +
                           r'<br>\n'.join(sourcelines) + r'</td></tr>\n')
            newhtml = newhtml + r'</table><em>Sources are presented in order of importation, not in order of importance.</em></div>\n'

        newhtml = newhtml + r'\n\1'

        html = re.sub(r'(\<\/body\>)', newhtml, html)

        fname_out = os.path.join(self.HTML_OUT_DIR, event_name + ".html")
        # self.touch(fname_out)
        self._save(fname_out, html)
        self._save_gzip(fname_out, html.encode())

        return
