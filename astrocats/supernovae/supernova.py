"""
"""
import json
import warnings
from collections import OrderedDict

from astropy.time import Time as astrotime

from astrocats.catalog.entry import KEYS, Entry
from astrocats.catalog.utils import (alias_priority, bandmetaf, bandrepf,
                                     get_event_filename, get_sig_digits,
                                     is_number, jd_to_mjd,
                                     make_date_string, pretty_num,
                                     read_json_dict, tprint, trim_str_arr,
                                     uniq_cdl, get_sig_digits)
from astrocats.supernovae.utils import (frame_priority, host_clean, name_clean,
                                        radec_clean, same_tag_num,
                                        same_tag_str)
from cdecimal import Decimal

from .. import SCHEMA
from .constants import FILENAME, MAX_BANDS, PREF_KINDS, REPR_BETTER_QUANTITY


class SN_KEYS(KEYS):
    CLAIMED_TYPE = 'claimedtype'
    DISCOVERY_DATE = 'discoverdate'
    ERRORS = 'errors'


class Supernova(Entry):
    """
    NOTE: OrderedDict data is just the `name` values from the JSON file.
          I.e. it does not include the highest nesting level
          { name: DATA }, it *just* includes DATA

    FIX: does this need to be `ordered`???
    FIX: check that no stored values are empty/invalid (delete key in that
         case?)
    FIX: distinguish between '.filename' and 'get_filename'

    sources
    -   All sources must have SN_KEYS.NAME and 'alias' parameters
    -   FIX: is url required if no bibcode???
    -   FIX: consider changing 'alias' for each source to 'src_num' or
             something
    -   FIX: Make source aliases integers (instead of strings of integers)??
    -   FIX: have list of allowed 'source' parameters??
    -   FIX: create class for 'errors'
    -   FIX: class or list of valid quantities and units

    """

    filename = ''
    _source_syns = {}

    def __init__(self, catalog, name, stub=False):
        super().__init__(catalog, name, stub=stub)

        # FIX: move this somewhere else (shouldnt be in each event)
        # Load source-name synonyms
        with open(FILENAME.SOURCE_SYNONYMS, 'r') as f:
            self._source_syns = json.loads(
                f.read(), object_pairs_hook=OrderedDict)
        return

    def add_source(self, srcname='', bibcode='', **src_kwargs):
        """Add a new source to this entry's SN_KEYS.SOURCES list.

        FIX: if source already exists, should dictionary be updated to any
             new values??

        Arguments
        ---------

        Returns
        -------
        src_alias : str (of integer)
            The alias number for this source.

        Notes
        -----
        Suggested `src_kwargs`:
            'url', 'secondary', 'acknowledgment', 'reference'

        """
        # Try to figure out each `srcname` or `bibcode` from the other, when
        # only one given
        if not srcname or not bibcode:
            srcname, bibcode = self._parse_srcname_bibcode(srcname, bibcode)

        self.catalog.log.debug("`srcname`: '{}', `bibcode`: '{}'".format(
            srcname, bibcode))

        # These are empty lists if no sources
        my_sources = self.get(SN_KEYS.SOURCES, [])
        my_src_aliases = [src[SN_KEYS.ALIAS] for src in my_sources]
        nsources = len(my_sources)

        # Try to find existing, matching source
        # -------------------------------------
        # If this source name already exists, return alias number
        try:
            my_src_names = [src[SN_KEYS.NAME] for src in my_sources]
            name_idx = my_src_names.index(srcname)
            return my_src_aliases[name_idx]
        # `KeyError` from `SN_KEYS.NAME` not existing, `ValueError` from
        # `srcname` not existing
        except (KeyError, ValueError):
            pass

        # If this bibcode already exists, return alias number
        try:
            my_src_bibs = [src[SN_KEYS.BIBCODE] for src in my_sources]
            bib_idx = my_src_bibs.index(bibcode)
            return my_src_aliases[bib_idx]
        # `KeyError` from `SN_KEYS.BIBCODE` not existing, `ValueError` from
        # `bibcode` not existing
        except (KeyError, ValueError):
            pass

        # Add new source that doesnt exist
        # --------------------------------
        source_alias = str(nsources + 1)
        new_src = OrderedDict()
        new_src[SN_KEYS.NAME] = srcname
        if bibcode:
            new_src[SN_KEYS.BIBCODE] = bibcode
        new_src[SN_KEYS.ALIAS] = source_alias
        # Add in any additional arguments passed (e.g. url, acknowledgment,
        # etc)
        new_src.update({k: v for (k, v) in src_kwargs.items() if k})
        self.setdefault(SN_KEYS.SOURCES, []).append(new_src)
        return source_alias

    def add_quantity(self, quantity, value, sources,
                     forcereplacebetter=False, derived='',
                     lowerlimit='', upperlimit='', error='', unit='',
                     kind='', extra='', probability=''):
        """
        """
        if not quantity:
            raise ValueError(self[SN_KEYS.NAME] +
                             "'s quantity must be specified for "
                             "add_quantity.")
        if not sources:
            raise ValueError(self[SN_KEYS.NAME] + "'s source must be specified for "
                             "quantity " +
                             quantity + ' before it is added.')
        if ((not isinstance(value, str) and
             (not isinstance(value, list) or not isinstance(value[0], str)))):
            raise ValueError(self[SN_KEYS.NAME] + "'s Quantity " + quantity +
                             " must be a string or an array of strings.")

        if self.is_erroneous(quantity, sources):
            return None

        my_quantity_list = self.get(quantity, [])

        svalue = value.strip()
        serror = error.strip()
        skind = kind.strip()
        sprob = probability.strip()
        sunit = ''

        if not svalue or svalue == '--' or svalue == '-':
            return
        if serror and (not is_number(serror) or float(serror) < 0):
            raise ValueError(self[SN_KEYS.NAME] + "'s quanta " + quantity +
                             ' error value must be a number and positive.')

        # Set default units
        if not unit and quantity == 'velocity':
            unit = 'KM/s'
        if not unit and quantity == 'ra':
            unit = 'hours'
        if not unit and quantity == 'dec':
            unit = 'degrees'
        if not unit and quantity in ['lumdist', 'comovingdist']:
            unit = 'Mpc'

        # Handle certain quantity
        if quantity == 'alias':
            svalue = name_clean(svalue)
            for df in self.get(SN_KEYS.DISTINCTS, []):
                if svalue == df['value']:
                    return

        if quantity in ['velocity', 'redshift', 'ebv', 'lumdist',
                        'comovingdist']:
            if not is_number(svalue):
                return
        if quantity == 'host':
            if is_number(svalue):
                return
            if svalue.lower() in ['anonymous', 'anon.', 'anon',
                                  'intergalactic']:
                return
            svalue = host_clean(svalue)
            if ((not skind and ((svalue.lower().startswith('abell') and
                                 is_number(svalue[5:].strip())) or
                                'cluster' in svalue.lower()))):
                skind = 'cluster'
        elif quantity == SN_KEYS.CLAIMED_TYPE:
            isq = False
            svalue = svalue.replace('young', '')
            if svalue.lower() in ['unknown', 'unk', '?', '-']:
                return
            if '?' in svalue:
                isq = True
                svalue = svalue.strip(' ?')
            for rep in self._source_syns:
                if svalue in self._source_syns[rep]:
                    svalue = rep
                    break
            if isq:
                svalue = svalue + '?'

        elif quantity in ['ra', 'dec', 'hostra', 'hostdec']:
            (svalue, sunit) = radec_clean(svalue, quantity, unit=unit)
        elif quantity == 'maxdate' or quantity == 'discoverdate':
            # Make sure month and day have leading zeroes
            sparts = svalue.split('/')
            if len(sparts[0]) > 4 and int(sparts[0]) > 0:
                raise ValueError('Date years limited to four digits.')
            if len(sparts) >= 2:
                svalue = sparts[0] + '/' + sparts[1].zfill(2)
            if len(sparts) == 3:
                svalue = svalue + '/' + sparts[2].zfill(2)

            for ii, ct in enumerate(my_quantity_list):
                # Only add dates if they have more information
                if len(ct['value'].split('/')) > len(svalue.split('/')):
                    return

        if is_number(svalue):
            svalue = '%g' % Decimal(svalue)
        if serror:
            serror = '%g' % Decimal(serror)

        for ii, ct in enumerate(my_quantity_list):
            if ct['value'] == svalue and sources:
                if 'kind' in ct and skind and ct['kind'] != skind:
                    return
                for source in sources.split(','):
                    if source not in my_quantity_list[ii]['source'].split(','):
                        my_quantity_list[ii]['source'] += ',' + source
                        if serror and 'error' not in my_quantity_list[ii]:
                            my_quantity_list[ii]['error'] = serror
                        if sprob and 'probability' not in my_quantity_list[ii]:
                            my_quantity_list[ii]['probability'] = sprob
                return

        if not sunit:
            sunit = unit

        quanta_entry = OrderedDict()
        quanta_entry['value'] = svalue
        if serror:
            quanta_entry['error'] = serror
        if sources:
            quanta_entry['source'] = sources
        if skind:
            quanta_entry['kind'] = skind
        if sprob:
            quanta_entry['probability'] = sprob
        if sunit:
            quanta_entry['unit'] = sunit
        if lowerlimit:
            quanta_entry['lowerlimit'] = lowerlimit
        if upperlimit:
            quanta_entry['upperlimit'] = upperlimit
        if derived:
            quanta_entry['derived'] = derived
        if extra:
            quanta_entry['extra'] = extra
        if (forcereplacebetter or quantity in REPR_BETTER_QUANTITY) and \
                len(my_quantity_list):
            newquantities = []
            isworse = True
            if quantity in ['discoverdate', 'maxdate']:
                for ct in my_quantity_list:
                    ctsplit = ct['value'].split('/')
                    svsplit = svalue.split('/')
                    if len(ctsplit) < len(svsplit):
                        isworse = False
                        continue
                    elif len(ctsplit) < len(svsplit) and len(svsplit) == 3:
                        val_one = max(2, get_sig_digits(
                            ctsplit[-1].lstrip('0')))
                        val_two = max(2, get_sig_digits(
                            svsplit[-1].lstrip('0')))
                        if val_one < val_two:
                            isworse = False
                            continue
                    newquantities.append(ct)
            else:
                newsig = get_sig_digits(svalue)
                for ct in my_quantity_list:
                    if 'error' in ct:
                        if serror:
                            if float(serror) < float(ct['error']):
                                isworse = False
                                continue
                        newquantities.append(ct)
                    else:
                        if serror:
                            isworse = False
                            continue
                        oldsig = get_sig_digits(ct['value'])
                        if oldsig >= newsig:
                            newquantities.append(ct)
                        if newsig >= oldsig:
                            isworse = False
            if not isworse:
                newquantities.append(quanta_entry)
            self[quantity] = newquantities
        else:
            self.setdefault(quantity, []).append(quanta_entry)
        return

    def is_erroneous(self, field, sources):
        if hasattr(self, SN_KEYS.ERRORS):
            my_errors = self['errors']
            for alias in sources.split(','):
                source = self.get_source_by_alias(alias)
                bib_err_values = [err['value'] for err in my_errors
                                  if err['kind'] == 'bibcode' and
                                  err['extra'] == field]
                if 'bibcode' in source and source['bibcode'] in bib_err_values:
                    return True

                name_err_values = [err['value'] for err in my_errors
                                   if err['kind'] == 'name' and
                                   err['extra'] == field]
                if 'name' in source and source['name'] in name_err_values:
                    return True

        return False

    def check(self):
        # Make sure there is a schema key in dict
        if SN_KEYS.SCHEMA not in self.keys():
            self[SN_KEYS.SCHEMA] = SCHEMA.URL
        # Make sure there is a name key in dict
        if SN_KEYS.NAME not in self.keys() or len(self[SN_KEYS.NAME]) == 0:
            raise ValueError("Supernova name is empty:\n\t{}".format(
                json.dumps(self, indent=2)))
        return

    def _get_save_path(self, bury=False):
        filename = get_event_filename(self[SN_KEYS.NAME])

        # Put non-SNe in the boneyard
        if bury:
            outdir = self.catalog.get_repo_boneyard()

        # Get normal repository save directory
        else:
            repo_folders = self.catalog.get_output_repo_folders()
            if SN_KEYS.DISCOVERY_DATE in self.keys():
                repo_years = self.catalog.get_repo_years()
                for r, year in enumerate(repo_years):
                    if int(self[SN_KEYS.DISCOVERY_DATE][0]['value'].
                           split('/')[0]) <= year:
                        outdir = repo_folders[r]
                        break
            else:
                outdir = repo_folders[0]

        return outdir, filename

    def sanitize(self):
        # Calculate some columns based on imported data, sanitize some fields
        name = self['name']

        aliases = self.get_aliases(includename=False)
        if name not in aliases:
            if 'sources' in self:
                self.add_quantity('alias', name, '1')
            else:
                source = self.add_source(
                    bibcode=self.catalog.OSC_BIBCODE,
                    srcname=self.catalog.OSC_NAME, url=self.catalog.OSC_URL,
                    secondary=True)
                self.add_quantity('alias', name, source)

        if ((name.startswith('SN') and is_number(name[2:6]) and
             'discoverdate' in self and
             int(self['discoverdate'][0]['value'].
                 split('/')[0]) >= 2016 and
             not any(['AT' in x for x in aliases]))):
            source = self.add_source(
                bibcode=self.catalog.OSC_BIBCODE, srcname=self.catalog.OSC_NAME,
                url=self.catalog.OSC_URL, secondary=True)
            self.add_quantity('alias', 'AT' + name[2:], source)

        self['alias'] = list(
            sorted(self['alias'],
                   key=lambda key: alias_priority(name, key)))
        aliases = self.get_aliases()

        if 'claimedtype' in self:
            # FIX: this is something that should be done completely internally
            #      i.e. add it to `clean` or something??
            self['claimedtype'] = self.ct_list_prioritized()
        if 'claimedtype' in self:
            self['claimedtype'][:] = [ct for ct in self[
                'claimedtype'] if (ct['value'] != '?' and ct['value'] != '-')]
            if not len(self['claimedtype']):
                del(self['claimedtype'])
        if 'claimedtype' not in self and name.startswith('AT'):
            source = self.add_source(
                bibcode=self.catalog.OSC_BIBCODE, srcname=self.catalog.OSC_NAME,
                url=self.catalog.OSC_URL, secondary=True)
            self.add_quantity('claimedtype', 'Candidate', source)

        if 'photometry' in self:
            self['photometry'].sort(
                key=lambda x: ((float(x['time']) if isinstance(x['time'], str)
                                else min([float(y) for y in x['time']])) if
                               'time' in x else 0.0,
                               x['band'] if 'band' in x else '',
                               float(x['magnitude']) if
                               'magnitude' in x else ''))
        if ('spectra' in self and
                list(filter(None, ['time' in x
                                   for x in self['spectra']]))):
            self['spectra'].sort(key=lambda x: (
                float(x['time']) if 'time' in x else 0.0))
        if 'sources' in self:
            for source in self['sources']:
                if 'bibcode' in source:
                    import urllib
                    from html import unescape
                    # First sanitize the bibcode
                    if len(source['bibcode']) != 19:
                        source['bibcode'] = urllib.parse.unquote(
                            unescape(source['bibcode'])).replace('A.A.', 'A&A')
                    if source['bibcode'] in self.catalog.biberror_dict:
                        source['bibcode'] = \
                            self.catalog.biberror_dict[source['bibcode']]

                    if source['bibcode'] not in self.catalog.bibauthor_dict:
                        bibcode = source['bibcode']
                        adsquery = (self.catalog.ADS_BIB_URL +
                                    urllib.parse.quote(bibcode) +
                                    '&data_type=Custom&format=%253m%20%25(y)')
                        response = urllib.request.urlopen(adsquery)
                        html = response.read().decode('utf-8')
                        hsplit = html.split("\n")
                        if len(hsplit) > 5:
                            bibcodeauthor = hsplit[5]
                        else:
                            bibcodeauthor = ''

                        if not bibcodeauthor:
                            warnings.warn(
                                "Bibcode didn't return authors, not converting"
                                "this bibcode.")

                        self.catalog.bibauthor_dict[bibcode] = unescape(
                            bibcodeauthor).strip()

            for source in self['sources']:
                if ('bibcode' in source and
                        source['bibcode'] in self.catalog.bibauthor_dict and
                        self.catalog.bibauthor_dict[source['bibcode']]):
                    source['reference'] = self.catalog.bibauthor_dict[
                        source['bibcode']]
                    if 'name' not in source and source['bibcode']:
                        source['name'] = source['bibcode']
        if 'redshift' in self:
            self['redshift'] = list(
                sorted(self['redshift'], key=lambda key:
                       frame_priority(key)))
        if 'velocity' in self:
            self['velocity'] = list(
                sorted(self['velocity'], key=lambda key:
                       frame_priority(key)))
        if 'claimedtype' in self:
            self['claimedtype'] = self.ct_list_prioritized()

    '''
    def save(self, empty=False, bury=False, final=False, gz=False):
        outdir, filename = self._get_save_path(bury=bury)

        if final:
            self.sanitize()

        # FIX: use 'dump' not 'dumps'
        jsonstring = json.dumps({self[SN_KEYS.NAME]: self},
                                indent='\t', separators=(',', ':'),
                                ensure_ascii=False)
        if not os.path.isdir(outdir):
            raise RuntimeError("Output directory '{}' for event '{}' does "
                               "not exist.".format(outdir, self[SN_KEYS.NAME]))
        save_name = os.path.join(outdir, filename + '.json')
        with codecs.open(save_name, 'w', encoding='utf8') as sf:
            sf.write(jsonstring)

        return save_name
    '''

    def get_source_by_alias(self, alias):
        for source in self.get(SN_KEYS.SOURCES, []):
            if source['alias'] == alias:
                return source
        raise ValueError(
            "Source '{}': alias '{}' not found!".format(
                self[SN_KEYS.NAME], alias))

    def _parse_srcname_bibcode(self, srcname, bibcode):
        # If no `srcname` is given, use `bibcode` after checking its validity
        if not srcname:
            if not bibcode:
                raise ValueError(
                    "`bibcode` must be specified if `srcname` is not.")
            if len(bibcode) != 19:
                raise ValueError(
                    "Bibcode '{}' must be exactly 19 characters "
                    "long".format(bibcode))
            srcname = bibcode

        # If a `srcname` is given, try to set a `bibcode`
        elif not bibcode:
            if srcname.upper().startswith('ATEL'):
                srcname = srcname.replace(
                    'ATEL', 'ATel').replace('Atel', 'ATel')
                srcname = srcname.replace(
                    'ATel #', 'ATel ').replace('ATel#', 'ATel')
                srcname = srcname.replace('ATel', 'ATel ')
                srcname = ' '.join(srcname.split())
                atelnum = srcname.split()[-1]
                if is_number(atelnum) and atelnum in self.catalog.atels_dict:
                    bibcode = atels_dict[atelnum]

            if srcname.upper().startswith('CBET'):
                srcname = srcname.replace('CBET', 'CBET ')
                srcname = ' '.join(srcname.split())
                cbetnum = srcname.split()[-1]
                if is_number(cbetnum) and cbetnum in self.catalog.cbets_dict:
                    bibcode = cbets_dict[cbetnum]

            if srcname.upper().startswith('IAUC'):
                srcname = srcname.replace('IAUC', 'IAUC ')
                srcname = ' '.join(srcname.split())
                iaucnum = srcname.split()[-1]
                if is_number(iaucnum) and iaucnum in self.catalog.iaucs_dict:
                    bibcode = iaucs_dict[iaucnum]

        for rep in self._source_syns:
            if srcname in self._source_syns[rep]:
                srcname = rep
                break

        return srcname, bibcode

    def clean_internal(self):
        """Clean input data from the 'Supernovae/input/internal' repository.

        Extends `Entry.clean_internal`.
        FIX: instead of making changes in place to `dirty_event`, should a new
             event be created, values filled, then returned??
        FIX: currently will fail if no bibcode and no url
        """
        # Call `Entry.clean_internal()` first
        super().clean_internal()

        bibcodes = []
        try:
            for ss, source in enumerate(self[SN_KEYS.SOURCES]):
                if SN_KEYS.BIBCODE in source:
                    bibcodes.append(source[SN_KEYS.BIBCODE])
        except KeyError:
            pass

        # Clean some legacy fields
        if 'aliases' in self and isinstance(self['aliases'], list):
            source = self.add_source(
                bibcode=self.catalog.OSC_BIBCODE, srcname=self.catalog.OSC_NAME,
                url=self.catalog.OSC_URL, secondary=True)
            for alias in self['aliases']:
                self.add_quantity('alias', alias, source)
            del self['aliases']

        # FIX: should this be an error if false??
        if ((SN_KEYS.DISTINCTS in self and
             isinstance(self[SN_KEYS.DISTINCTS], list) and
             isinstance(self[SN_KEYS.DISTINCTS][0], str))):
            distinctfroms = [x for x in self[SN_KEYS.DISTINCTS]]
            del self[SN_KEYS.DISTINCTS]
            source = self.add_source(
                bibcode=self.catalog.OSC_BIBCODE, srcname=self.catalog.OSC_NAME,
                url=self.catalog.OSC_URL, secondary=True)
            for df in distinctfroms:
                self.add_quantity(SN_KEYS.DISTINCTS, df, source)

        if 'errors' in self and \
                isinstance(self['errors'], list) and \
                'sourcekind' in self['errors'][0]:
            source = self.add_source(
                bibcode=self.catalog.OSC_BIBCODE, srcname=self.catalog.OSC_NAME,
                url=self.catalog.OSC_URL, secondary=True)
            for err in self['errors']:
                self.add_quantity('error', err['quantity'], source,
                                  kind=err['sourcekind'], extra=err['id'])
            del self['errors']

        if not bibcodes:
            self.add_source(bibcode=self.catalog.OSC_BIBCODE,
                            srcname=self.catalog.OSC_NAME,
                            url=self.catalog.OSC_URL, secondary=True)
            bibcodes = [self.catalog.OSC_BIBCODE]

        # Go through all keys in 'dirty' event
        for key in self.keys():
            if key in [SN_KEYS.NAME, SN_KEYS.SCHEMA,
                       SN_KEYS.SOURCES, SN_KEYS.ERRORS]:
                pass
            elif key == 'photometry':
                for p, photo in enumerate(self['photometry']):
                    if photo['u_time'] == 'JD':
                        self['photometry'][p]['u_time'] = 'MJD'
                        self['photometry'][p]['time'] = str(
                            jd_to_mjd(Decimal(photo['time'])))
                    if bibcodes and 'source' not in photo:
                        source = self.add_source(bibcode=bibcodes[0])
                        self['photometry'][p]['source'] = source
            else:
                for qi, quantity in enumerate(self[key]):
                    if bibcodes and 'source' not in quantity:
                        source = self.add_source(bibcode=bibcodes[0])
                        self[key][qi]['source'] = source

        self.check()
        return

    def add_photometry(self, time="", u_time="MJD", e_time="",
                       telescope="", instrument="", band="", magnitude="",
                       e_magnitude="", source="", upperlimit=False, system="",
                       scorrected="", observatory="", observer="", host=False,
                       includeshost=False, survey="", kcorrected="", flux="",
                       fluxdensity="", e_flux="", e_fluxdensity="", u_flux="",
                       u_fluxdensity="", frequency="", u_frequency="",
                       counts="",
                       e_counts="", nhmw="", photonindex="", unabsorbedflux="",
                       e_unabsorbedflux="", energy="", u_energy="",
                       e_lower_magnitude="", e_upper_magnitude="",
                       e_lower_time="", e_upper_time="", mcorrected=""):
        name = self[SN_KEYS.NAME]
        if (not time and not host) or (not magnitude and not flux and not
                                       fluxdensity and not counts and not
                                       unabsorbedflux):
            warnings.warn(
                "Time or brightness not specified when adding photometry, not "
                "adding.")
            tprint('Name : "' + name + '", Time: "' + time + '", Band: "' +
                   band + '", AB magnitude: "' + magnitude + '"')
            return

        if ((not host and not is_number(time)) or
            (not is_number(magnitude) and not is_number(flux) and not
             is_number(fluxdensity) and not is_number(counts))):
            warnings.warn('Time or brightness not numerical, not adding.')
            tprint('Name : "' + name + '", Time: "' + time + '", Band: "' +
                   band + '", AB magnitude: "' + magnitude + '"')
            return

        if (((e_magnitude and not is_number(e_magnitude)) or
             (e_flux and not is_number(e_flux)) or
             (e_fluxdensity and not is_number(e_fluxdensity)) or
             (e_counts and not is_number(e_counts)))):
            warnings.warn('Brightness error not numerical, not adding.')
            tprint('Name : "' + name + '", Time: "' + time +
                   '", Band: "' + band + '", AB error: "' + e_magnitude + '"')
            return

        if e_time and not is_number(e_time):
            warnings.warn('Time error not numerical, not adding.')
            tprint('Name : "' + name + '", Time: "' +
                   time + '", Time error: "' + e_time + '"')
            return

        if ((flux or fluxdensity) and ((not u_flux and not u_fluxdensity) or
                                       (not frequency and not band and not
                                        energy))):
            warnings.warn(
                "Unit and band/frequency must be set when adding photometry by "
                "flux or flux density, not adding.")
            tprint('Name : "' + name + '", Time: "' + time)
            return

        if not source:
            ValueError('Photometry must have source before being added!')

        if self.is_erroneous('photometry', source):
            return

        # Do some basic homogenization
        sband = bandrepf(band)

        sinstrument = instrument
        ssystem = system
        stelescope = telescope

        if not sinstrument:
            sinstrument = bandmetaf(sband, 'instrument')
        if not stelescope:
            stelescope = bandmetaf(sband, 'telescope')
        if not ssystem:
            ssystem = bandmetaf(sband, 'system')

        # Look for duplicate data and don't add if duplicate
        if 'photometry' in self:
            for photo in self['photometry']:
                if ((same_tag_str(photo, sband, 'band') and
                     same_tag_str(photo, u_time, 'u_time') and
                     same_tag_num(photo, time, 'time', canbelist=True) and
                     same_tag_num(photo, magnitude, 'magnitude') and
                     (('host' not in photo and not host) or
                      ('host' in photo and host)) and
                     same_tag_num(photo, flux, 'flux') and
                     same_tag_num(photo, unabsorbedflux, 'unabsorbedflux') and
                     same_tag_num(photo, fluxdensity, 'fluxdensity') and
                     same_tag_num(photo, counts, 'counts') and
                     same_tag_num(photo, energy, 'energy', canbelist=True) and
                     same_tag_num(photo, frequency, 'frequency') and
                     same_tag_num(photo, photonindex, 'photonindex') and
                     same_tag_num(photo, e_magnitude, 'e_magnitude') and
                     same_tag_num(photo, e_lower_time, 'e_lower_time') and
                     same_tag_num(photo, e_upper_time, 'e_upper_time') and
                     same_tag_num(photo, e_lower_magnitude,
                                  'e_lower_magnitude') and
                     same_tag_num(photo, e_upper_magnitude,
                                  'e_upper_magnitude') and
                     same_tag_num(photo, e_flux, 'e_flux') and
                     same_tag_num(photo, e_unabsorbedflux, 'e_unabsorbedflux') and
                     same_tag_num(photo, e_fluxdensity, 'e_fluxdensity') and
                     same_tag_num(photo, e_counts, 'e_counts') and
                     same_tag_str(photo, u_flux, 'u_flux') and
                     same_tag_str(photo, u_fluxdensity, 'u_fluxdensity') and
                     same_tag_str(photo, u_frequency, 'u_frequency') and
                     same_tag_str(photo, u_energy, 'u_energy') and
                     same_tag_num(photo, u_flux, 'u_flux') and
                     same_tag_num(photo, u_fluxdensity, 'u_fluxdensity') and
                     same_tag_num(photo, u_frequency, 'u_frequency') and
                     same_tag_num(photo, u_energy, 'u_energy') and
                     same_tag_str(photo, ssystem, 'system'))):
                    return

        photoentry = OrderedDict()
        if time:
            photoentry['time'] = time if isinstance(
                time, list) or isinstance(time, str) else str(time)
        if e_time:
            photoentry['e_time'] = str(e_time)
        if e_lower_time:
            photoentry['e_lower_time'] = str(e_lower_time)
        if e_upper_time:
            photoentry['e_upper_time'] = str(e_upper_time)
        if u_time:
            photoentry['u_time'] = u_time
        if sband:
            photoentry['band'] = sband
        if ssystem:
            photoentry['system'] = ssystem
        if magnitude:
            photoentry['magnitude'] = str(magnitude)
        if e_magnitude:
            photoentry['e_magnitude'] = str(e_magnitude)
        if e_lower_magnitude:
            photoentry['e_lower_magnitude'] = str(e_lower_magnitude)
        if e_upper_magnitude:
            photoentry['e_upper_magnitude'] = str(e_upper_magnitude)
        if frequency:
            photoentry['frequency'] = frequency if isinstance(
                frequency, list) or isinstance(frequency, str) else str(frequency)
        if u_frequency:
            photoentry['u_frequency'] = u_frequency
        if energy:
            photoentry['energy'] = energy if isinstance(
                energy, list) or isinstance(energy, str) else str(energy)
        if u_energy:
            photoentry['u_energy'] = u_energy
        if flux:
            photoentry['flux'] = str(flux)
        if e_flux:
            photoentry['e_flux'] = str(e_flux)
        if unabsorbedflux:
            photoentry['unabsorbedflux'] = str(unabsorbedflux)
        if e_unabsorbedflux:
            photoentry['e_unabsorbedflux'] = str(e_unabsorbedflux)
        if u_flux:
            photoentry['u_flux'] = str(u_flux)
        if photonindex:
            photoentry['photonindex'] = str(photonindex)
        if fluxdensity:
            photoentry['fluxdensity'] = str(fluxdensity)
        if e_fluxdensity:
            photoentry['e_fluxdensity'] = str(e_fluxdensity)
        if u_fluxdensity:
            photoentry['u_fluxdensity'] = str(u_fluxdensity)
        if counts:
            photoentry['counts'] = str(counts)
        if e_counts:
            photoentry['e_counts'] = str(e_counts)
        if upperlimit:
            photoentry['upperlimit'] = upperlimit
        if host:
            photoentry['host'] = host
        if includeshost:
            photoentry['includeshost'] = includeshost
        if kcorrected:
            photoentry['kcorrected'] = kcorrected
        if scorrected:
            photoentry['scorrected'] = scorrected
        if mcorrected:
            photoentry['mcorrected'] = mcorrected
        if observer:
            photoentry['observer'] = observer
        if survey:
            photoentry['survey'] = survey
        if observatory:
            photoentry['observatory'] = observatory
        if stelescope:
            photoentry['telescope'] = stelescope
        if sinstrument:
            photoentry['instrument'] = sinstrument
        if nhmw:
            photoentry['nhmw'] = nhmw
        if source:
            photoentry['source'] = source
        self.setdefault('photometry', []).append(photoentry)
        return

    def add_spectrum(self, waveunit, fluxunit, wavelengths="", fluxes="",
                     u_time="", time="", instrument="", deredshifted="",
                     dereddened="", errorunit="", errors="", source="",
                     snr="", telescope="", observer="", survey="", reducer="",
                     filename="", observatory="", data=""):

        if self.is_erroneous('spectra', source):
            return

        spectrumentry = OrderedDict()

        if 'spectra' in self:
            for si, spectrum in enumerate(self['spectra']):
                if 'filename' in spectrum and spectrum['filename'] == filename:
                    # Copy exclude info
                    if 'exclude' in spectrum:
                        spectrumentry['exclude'] = spectrum['exclude']
                    # Don't add duplicate spectra
                    if 'data' in spectrum:
                        return
                    del self['spectra'][si]
                    break

        if not waveunit:
            warnings.warn('No error unit specified, not adding spectrum.')
            return
        if not fluxunit:
            warnings.warn('No flux unit specified, not adding spectrum.')
            return

        if not data or (not wavelengths or not fluxes):
            ValueError("Spectrum must have wavelengths and fluxes set, or data "
                       "set.")

        if not source:
            ValueError('Spectrum must have source before being added!')

        if deredshifted != '':
            spectrumentry['deredshifted'] = deredshifted
        if dereddened != '':
            spectrumentry['dereddened'] = dereddened
        if instrument:
            spectrumentry['instrument'] = instrument
        if telescope:
            spectrumentry['telescope'] = telescope
        if observatory:
            spectrumentry['observatory'] = observatory
        if u_time:
            spectrumentry['u_time'] = u_time
        if time:
            spectrumentry['time'] = time
        if snr:
            spectrumentry['snr'] = snr
        if observer:
            spectrumentry['observer'] = observer
        if reducer:
            spectrumentry['reducer'] = reducer
        if survey:
            spectrumentry['survey'] = survey
        if filename:
            spectrumentry['filename'] = filename

        spectrumentry['waveunit'] = waveunit
        spectrumentry['fluxunit'] = fluxunit
        if data:
            spectrumentry['data'] = data
        else:
            if errors and max([float(x) for x in errors]) > 0.:
                if not errorunit:
                    warnings.warn(
                        'No error unit specified, not adding spectrum.')
                    return
                spectrumentry['errorunit'] = errorunit
                data = [trim_str_arr(wavelengths), trim_str_arr(
                    fluxes), trim_str_arr(errors)]
            else:
                data = [trim_str_arr(wavelengths), trim_str_arr(fluxes)]
            spectrumentry['data'] = [list(i) for i in zip(*data)]
        if source:
            spectrumentry['source'] = source
        self.setdefault('spectra', []).append(spectrumentry)
        return

    def _get_max_light(self):
        if 'photometry' not in self:
            return (None, None, None, None)

        # FIX: THIS
        eventphoto = [(x['u_time'], x['time'],
                       Decimal(x['magnitude']), x[
            'band'] if 'band' in x else '',
                       x['source']) for x in self['photometry'] if
            ('magnitude' in x and 'time' in x and 'u_time' in x and
             'upperlimit' not in x)]
        if not eventphoto:
            return None, None, None, None

        mlmag = None
        for mb in MAX_BANDS:
            leventphoto = [x for x in eventphoto if x[3] in mb]
            if leventphoto:
                mlmag = min([x[2] for x in leventphoto])
                eventphoto = leventphoto
                break

        if not mlmag:
            mlmag = min([x[2] for x in eventphoto])

        mlindex = [x[2] for x in eventphoto].index(mlmag)
        mlband = eventphoto[mlindex][3]
        mlsource = eventphoto[mlindex][4]

        if eventphoto[mlindex][0] == 'MJD':
            mlmjd = float(eventphoto[mlindex][1])
            mlmjd = astrotime(mlmjd, format='mjd').datetime
            return mlmjd, mlmag, mlband, mlsource
        else:
            return None, mlmag, mlband, mlsource

    def _get_first_light(self):
        if 'photometry' not in self:
            return None, None

        # FIX THIS
        eventphoto = [(Decimal(x['time']) if isinstance(x['time'], str) else
                       Decimal(min(float(y) for y in x['time'])),
                       x['source']) for x in self['photometry'] if
                      'upperlimit' not in x and
                      'time' in x and 'u_time' in x and x['u_time'] == 'MJD']
        if not eventphoto:
            return None, None
        flmjd = min([x[0] for x in eventphoto])
        flindex = [x[0] for x in eventphoto].index(flmjd)
        flmjd = astrotime(float(flmjd), format='mjd').datetime
        flsource = eventphoto[flindex][1]
        return flmjd, flsource

    def set_first_max_light(self):
        if 'maxappmag' not in self:
            mldt, mlmag, mlband, mlsource = self._get_max_light()
            if mldt:
                source = self.add_source(
                    bibcode=self.catalog.OSC_BIBCODE,
                    srcname=self.catalog.OSC_NAME, url=self.catalog.OSC_URL,
                    secondary=True)
                max_date = make_date_string(mldt.year, mldt.month, mldt.day)
                self.add_quantity(
                    'maxdate', max_date,
                    uniq_cdl([source] + mlsource.split(',')),
                    derived=True)
            if mlmag:
                source = self.add_source(
                    bibcode=self.catalog.OSC_BIBCODE,
                    srcname=self.catalog.OSC_NAME, url=self.catalog.OSC_URL,
                    secondary=True)
                self.add_quantity(
                    'maxappmag', pretty_num(mlmag),
                    uniq_cdl([source] + mlsource.split(',')),
                    derived=True)
            if mlband:
                source = self.add_source(
                    bibcode=self.catalog.OSC_BIBCODE,
                    srcname=self.catalog.OSC_NAME, url=self.catalog.OSC_URL,
                    secondary=True)
                (self
                 .add_quantity('maxband',
                               mlband,
                               uniq_cdl([source] + mlsource.split(',')),
                               derived=True))

        if ('discoverdate' not in self or
                max([len(x['value'].split('/')) for x in
                     self['discoverdate']]) < 3):
            fldt, flsource = self._get_first_light()
            if fldt:
                source = self.add_source(
                    bibcode=self.catalog.OSC_BIBCODE,
                    srcname=self.catalog.OSC_NAME, url=self.catalog.OSC_URL,
                    secondary=True)
                disc_date = make_date_string(fldt.year, fldt.month, fldt.day)
                self.add_quantity(
                    'discoverdate', disc_date,
                    uniq_cdl([source] + flsource.split(',')),
                    derived=True)

        if 'discoverdate' not in self and 'spectra' in self:
            minspecmjd = float("+inf")
            for spectrum in self['spectra']:
                if 'time' in spectrum and 'u_time' in spectrum:
                    if spectrum['u_time'] == 'MJD':
                        mjd = float(spectrum['time'])
                    elif spectrum['u_time'] == 'JD':
                        mjd = float(jd_to_mjd(Decimal(spectrum['time'])))
                    else:
                        continue

                    if mjd < minspecmjd:
                        minspecmjd = mjd
                        minspecsource = spectrum['source']

            if minspecmjd < float("+inf"):
                fldt = astrotime(minspecmjd, format='mjd').datetime
                source = self.add_source(
                    bibcode=self.catalog.OSC_BIBCODE,
                    srcname=self.catalog.OSC_NAME, url=self.catalog.OSC_URL,
                    secondary=True)
                disc_date = make_date_string(fldt.year, fldt.month, fldt.day)
                self.add_quantity(
                    'discoverdate', disc_date,
                    uniq_cdl([source] + minspecsource.split(',')),
                    derived=True)
        return

    def get_best_redshift(self):
        bestsig = -1
        bestkind = 10
        for z in self['redshift']:
            kind = PREF_KINDS.index(z['kind'] if 'kind' in z else '')
            sig = get_sig_digits(z['value'])
            if sig > bestsig and kind <= bestkind:
                bestz = z['value']
                bestkind = kind
                bestsig = sig
                bestsrc = z['source']

        return bestz, bestkind, bestsig, bestsrc

    def ct_list_prioritized(self):
        ct_list = list(sorted(
            self[SN_KEYS.CLAIMED_TYPE], key=lambda key: self._ct_priority(key)))
        return ct_list

    def _ct_priority(self, attr):
        aliases = attr['source'].split(',')
        max_source_year = -10000
        vaguetypes = ['CC', 'I']
        if attr['value'] in vaguetypes:
            return -max_source_year
        for alias in aliases:
            if alias == 'D':
                continue
            source = self.get_source_by_alias(alias)
            if 'bibcode' in source:
                source_year = self.get_source_year(source)
                if source_year > max_source_year:
                    max_source_year = source_year
        return -max_source_year
