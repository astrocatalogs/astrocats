"""
"""
import codecs
import json
import os
import warnings
from collections import OrderedDict

from cdecimal import Decimal
from scripts import FILENAME, PATH, SCHEMA

from ..utils import (bandmetaf, bandrepf,
                     get_repo_folders, get_repo_paths, get_repo_years,
                     get_sig_digits, is_number, tprint)
from .constants import (OSC_BIBCODE, OSC_NAME, OSC_URL, REPR_BETTER_QUANTITY)
from .funcs import (get_atels_dict, get_cbets_dict,
                    get_iaucs_dict, host_clean, jd_to_mjd, name_clean,
                    radec_clean, same_tag_num, same_tag_str, trim_str_arr)


class KEYS:
    ALIAS = 'alias'
    BIBCODE = 'bibcode'
    CLAIMED_TYPE = 'claimedtype'
    DISTINCTS = 'distinctfrom'
    DISCOVERY_DATE = 'discoverdate'
    ERRORS = 'errors'
    NAME = 'name'
    SOURCES = 'sources'
    SCHEMA = 'schema'
    URL = 'url'


class Entry(OrderedDict):

    # Whether or not this entry is a 'stub'.  Assume False
    _stub = False
    filename = None

    def __init__(self, name, stub=False):
        """Create a new `Entry` object with the given `name`.
        """
        # if not name:
        #     raise ValueError("New `Entry` objects must have a valid name!")
        self[KEYS.NAME] = name
        self._stub = stub
        self.filename = None
        return

    @classmethod
    def init_from_file(cls, name=None, path=None, clean=False):
        if name is None and path is None:
            raise ValueError("Either event `name` or `path` must be specified "
                             "to load event.")
        if name is not None and path is not None:
            raise ValueError("Either event `name` or `path` should be "
                             "specified, not both.")

        # If the path is given, use that to load from
        load_path = ''
        if path is not None:
            load_path = path
            name = ''
        # If the name is given, try to find a path for it
        else:
            repo_paths = get_repo_paths()
            for rep in repo_paths:
                filename = get_event_filename(name)
                newpath = os.path.join(rep, filename + '.json')
                if os.path.isfile(newpath):
                    load_path = newpath
                    break

        if load_path is None or not os.path.isfile(load_path):
            # FIX: is this warning worthy?
            return None

        # Create a new `EVENT` instance
        new_event = cls(name)
        # Fill it with data from json file
        new_event._load_data_from_json(load_path)

        if clean:
            new_event.clean()

        return new_event

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
        To convert a normal event into a stub (for example), overwrite the event
        in place, i.e.
        >>> events[name] = events[name].get_stub()

        """
        stub = EVENT(self[KEYS.NAME], stub=True)
        if KEYS.ALIAS in self.keys():
            stub[KEYS.ALIAS] = self[KEYS.ALIAS]
        return stub


class EVENT(Entry):
    """
    NOTE: OrderedDict data is just the `name` values from the JSON file.
          I.e. it does not include the highest nesting level
          { name: DATA }, it *just* includes DATA

    FIX: does this need to be `ordered`???
    FIX: check that no stored values are empty/invalid (delete key in that
         case?)
    FIX: distinguish between '.filename' and 'get_filename'

    sources
    -   All sources must have KEYS.NAME and 'alias' parameters
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

    def __init__(self, name, stub=False):
        super().__init__(name, stub=stub)

        # FIX: move this somewhere else (shouldnt be in each event)
        # Load source-name synonyms
        with open(FILENAME.SOURCE_SYNONYMS, 'r') as f:
            self._source_syns = json.loads(
                f.read(), object_pairs_hook=OrderedDict)
        return

    def add_source(self, srcname='', bibcode='', **src_kwargs):
        """Add a new source to this events KEYS.SOURCES list.

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

        # These are empty lists if no sources
        my_sources = self.get(KEYS.SOURCES, [])
        my_src_aliases = [src[KEYS.ALIAS] for src in my_sources]
        nsources = len(my_sources)

        # Try to find existing, matching source
        # -------------------------------------
        # If this source name already exists, return alias number
        try:
            my_src_names = [src[KEYS.NAME] for src in my_sources]
            name_idx = my_src_names.index(srcname)
            return my_src_aliases[name_idx]
        # `KeyError` from `KEYS.NAME` not existing, `ValueError` from
        # `srcname` not existing
        except (KeyError, ValueError):
            pass

        # If this bibcode already exists, return alias number
        try:
            my_src_bibs = [src[KEYS.BIBCODE] for src in my_sources]
            bib_idx = my_src_bibs.index(bibcode)
            return my_src_aliases[bib_idx]
        # `KeyError` from `KEYS.BIBCODE` not existing, `ValueError` from
        # `bibcode` not existing
        except (KeyError, ValueError):
            pass

        # Add new source that doesnt exist
        # --------------------------------
        source_alias = str(nsources + 1)
        new_src = OrderedDict()
        new_src[KEYS.NAME] = srcname
        if bibcode:
            new_src[KEYS.BIBCODE] = bibcode
        new_src[KEYS.ALIAS] = source_alias
        # Add in any additional arguments passed (e.g. url, acknowledgment,
        # etc)
        new_src.update({k: v for (k, v) in src_kwargs.items() if k})
        self.setdefault(KEYS.SOURCES, []).append(new_src)

        return source_alias

    def add_quantity(self, quantity, value, sources,
                     forcereplacebetter=False, derived='',
                     lowerlimit='', upperlimit='', error='', unit='',
                     kind='', extra='', probability=''):
        """
        """
        if not quantity:
            raise ValueError(self[KEYS.NAME] +
                             "'s quantity must be specified for "
                             "add_quantity.")
        if not sources:
            raise ValueError(self[KEYS.NAME] + "'s source must be specified for "
                             "quantity " +
                             quantity + ' before it is added.')
        if ((not isinstance(value, str) and
             (not isinstance(value, list) or not isinstance(value[0], str)))):
            raise ValueError(self[KEYS.NAME] + "'s Quantity " + quantity +
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
            raise ValueError(self[KEYS.NAME] + "'s quanta " + quantity +
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
            for df in self.get(KEYS.DISTINCTS, []):
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
        elif quantity == KEYS.CLAIMED_TYPE:
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
        if hasattr(self, KEYS.ERRORS):
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
        if KEYS.SCHEMA not in self.keys():
            self[KEYS.SCHEMA] = SCHEMA.URL
        # Make sure there is a name key in dict
        if KEYS.NAME not in self.keys() or len(self[KEYS.NAME]) == 0:
            raise ValueError("Event name is empty:\n\t{}".format(
                json.dumps(self, indent=2)))
        return

    def _get_save_path(self, bury=False):
        filename = get_event_filename(self[KEYS.NAME])

        # Put non-SNe in the boneyard
        if bury:
            outdir = str(PATH.REPO_BONEYARD)

        # Get normal repository save directory
        else:
            repo_folders = get_repo_folders()
            outdir = str(PATH.ROOT)
            if KEYS.DISCOVERY_DATE in self.keys():
                repo_years = get_repo_years(repo_folders)
                for r, year in enumerate(repo_years):
                    if int(self[KEYS.DISCOVERY_DATE][0]['value'].
                           split('/')[0]) <= year:
                        outdir = os.path.join(outdir, repo_folders[r])
                        break
            else:
                outdir = os.path.join(outdir, repo_folders[0])

        return outdir, filename

    def save(self, empty=False, bury=False, gz=False):
        outdir, filename = self._get_save_path(bury=bury)

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

    def get_source_by_alias(self, alias):
        for source in self.get(KEYS.SOURCES, []):
            if source['alias'] == alias:
                return source
        raise ValueError(
            "Source '{}': alias '{}' not found!".format(self[KEYS.NAME], alias))

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
                atels_dict = get_atels_dict()
                srcname = srcname.replace(
                    'ATEL', 'ATel').replace('Atel', 'ATel')
                srcname = srcname.replace(
                    'ATel #', 'ATel ').replace('ATel#', 'ATel')
                srcname = srcname.replace('ATel', 'ATel ')
                srcname = ' '.join(srcname.split())
                atelnum = srcname.split()[-1]
                if is_number(atelnum) and atelnum in atels_dict:
                    bibcode = atels_dict[atelnum]

            if srcname.upper().startswith('CBET'):
                cbets_dict = get_cbets_dict()
                srcname = srcname.replace('CBET', 'CBET ')
                srcname = ' '.join(srcname.split())
                cbetnum = srcname.split()[-1]
                if is_number(cbetnum) and cbetnum in cbets_dict:
                    bibcode = cbets_dict[cbetnum]

            if srcname.upper().startswith('IAUC'):
                iaucs_dict = get_iaucs_dict()
                srcname = srcname.replace('IAUC', 'IAUC ')
                srcname = ' '.join(srcname.split())
                iaucnum = srcname.split()[-1]
                if is_number(iaucnum) and iaucnum in iaucs_dict:
                    bibcode = iaucs_dict[iaucnum]

        for rep in self._source_syns:
            if srcname in self._source_syns[rep]:
                srcname = rep
                break

        return srcname, bibcode

    def clean(self):
        """

        FIX: instead of making changes in place to `dirty_event`, should a new
             event be created, values filled, then returned??
        FIX: currently will fail if no bibcode and no url
        """
        bibcodes = []

        # Rebuild the sources if they exist
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

        # Clean some legacy fields
        if 'aliases' in self and isinstance(self['aliases'], list):
            source = self.add_source(
                bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
            for alias in self['aliases']:
                self.add_quantity('alias', alias, source)
            del self['aliases']

        # FIX: should this be an error if false??
        if ((KEYS.DISTINCTS in self and
             isinstance(self[KEYS.DISTINCTS], list) and
             isinstance(self[KEYS.DISTINCTS][0], str))):
            distinctfroms = [x for x in self[KEYS.DISTINCTS]]
            del self[KEYS.DISTINCTS]
            source = self.add_source(
                bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
            for df in distinctfroms:
                self.add_quantity(KEYS.DISTINCTS, df, source)

        if 'errors' in self and \
                isinstance(self['errors'], list) and \
                'sourcekind' in self['errors'][0]:
            source = self.add_source(
                bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
            for err in self['errors']:
                self.add_quantity('error', err['quantity'], source,
                                  kind=err['sourcekind'], extra=err['id'])
            del self['errors']

        if not bibcodes:
            self.add_source(bibcode=OSC_BIBCODE,
                            srcname=OSC_NAME, url=OSC_URL, secondary=True)
            bibcodes = [OSC_BIBCODE]

        # Go through all keys in 'dirty' event
        for key in self.keys():
            if key in [KEYS.NAME, KEYS.SCHEMA, KEYS.SOURCES, KEYS.ERRORS]:
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
                       u_fluxdensity="", frequency="", u_frequency="", counts="",
                       e_counts="", nhmw="", photonindex="", unabsorbedflux="",
                       e_unabsorbedflux="", energy="", u_energy="",
                       e_lower_magnitude="", e_upper_magnitude="",
                       e_lower_time="", e_upper_time="", mcorrected=""):
        name = self[KEYS.NAME]
        if (not time and not host) or (not magnitude and not flux and not
                                       fluxdensity and not counts and not
                                       unabsorbedflux):
            warnings.warn(
                "Time or brightness not specified when adding photometry, not "
                "adding.")
            tprint('Name : "' + name + '", Time: "' + time + '", Band: "' +
                   band + '", AB magnitude: "' + magnitude + '"')
            return

        if (not host and not is_number(time)) or (not is_number(magnitude) and not
                                                  is_number(flux) and not
                                                  is_number(fluxdensity) and not
                                                  is_number(counts)):
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
                    warnings.warn('No error unit specified, not adding spectrum.')
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


def get_event_filename(name):
    return name.replace('/', '_')


def get_event_text(eventfile):
    import gzip
    if eventfile.split('.')[-1] == 'gz':
        with gzip.open(eventfile, 'rt') as f:
            filetext = f.read()
    else:
        with open(eventfile, 'r') as f:
            filetext = f.read()
    return filetext
