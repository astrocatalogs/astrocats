"""
"""
import codecs
import json
import os
import warnings
from collections import OrderedDict

from cdecimal import Decimal
from scripts import FILENAME, PATH, SCHEMA

from ..utils import (get_repo_folders, get_repo_paths, get_repo_years,
                     get_sig_digits, is_number, pbar, repo_file_list)
from .constants import (COMPRESS_ABOVE_FILESIZE, NON_SNE_PREFIXES, OSC_BIBCODE,
                        OSC_NAME, OSC_URL, REPR_BETTER_QUANTITY,
                        TRAVIS_QUERY_LIMIT)
from .funcs import (event_attr_priority, get_atels_dict, get_cbets_dict,
                    get_iaucs_dict, host_clean, jd_to_mjd, name_clean,
                    null_field, radec_clean, uniq_cdl)


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
    def init_from_file(cls, name=None, path=None,
                       clean=False, delete=True, append=False)
        if not name and not path:
            raise ValueError("Either event `name` or `path` must be specified "
                             "to load event.")
        if name and path:
            raise ValueError("Either event `name` or `path` should be "
                             "specified, not both.")

        # If the path is given, use that to load from
        path_from_name = ''
        if path:
            load_path = path
        # If the name is given, try to find a path for it
        else:
            repo_paths = get_repo_paths()
            for rep in repo_paths:
                filename = get_event_filename(name)
                newpath = os.path.join(rep, filename + '.json')
                if os.path.isfile(path):
                    path_from_name = newpath
                    break

            load_path = path_from_name

        if not load_path or not os.path.isfile(load_path):
            warnings.warn("No path found for name: '{}', path: '{}'".format(
                name, path))
            return None

        # Create a new `EVENT` instance
        new_event = cls(name)
        # Fill it with data from json file
        new_event.load_data_from_json(load_path, log)

        if clean:
            new_event.clean()

        return new_event

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

    def load_data_from_json(self, fhand, log):
        """FIX: check for overwrite??
        """
        log.debug("Events.EVENT.load_data_from_json()")
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


def load_event_from_file(events, args, tasks, log, name='', path='',
                         clean=False, delete=True, append=False):
    """

    FIX: currently will error if cant find a path from `name`
    """
    log.debug("Events.load_event_from_file()")
    if not name and not path:
        raise ValueError(
            'Either event `name` or `path` must be specified to load event.')

    if name and path:
        raise ValueError(
            'Either event `name` or `path` should be specified, not both.')

    # If the path is given, use that to load from
    path_from_name = ''
    if path:
        load_path = path
    # If the name is given, try to find a path for it
    else:
        repo_paths = get_repo_paths()
        for rep in repo_paths:
            filename = get_event_filename(name)
            newpath = os.path.join(rep, filename + '.json')
            if os.path.isfile(path):
                path_from_name = newpath
                break

        load_path = path_from_name

    if not load_path or not os.path.isfile(load_path):
        log.debug("No path found for name: '{}', path: '{}'".
                  format(name, path))
        return None

    new_event = EVENT(name)
    new_event.load_data_from_json(load_path, log)

    # Delete old version
    if name in events:
        del events[name]

    if clean:
        new_event = clean_event(new_event)

    name = new_event[KEYS.NAME]
    log.log(log._LOADED, "Loaded {} from '{}'".format(
        name.ljust(20), load_path))

    # If this event loaded from an existing repo path and we will resave
    # later, delete that version
    # FIX: have this check done to determine if `delete` is passed as True,
    # when calling this func
    if 'writeevents' in tasks and delete and path_from_name:
        os.remove(path_from_name)
        log.debug("Deleted '{}'".format(path_from_name))

    return new_event
