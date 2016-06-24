"""
"""
from cdecimal import Decimal
import codecs
from collections import OrderedDict
import json
from math import floor
import os
import warnings

from scripts import FILENAME, PATH, SCHEMA
from .constants import COMPRESS_ABOVE_FILESIZE, NON_SNE_PREFIXES, \
    OSC_BIBCODE, OSC_NAME, OSC_URL, REPR_BETTER_QUANTITY, TRAVIS_QUERY_LIMIT
from .funcs import copy_to_event, get_atels_dict, get_cbets_dict, get_iaucs_dict, \
    jd_to_mjd, name_clean, host_clean, radec_clean
from ..utils import get_repo_folders, get_repo_years, get_repo_paths, get_sig_digits, is_number, \
    pbar, pretty_num, tprint, zpad, repo_file_list


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


class EVENT(OrderedDict):
    """
    NOTE: OrderedDict data is just the `name` values from the JSON file.
          I.e. it does not include the highest nesting level { name: DATA }, it *just* includes
               DATA

    FIX: does this need to be `ordered`???
    FIX: check that no stored values are empty/invalid (delete key in that case?)
    FIX: be careful / remove duplicity between EVENT.name and EVENT['name'].
    FIX: distinguish between '.filename' and 'get_filename'

    sources
    -   All sources must have KEYS.NAME and 'alias' parameters
    -   FIX: is url required if no bibcode???
    -   FIX: consider changing 'alias' for each source to 'src_num' or something
    -   FIX: Make source aliases integers (instead of strings of integers)??
    -   FIX: have list of allowed 'source' parameters??
    -   FIX: create class for 'errors'
    -   FIX: class or list of valid quantities and units

    """

    name = ''
    filename = ''
    _source_syns = {}

    def __init__(self, name):
        self.name = name
        self[KEYS.NAME] = name
        # FIX: move this somewhere else (shouldnt be in each event)
        # Load source-name synonyms
        with open(FILENAME.SOURCE_SYNONYMS, 'r') as f:
            self._source_syns = json.loads(f.read(), object_pairs_hook=OrderedDict)
        return

    def load_data_from_json(self, fhand):
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
        if len(self.name) == 0:
            self.name = name
        # Warn if there is a name mismatch
        elif self.name.lower().strip() != name.lower().strip():
            warnings.warn("Object name '{}' does not match name in json: '{}'".format(
                self.name, name))

        self.check()
        return

    def add_source(self, srcname='', bibcode='', **src_kwargs):
        """Add a new source to this events KEYS.SOURCES list.

        FIX: if source already exists, should dictionary be updated to any new values??

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
        # Try to figure out each `srcname` or `bibcode` from the other, when only one given
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
        # `KeyError` from `KEYS.NAME` not existing, `ValueError` from `srcname` not existing
        except (KeyError, ValueError):
            pass

        # If this bibcode already exists, return alias number
        try:
            my_src_bibs = [src[KEYS.BIBCODE] for src in my_sources]
            bib_idx = my_src_bibs.index(bibcode)
            return my_src_aliases[bib_idx]
        # `KeyError` from `KEYS.BIBCODE` not existing, `ValueError` from `bibcode` not existing
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
        # Add in any additional arguments passed (e.g. url, acknowledgment, etc)
        new_src.update(src_kwargs)
        self.setdefault(KEYS.SOURCES, []).append(new_src)

        return source_alias

    def add_quantity(self, quantity, value, sources, forcereplacebetter=False, derived='',
                     lowerlimit='', upperlimit='', error='', unit='', kind='', extra='', probability=''):
        """
        """
        if not quantity:
            raise(ValueError(self.name + "'s quantity must be specified for add_quantity."))
        if not sources:
            raise(ValueError(self.name + "'s source must be specified for quantity " +
                quantity + ' before it is added.'))
        if ((not isinstance(value, str) and
             (not isinstance(value, list) or not isinstance(value[0], str)))):
            raise(ValueError(self.name + "'s Quantity " + quantity + " must be a string or an array of strings."))

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
            raise(ValueError(self.name + "'s quanta " + quantity + ' error value must be a number and positive.'))

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

        if quantity in ['velocity', 'redshift', 'ebv', 'lumdist', 'comovingdist']:
            if not is_number(svalue):
                return
        if quantity == 'host':
            if is_number(svalue):
                return
            if svalue.lower() in ['anonymous', 'anon.', 'anon', 'intergalactic']:
                return

            svalue = host_clean(svalue)

            if (not skind and ((svalue.lower().startswith('abell') and is_number(svalue[5:].strip())) or
                'cluster' in svalue.lower())):
                skind = 'cluster'
                if not skind and (is_abell or 'cluster' in svalue.lower()):
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
            (svalue, sunit) = radec_clean(svalue, quantity, unit = unit)
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
            quantaentry['probability'] = sprob
        if sunit:
            quanta_entry['unit'] = sunit
        if lowerlimit:
            quanta_entry['lowerlimit'] = lowerlimit
        if upperlimit:
            quanta_entry['upperlimit'] = upperlimit
        if derived:
            quantaentry['derived'] = derived
        if extra:
            quanta_entry['extra'] = extra
        if (forcereplacebetter or quantity in REPR_BETTER_QUANTITY) and len(my_quantity_list):
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
                        val_one = max(2, get_sig_digits(ctsplit[-1].lstrip('0')))
                        val_two = max(2, get_sig_digits(svsplit[-1].lstrip('0')))
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

    def _parse_srcname_bibcode(self, srcname, bibcode):
        # If no `srcname` is given, use `bibcode` after checking its validity
        if not srcname:
            if not bibcode:
                raise ValueError("`bibcode` must be specified if `srcname` is not.")
            if len(bibcode) != 19:
                raise ValueError("Bibcode '{}' must be exactly 19 characters long".format(bibcode))
            srcname = bibcode

        # If a `srcname` is given, try to set a `bibcode`
        elif not bibcode:
            if srcname.upper().startswith('ATEL'):
                atels_dict = get_atels_dict()
                srcname = srcname.replace('ATEL', 'ATel').replace('Atel', 'ATel')
                srcname = srcname.replace('ATel #', 'ATel ').replace('ATel#', 'ATel')
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

    def is_erroneous(self, field, sources):
        if hasattr(self, KEYS.ERRORS):
            my_errors = self['errors']
            for alias in sources.split(','):
                source = self.get_source_by_alias(alias)
                bib_err_values = [err['value'] for err in my_errors
                                  if err['kind'] == 'bibcode' and err['extra'] == field]
                if 'bibcode' in source and source['bibcode'] in bib_err_values:
                    return True

                name_err_values = [err['value'] for err in my_errors
                                   if err['kind'] == 'name' and err['extra'] == field]
                if 'name' in source and source['name'] in name_err_values:
                    return True

        return False

    def get_aliases(self, includename=True):
        # empty list if doesnt exist
        aliases = self.get(KEYS.ALIAS, [])
        if includename and self.name not in aliases:
            aliases = [self.name] + aliases
        return aliases

    def get_source_by_alias(self, alias):
        for source in self.get(KEYS.SOURCES, []):
            if source['alias'] == alias:
                return source
        raise ValueError("Source '{}': alias '{}' not found!".format(self.name, alias))

    def check(self):
        # Make sure there is a schema key in dict
        if KEYS.SCHEMA not in self.keys():
            self[KEYS.SCHEMA] = SCHEMA.URL
        # Make sure there is a name key in dict
        if KEYS.NAME not in self.keys():
            self[KEYS.NAME] = self.name
            if len(self[KEYS.NAME]) == 0:
                raise ValueError("Event name is empty:\n\t{}".format(json.dumps(self, indent=2)))
        # Make sure there is a name attribute in object
        if len(self.name) == '' and len(self[KEYS.NAME]) > 0:
            self.name = str(self[KEYS.NAME])

        return

    def _get_save_path(self, bury=False):
        filename = get_event_filename(self.name)

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
                    if int(self[KEYS.DISCOVERY_DATE][0]['value'].split('/')[0]) <= year:
                        outdir = os.path.join(outdir, repo_folders[r])
                        break
            else:
                outdir = os.path.join(outdir, repo_folders[0])

        return outdir, filename

    def save(self, empty=False, bury=False, gz=False):
        warnings.warn("Supressing 'stub' behavior in `EVENT.save`!")
        # if 'stub' in event_obj.keys():
        #     if not empty:
        #         continue
        #     else:
        #         del event_obj['stub']
        outdir, filename = self._get_save_path(bury=bury)

        # FIX: use 'dump' not 'dumps'
        jsonstring = json.dumps({self.name: self},
                                indent='\t', separators=(',', ':'), ensure_ascii=False)
        save_name = os.path.join(outdir, filename + '.json')
        with codecs.open(save_name, 'w', encoding='utf8') as sf:
            sf.write(jsonstring)

        return save_name

    def get_stub(self):
        stub = OrderedDict({KEYS.NAME: self[KEYS.NAME]})
        if KEYS.ALIAS in self.keys():
            stub[KEYS.ALIAS] = self[KEYS.ALIAS]
        return stub


def add_event(tasks, args, events, name, log, load=True, delete=True):
    """Find an existing event in, or add a new one to, the `events` dict.

    FIX: rename to `create_event`???

    Returns
    -------
    events : OrderedDict of EVENT objects
    newname : str
        Name of matching event found in `events`, or new event added to `events`
    """
    newname = name_clean(name)
    # If event already exists, return
    if newname in events:
        log.debug("`newname`: '{}' (name: '{}') already exists.".format(newname, name))
        return events, newname

    # If event is alias of another event, find and return that
    match_name = find_event_name_of_alias(events, newname)
    if match_name is not None:
        log.debug("`newname`: '{}' (name: '{}') already exist as alias for '{}'.".format(
            newname, name, match_name))
        return events, match_name

    # Load Event from file
    if load:
        loaded_event = load_event_from_file(events, args, tasks, log, name=newname, delete=delete)
        if loaded_event is not None:
            events[newname] = loaded_event
            log.debug("Added '{}', from '{}', to `events`".format(newname, loaded_event.filename))
            return events, loaded_event

    # Create new event
    new_event = EVENT(newname)
    new_event['schema'] = SCHEMA.URL
    log.log(log._LOADED, "Created new, empty event for '{}'".format(newname))
    # Add event to dictionary
    events[newname] = new_event

    return events, newname


def new_event(tasks, args, events, name, log, load = True, delete = True, loadifempty = True,
              refname = '', reference = '', url = '', bibcode = '', secondary = '', acknowledgment = ''):
    oldname = name
    events, name = add_event(tasks, args, events, name, load = load, delete = delete)
    source = events[name].add_source(name, bibcode=bibcode, refname=refname, reference=reference,
                                     url=url, secondary = secondary, acknowledgment = acknowledgment)
    events[name].add_quantity('alias', oldname, source)
    return events, name, source



def find_event_name_of_alias(events, alias):
    """Return the first event name with the given 'alias' included in its list of aliases.

    Returns
    -------
    name of matching event (str) or 'None' if no matches

    """
    for name, event in events.items():
        aliases = event.get_aliases()
        if alias in aliases:
            if (KEYS.DISTINCS not in event.keys()) or (alias not in event[KEYS.DISTINCS]):
                return name

    return None


def clean_event(dirty_event):
    """

    FIX: instead of making changes in place to `dirty_event`, should a new event be created,
         values filled, then returned??
    FIX: currently will fail if no bibcode and no url
    """
    bibcodes = []

    # Rebuild the sources
    try:
        old_sources = dirty_event[KEYS.SOURCES]
        del dirty_event[KEYS.SOURCES]
        for ss, source in enumerate(old_sources):
            if KEYS.BIBCODE in source:
                bibcodes.append(source[KEYS.BIBCODE])
                dirty_event.add_source(bibcode=source[KEYS.BIBCODE])
            else:
                dirty_event.add_source(srcname=source[KEYS.NAME], url=source[KEYS.URL])
    except KeyError:
        pass

    # Clean some legacy fields
    if 'aliases' in dirty_event and isinstance(dirty_event['aliases'], list):
        source = dirty_event.add_source(
            bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
        for alias in dirty_event['aliases']:
            dirty_event.add_quantity('alias', alias, source)
        del dirty_event['aliases']

    # FIX: should this be an error if false??
    if ((KEYS.DISTINCTS in dirty_event and isinstance(dirty_event[KEYS.DISTINCTS], list) and
         isinstance(dirty_event[KEYS.DISTINCTS][0], str))):
            distinctfroms = [x for x in dirty_event[KEYS.DISTINCTS]]
            del dirty_event[KEYS.DISTINCTS]
            source = dirty_event.add_source(
                bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
            for df in distinctfroms:
                dirty_event.add_quantity(KEYS.DISTINCTS, df, source)

    if (('errors' in dirty_event and isinstance(dirty_event['errors'], list) and
         'sourcekind' in dirty_event['errors'][0])):
            source = dirty_event.add_source(
                bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
            for err in dirty_event['errors']:
                dirty_event.add_quantity(
                    'error', err['quantity'], source, kind=err['sourcekind'], extra=err['id'])
            del dirty_event['errors']

    if not bibcodes:
        dirty_event.add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
        bibcodes = [OSC_BIBCODE]

    # Go through all keys in 'dirty' event
    for key in dirty_event.keys():
        if key in [KEYS.NAME, KEYS.SCHEMA, KEYS.SOURCES, KEYS.ERRORS]:
            pass
        elif key == 'photometry':
            for p, photo in enumerate(dirty_event['photometry']):
                if photo['u_time'] == 'JD':
                    dirty_event['photometry'][p]['u_time'] = 'MJD'
                    dirty_event['photometry'][p]['time'] = str(jd_to_mjd(Decimal(photo['time'])))
                if bibcodes and 'source' not in photo:
                    source = dirty_event.add_source(bibcode=bibcodes[0])
                    dirty_event['photometry'][p]['source'] = source
        else:
            for qi, quantity in enumerate(dirty_event[key]):
                if bibcodes and 'source' not in quantity:
                    source = dirty_event.add_source(bibcode=bibcodes[0])
                    dirty_event[key][qi]['source'] = source

    dirty_event.check()
    return dirty_event


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


def journal_events(tasks, args, events, stubs, log, clear=True, gz=False, bury=False):
    """Write all events in `events` to files, and clear.  Depending on arguments and `tasks`.
    """
    # FIX: store this somewhere instead of re-loading each time
    with open(FILENAME.NON_SNE_TYPES, 'r') as f:
        non_sne_types = json.loads(f.read(), object_pairs_hook=OrderedDict)
        non_sne_types = [x.upper() for x in non_sne_types]

    # Write it all out!
    # NOTE: this needs to use a `list` wrapper to allow modification of dictionary
    for name in list(events.keys()):
        if 'writeevents' in tasks:
            # See if this event should be buried

            # Bury non-SN events here if only claimed type is non-SN type, or if primary
            # name starts with a non-SN prefix.
            buryevent = False
            save_event = True
            ct_val = None
            if bury:
                if name.startswith(NON_SNE_PREFIXES):
                    log.debug("Killing '{}', non-SNe prefix.".format(name))
                    save_event = False
                else:
                    if KEYS.CLAIMED_TYPE in events[name]:
                        for ct in events[name][KEYS.CLAIMED_TYPE]:
                            up_val = ct['value'].upper()
                            if up_val not in non_sne_types and up_val != 'CANDIDATE':
                                buryevent = False
                                break
                            if up_val in non_sne_types:
                                buryevent = True
                                ct_val = ct['value']

                    if buryevent:
                        log.debug("Burying '{}', {}.".format(name, ct_val))

            if save_event:
                save_name = events[name].save(bury=buryevent)
                log.info("Saved {} to '{}'.".format(name.ljust(20), save_name))
                if gz and os.path.getsize(save_name) > COMPRESS_ABOVE_FILESIZE:
                    save_name = _compress_gz(save_name)
                    log.debug("Compressed '{}' to '{}'".format(name, save_name))
                    # FIX: use subprocess
                    outdir, filename = os.path.split(save_name)
                    filename = filename.split('.')[:-1]
                    os.system('cd ' + outdir + '; git rm ' + filename + '.json;' +
                              ' git add -f ' + filename + '.json.gz; cd ' + '../scripts')

        if clear:
            # Store stub of this object
            stubs.update({name: events[name].get_stub()})
            # Delete object
            del events[name]
            log.debug("Added stub for {}, deleted event.".format(name))

    return events, stubs


def load_event_from_file(events, args, tasks, log, name='', path='',
                         clean=False, delete=True, append=False):
    """

    FIX: currently will error if cant find a path from `name`
    """
    if not name and not path:
        raise ValueError('Either event `name` or `path` must be specified to load event.')

    if name and path:
        raise ValueError('Either event `name` or `path` should be specified, not both.')

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
        log.debug("No path found for name: '{}', path: '{}'".format(name, path))
        return None

    new_event = EVENT(name)
    new_event.load_data_from_json(load_path)

    # Delete old version
    if name in events:
        del events[name]

    if clean:
        new_event = clean_event(new_event)

    name = new_event[KEYS.NAME]
    log.log(log._LOADED, "Loaded {} from '{}'".format(name.ljust(20), load_path))

    # If this event loaded from an existing repo path and we will resave later, delete that version
    # FIX: have this check done to determine if `delete` is passed as True, when calling this func
    if 'writeevents' in tasks and delete and path_from_name:
        os.remove(path_from_name)
        log.debug("Deleted '{}'".format(path_from_name))

    return new_event


def load_stubs(tasks, args, events, log):
    """Load a 'stub' dictionary for all repository files.

    'stubs' only contain an event's name and aliases for cross-referencing.

    Returns
    -------
    stubs : OrderedDict,
        Dictionary of (name: stub), pairs where the stub contains the `KEYS.NAME` and `KEYS.ALIAS`
        parameters only.

    """
    currenttask = 'Loading event stubs'
    files = repo_file_list()
    stubs = OrderedDict()
    for fi in pbar(files, currenttask):
        fname = fi
        if '.gz' in fi:
            fname = _uncompress_gz(fi)
        name = os.path.basename(os.path.splitext(fname)[0]).replace('.json', '')
        new_event = load_event_from_file(events, args, tasks, log, path=fname, delete=False)
        stubs[name] = new_event.get_stub()
        log.debug("Added stub for '{}'".format(name))

    return stubs


def merge_duplicates(tasks, args, events):
    """Merge and remove duplicate events
    """
    if len(events) == 0:
        if args.update:
            import sys
            tprint('No sources changed, event files unchanged in update.')
            sys.exit(1)

        load_stubs(tasks, args, events)

    currenttask = 'Merging duplicate events'
    keys = list(sorted(list(events.keys())))
    for n1, name1 in enumerate(pbar(keys[:], currenttask)):
        if name1 not in events:
            continue
        # allnames1 = events[name1].get_aliases() + (['AT' + name1[2:]] if
        #     (name1.startswith('SN') and is_number(name1[2:6])) else [])
        allnames1 = set(events[name1].get_aliases(name1) + (['AT' + name1[2:]] if (name1.startswith('SN') and is_number(name1[2:6])) else []))

        for name2 in keys[n1+1:]:
            if name2 not in events or name1 == name2:
                continue
            # allnames2 = events[name2].get_aliases() + (['AT' + name2[2:]]
            #    if (name2.startswith('SN') and is_number(name2[2:6])) else [])
            allnames2 = set(events[name2].get_aliases(name2) + (['AT' + name2[2:]] if (name2.startswith('SN') and is_number(name2[2:6])) else []))
            if bool(allnames1 & allnames2):
                tprint("Found single event with multiple entries ('{}' and '{}'), merging.".format(
                    name1, name2))

                load1 = load_event_from_file(events, args, tasks, name1, delete=True)
                load2 = load_event_from_file(events, args, tasks, name2, delete=True)
                if load1 and load2:
                    priority1 = 0
                    priority2 = 0
                    for an in list(allnames1):
                        if len(an) >= 2 and an.startswith(('SN', 'AT')):
                            priority1 = priority1 + 1
                    for an in list(allnames2):
                        if len(an) >= 2 and an.startswith(('SN', 'AT')):
                            priority2 = priority2 + 1

                    if priority1 > priority2:
                        copy_to_event(events, name2, name1)
                        keys.append(name1)
                        del events[name2]
                    else:
                        copy_to_event(events, name1, name2)
                        keys.append(name2)
                        del events[name1]
                else:
                    print('Duplicate already deleted')
                journal_events(tasks, args, events)
        if args.travis and n1 > TRAVIS_QUERY_LIMIT:
            break

    return events


def set_preferred_names(tasks, args, events):
    currenttask = 'Setting preferred names'
    if not len(events):
        load_stubs(tasks, args, events)

    for ni, name in pbar(list(sorted(list(events.keys()))), currenttask):
        if name not in events:
            continue
        newname = ''
        aliases = events[name].get_aliases()
        if len(aliases) <= 1:
            continue
        if (name.startswith('SN') and ((is_number(name[2:6]) and not is_number(name[6:])) or
                                       (is_number(name[2:5]) and not is_number(name[5:])))):
            continue
        for alias in aliases:
            if (alias[:2] == 'SN' and ((is_number(alias[2:6]) and not is_number(alias[6:])) or
                                       (is_number(alias[2:5]) and not is_number(alias[5:])))):
                newname = alias
                break
        if not newname and 'discoverer' in events[name]:
            discoverer = ','.join([x['value'].upper() for x in events[name]['discoverer']])
            if 'ASAS' in discoverer:
                for alias in aliases:
                    if 'ASASSN' in alias.upper():
                        newname = alias
                        break
            if not newname and 'OGLE' in discoverer:
                for alias in aliases:
                    if 'OGLE' in alias.upper():
                        newname = alias
                        break
            if not newname and 'CRTS' in discoverer:
                for alias in aliases:
                    if True in [x in alias.upper() for x in ['CSS', 'MLS', 'SSS', 'SNHUNT']]:
                        newname = alias
                        break
            if not newname and 'PS1' in discoverer:
                for alias in aliases:
                    if 'PS1' in alias.upper():
                        newname = alias
                        break
            if not newname and 'PTF' in discoverer:
                for alias in aliases:
                    if 'PTF' in alias.upper():
                        newname = alias
                        break
            if not newname and 'GAIA' in discoverer:
                for alias in aliases:
                    if 'GAIA' in alias.upper():
                        newname = alias
                        break
        if not newname:
            for alias in aliases:
                # Always prefer another alias over PSN
                if name.startswith('PSN'):
                    newname = alias
                    break
        if newname and name != newname:
            # Make sure new name doesn't already exist
            if load_event_from_file(events, args, tasks, newname):
                continue
            if load_event_from_file(events, args, tasks, name, delete=True):
                tprint('Changing event name (' + name + ') to preferred name (' + newname + ').')
                events[newname] = events[name]
                events[newname][KEYS.NAME] = newname
                del events[name]
                journal_events(tasks, args, events)
        if args.travis and ni > TRAVIS_QUERY_LIMIT:
            break

    return events


def _compress_gz(fname):
    import shutil
    import gzip
    comp_fname = fname + '.gz'
    with open(fname, 'rb') as f_in, gzip.open(comp_fname, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(fname)
    return comp_fname


def _uncompress_gz(fname):
    import shutil
    import gzip
    uncomp_name = fname.replace('.gz', '')
    with gzip.open(fname, 'rb') as f_in, open(uncomp_name, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(fname)
    return uncomp_name
