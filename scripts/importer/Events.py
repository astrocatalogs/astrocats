"""
"""
from cdecimal import Decimal
from collections import OrderedDict
import json
from math import floor
import os
import warnings

from scripts import FILENAME
from .constants import OSC_BIBCODE, OSC_NAME, OSC_URL, REPR_BETTER_QUANTITY
from .funcs import copy_to_event, get_aliases, get_atels_dict, get_cbets_dict, get_iaucs_dict, \
    jd_to_mjd, load_stubs, name_clean
from ..utils import get_repo_paths, get_event_filename, get_sig_digits, is_number, \
    pbar, \
    pretty_num, tprint, zpad


class KEYS:
    ALIAS = 'alias'
    BIBCODE = 'bibcode'
    DISTINCTS = 'distinctfrom'
    ERRORS = 'errors'
    NAME = 'name'
    SOURCES = 'sourcs'
    URL = 'url'


class EVENT(OrderedDict):
    """
    NOTE: OrderedDict data is just the `name` values from the JSON file.
          I.e. it does not include the highest nesting level { name: DATA }, it *just* includes
               DATA

    FIX: does this need to be `ordered`???
    FIX: check that no stored values are empty/invalid (delete key in that case?)
    FIX: be careful / remove duplicity between EVENT.name and EVENT['name'].

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
        # Load source-name synonyms
        with open(FILENAME.SOURCE_SYNONYMS, 'r') as f:
            self._source_syns = json.loads(f.read(), object_pairs_hook=OrderedDict)
        return

    def load_data_from_json(self, fhand):
        """FIX: check for overwrite??
        """
        with open(fhand, 'r') as jfil:
            data = json.load(jfil, object_pairs_hook=OrderedDict)
            json.dump(data, open('dirty_event_raw', 'w'), indent=2)
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
        my_src_names = [src[KEYS.NAME] for src in my_sources]
        my_src_bibs = [src[KEYS.BIBCODE] for src in my_sources]
        nsources = len(my_sources)

        # Try to find existing, matching source
        # -------------------------------------
        # If this source name already exists, return alias number
        try:
            name_idx = my_src_names.index(srcname)
            return my_src_aliases[name_idx]
        except ValueError:
            pass

        # If this bibcode already exists, return alias number
        try:
            bib_idx = my_src_bibs.index(bibcode)
            return my_src_aliases[bib_idx]
        except ValueError:
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

    def add_quantity(self, quantity, value, sources, forcereplacebetter=False,
                     lowerlimit='', upperlimit='', error='', unit='', kind='', extra=''):
        if not quantity:
            raise ValueError('Quantity must be specified for add_quantity.')
        if not sources:
            raise ValueError('Source must be specified for quantity before it is added.')
        if ((not isinstance(value, str) and
             (not isinstance(value, list) or not isinstance(value[0], str)))):
            raise ValueError('Quantity must be a string or an array of strings.')

        if self.is_erroneous(quantity, sources):
            return None

        my_quantity_list = self.get(quantity, [])

        svalue = value.strip()
        serror = error.strip()
        skind = kind.strip()
        sunit = ''

        if not svalue or svalue == '--' or svalue == '-':
            return
        if serror and (not is_number(serror) or float(serror) < 0):
            raise ValueError('Quanta error value must be a number and positive.')

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
            if svalue.startswith('M ') and is_number(svalue[2:]):
                svalue.replace('M ', 'M', 1)
            svalue = svalue.strip("()").replace('  ', ' ', 1)
            svalue = svalue.replace("Abell", "Abell ", 1)
            svalue = svalue.replace("APMUKS(BJ)", "APMUKS(BJ) ", 1)
            svalue = svalue.replace("ARP", "ARP ", 1)
            svalue = svalue.replace("CGCG", "CGCG ", 1)
            svalue = svalue.replace("HOLM", "HOLM ", 1)
            svalue = svalue.replace("IC", "IC ", 1)
            svalue = svalue.replace("Intergal.", "Intergalactic", 1)
            svalue = svalue.replace("MCG+", "MCG +", 1)
            svalue = svalue.replace("MCG-", "MCG -", 1)
            svalue = svalue.replace("M+", "MCG +", 1)
            svalue = svalue.replace("M-", "MCG -", 1)
            svalue = svalue.replace("MGC ", "MCG ", 1)
            svalue = svalue.replace("Mrk", "MRK", 1)
            svalue = svalue.replace("MRK", "MRK ", 1)
            svalue = svalue.replace("NGC", "NGC ", 1)
            svalue = svalue.replace("PGC", "PGC ", 1)
            svalue = svalue.replace("SDSS", "SDSS ", 1)
            svalue = svalue.replace("UGC", "UGC ", 1)
            if len(svalue) > 4 and svalue.startswith("PGC "):
                svalue = svalue[:4] + svalue[4:].lstrip(" 0")
            if len(svalue) > 4 and svalue.startswith("UGC "):
                svalue = svalue[:4] + svalue[4:].lstrip(" 0")
            if len(svalue) > 5 and svalue.startswith(("MCG +", "MCG -")):
                svalue = svalue[:5] + '-'.join([x.zfill(2) for x in svalue[5:].strip().split("-")])
            if len(svalue) > 5 and svalue.startswith("CGCG "):
                svalue = svalue[:5] + '-'.join([x.zfill(3) for x in svalue[5:].strip().split("-")])
            if (((len(svalue) > 1 and svalue.startswith("E")) or
                 (len(svalue) > 3 and svalue.startswith('ESO')))):
                if svalue[0] == "E":
                    esplit = svalue[1:].split("-")
                else:
                    esplit = svalue[3:].split("-")
                if len(esplit) == 2 and is_number(esplit[0].strip()):
                    if esplit[1].strip()[0] == 'G':
                        parttwo = esplit[1][1:].strip()
                    else:
                        parttwo = esplit[1].strip()
                    if is_number(parttwo.strip()):
                        svalue = 'ESO ' + esplit[0].lstrip('0') + '-G' + parttwo.lstrip('0')
            svalue = ' '.join(svalue.split())

            is_abell = svalue.lower().startswith('abell') and is_number(svalue[5:].strip())
            if not skind and (is_abell or 'cluster' in svalue.lower()):
                skind = 'cluster'

        elif quantity == 'claimedtype':
            isq = False
            svalue = svalue.replace('young', '')
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
            if unit == 'floatdegrees':
                deg = float('%g' % Decimal(svalue))
                sig = get_sig_digits(svalue)
                if 'ra' in quantity:
                    flhours = deg / 360.0 * 24.0
                    hours = floor(flhours)
                    minutes = floor((flhours - hours) * 60.0)
                    seconds = (flhours * 60.0 - (hours * 60.0 + minutes)) * 60.0
                    if seconds > 60.0:
                        raise ValueError('Invalid seconds value for ' + quantity)
                    svalue = (str(hours).zfill(2) + ':' + str(minutes).zfill(2) + ':' +
                              zpad(pretty_num(seconds, sig=sig-1)))
                elif 'dec' in quantity:
                    fldeg = abs(deg)
                    degree = floor(fldeg)
                    minutes = floor((fldeg - degree) * 60.0)
                    seconds = (fldeg * 60.0 - (degree * 60.0 + minutes)) * 60.0
                    if seconds > 60.0:
                        raise ValueError('Invalid seconds value for ' + quantity)
                    svalue = '+' if deg >= 0.0 else '-'
                    svalue += (str(degree).strip('+-').zfill(2) + ':' +
                               str(minutes).zfill(2) + ':' + zpad(pretty_num(seconds, sig=sig-1)))
            elif unit == 'nospace' and 'ra' in quantity:
                svalue = (svalue[:2] + ':' + svalue[2:4] +
                          ((':' + zpad(svalue[4:])) if len(svalue) > 4 else ''))
            elif unit == 'nospace' and 'dec' in quantity:
                if svalue.startswith(('+', '-')):
                    svalue = (svalue[:3] + ':' + svalue[3:5] +
                              ((':' + zpad(svalue[5:])) if len(svalue) > 5 else ''))
                else:
                    svalue = ('+' + svalue[:2] + ':' + svalue[2:4] +
                              ((':' + zpad(svalue[4:])) if len(svalue) > 4 else ''))
            else:
                svalue = svalue.replace(' ', ':')
                if 'dec' in quantity:
                    valuesplit = svalue.split(':')
                    svalue = (('-' if valuesplit[0].startswith('-') else '+') +
                              valuesplit[0].strip('+-').zfill(2) +
                              (':' + valuesplit[1].zfill(2) if len(valuesplit) > 1 else '') +
                              (':' + zpad(valuesplit[2]) if len(valuesplit) > 2 else ''))

            if 'ra' in quantity:
                sunit = 'hours'
            elif 'dec' in quantity:
                sunit = 'degrees'

            # Correct case of arcseconds = 60.0.
            valuesplit = svalue.split(':')
            if len(valuesplit) == 3 and valuesplit[-1] in ["60.0", "60.", "60"]:
                svalue = (valuesplit[0] + ':' + str(Decimal(valuesplit[1]) +
                          Decimal(1.0)) + ':' + "00.0")

            # Strip trailing dots.
            svalue = svalue.rstrip('.')

        elif quantity == 'maxdate' or quantity == 'discoverdate':
            # Make sure month and day have leading zeroes
            sparts = svalue.split('/')
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
        if sunit:
            quanta_entry['unit'] = sunit
        if lowerlimit:
            quanta_entry['lowerlimit'] = lowerlimit
        if upperlimit:
            quanta_entry['upperlimit'] = upperlimit
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
                        if max(2, get_sig_digits(ctsplit[-1].lstrip('0'))) < max(2, get_sig_digits(svsplit[-1].lstrip('0'))):
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
                srcname = srcname.replace('ATel #', 'ATel ').replace('ATel#', 'ATel').replace('ATel', 'ATel ')
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

    def get_source_by_alias(self, alias):
        for source in self.get(KEYS.SOURCES, []):
            if source['alias'] == alias:
                return source
        raise ValueError("Source '{}': alias '{}' not found!".format(self.name, alias))

    def check(self):
        # Make sure there is a name key in dict
        if KEYS.NAME not in self.keys():
            self[KEYS.NAME] = self.name
            if len(self[KEYS.NAME]) == 0:
                raise ValueError("Event name is empty:\n\t{}".format(json.dumps(self, indent=2)))
        # Make sure there is a name attribute in object
        if len(self.name) == '' and len(self[KEYS.NAME]) > 0:
            self.name = str(self[KEYS.NAME])
        return


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
        if key in [KEYS.NAME, KEYS.SOURCES, KEYS.ERRORS]:
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


def load_event_from_file(events, args, tasks, log, name='', path='',
                         clean=False, delete=True, append=False):
    """

    FIX: currently will error if cant find a path from `name`
    """
    if not name and not path:
        raise ValueError('Either event `name` or `path` must be specified to load event.')

    if name and path:
        raise ValueError('Either event `name` or `path` should be specified, not both.')

    # If the path is given, load event from it
    path_from_name = ''
    if path:
        load_path = path
    else:
        repo_paths = get_repo_paths()
        for rep in repo_paths:
            filename = get_event_filename(name)
            newpath = os.path.join(rep, filename + '.json')
            if os.path.isfile(path):
                path_from_name = newpath
                break

        load_path = path_from_name

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
        # allnames1 = get_aliases(events, name1) + (['AT' + name1[2:]] if
        #     (name1.startswith('SN') and is_number(name1[2:6])) else [])
        allnames1 = get_aliases(events, name1)
        if name1.startswith('SN') and is_number(name1[2:6]):
            allnames1 += ['AT' + name1[2:]]

        for name2 in keys[n1+1:]:
            if name2 not in events or name1 == name2:
                continue
            # allnames2 = get_aliases(events, name2) + (['AT' + name2[2:]]
            #    if (name2.startswith('SN') and is_number(name2[2:6])) else [])
            allnames2 = get_aliases(events, name2)
            if name2.startswith('SN') and is_number(name2[2:6]):
                allnames2 += ['AT' + name2[2:]]
            if any(ii in allnames1 for ii in allnames2):
                tprint("Found single event with multiple entries ('{}' and '{}'), merging.".format(
                    name1, name2))

                load1 = load_event_from_file(events, args, tasks, name1, delete=True)
                load2 = load_event_from_file(events, args, tasks, name2, delete=True)
                if load1 and load2:
                    priority1 = 0
                    priority2 = 0
                    for an in allnames1:
                        if len(an) >= 2 and an.startswith(('SN', 'AT')):
                            priority1 = priority1 + 1
                    for an in allnames2:
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

    return events


def set_preferred_names(tasks, args, events):
    if not len(events):
        load_stubs(tasks, args, events)

    for name in list(sorted(list(events.keys()))):
        if name not in events:
            continue
        newname = ''
        aliases = get_aliases(events, name)
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

    return events


def write_all_events(events, args, empty=False, gz=False, bury=False):
    """Save all `events` to files.
    """
    import codecs
    from scripts import PATH
    from .. utils import get_event_filename, get_repo_folders, get_repo_years, is_number, tprint
    repo_folders = get_repo_folders()
    non_sne_types = None
    if bury:
        with open(FILENAME.NON_SNE_TYPES, 'r') as f:
            non_sne_types = json.loads(f.read(), object_pairs_hook=OrderedDict)
            non_sne_types = [x.upper() for x in non_sne_types]

    # Write it all out!
    for name in events:
        if 'stub' in events[name]:
            if not empty:
                continue
            else:
                del events[name]['stub']
        if args.verbose and not args.travis:
            tprint('Writing ' + name)
        filename = get_event_filename(name)

        # outdir = '../'
        outdir = str(PATH.ROOT)
        if 'discoverdate' in events[name]:
            repo_years = get_repo_years(repo_folders)
            for r, year in enumerate(repo_years):
                if int(events[name]['discoverdate'][0]['value'].split('/')[0]) <= year:
                    # outdir += repo_folders[r]
                    outdir = os.path.join(outdir, repo_folders[r])
                    break
        else:
            # outdir += str(repo_folders[0])
            outdir = os.path.join(outdir, repo_folders[0])

        # Delete non-SN events here without IAU designations (those with only banned types)
        if bury:
            buryevent = False
            nonsneprefixes = ('PNVJ', 'PNV J', 'OGLE-2013-NOVA')
            if name.startswith(nonsneprefixes):
                tprint('Burying ' + name + ', non-SNe prefix.')
                continue
            if 'claimedtype' in events[name] and not (name.startswith('SN') and is_number(name[2:6])):
                for ct in events[name]['claimedtype']:
                    if ct['value'].upper() not in non_sne_types and ct['value'].upper() != 'CANDIDATE':
                        buryevent = False
                        break
                    if ct['value'].upper() in non_sne_types:
                        buryevent = True
                if buryevent:
                    tprint('Burying ' + name + ' (' + ct['value'] + ').')
                    # outdir = '../sne-boneyard'
                    outdir = str(PATH.REPO_BONEYARD)  # os.path.join(PATH.ROOT, 'sne-boneyard')

        jsonstring = json.dumps({name: events[name]}, indent='\t', separators=(',', ':'), ensure_ascii=False)

        # path = outdir + '/' + filename + '.json'
        path = os.path.join(outdir, filename + '.json')
        with codecs.open(path, 'w', encoding='utf8') as f:
            f.write(jsonstring)

        if gz:
            if os.path.getsize(path) > 90000000:
                import shutil
                import gzip
                if not args.travis:
                    tprint('Compressing ' + name)
                with open(path, 'rb') as f_in, gzip.open(path + '.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                os.remove(path)
                os.system('cd ' + outdir + '; git rm ' + filename + '.json; git add -f ' + filename + '.json.gz; cd ' + '../scripts')

    return
