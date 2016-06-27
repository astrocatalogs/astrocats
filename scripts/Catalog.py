import os
import sys
from collections import OrderedDict
import json

from git import Repo

from scripts import PATH, SCHEMA, FILENAME

from .utils import logger, pbar
from .importer.funcs import (get_bibauthor_dict, get_biberror_dict,
                             get_extinctions_dict, get_repos_dict)

from .utils import (is_number, repo_file_list)
from .importer.constants import (COMPRESS_ABOVE_FILESIZE, NON_SNE_PREFIXES,
                                 TRAVIS_QUERY_LIMIT)
from .importer.funcs import (add_photometry, add_spectrum, event_attr_priority,
                             name_clean, null_field, uniq_cdl)
from .importer import Events
from .importer.Events import KEYS, EVENT


class Catalog():
    """
    Object to hold the main catalog dictionary and other catalog globals.
    """

    def __init__(self, args):
        # Store runtime arguments
        self.args = args
        # Load a logger object
        # Determine verbosity ('None' means use default)
        log_stream_level = None
        if args.debug:
            log_stream_level = logger.DEBUG
        elif args.verbose:
            log_stream_level = logger.INFO

        # Destination of log-file ('None' means no file)
        self.log = logger.get_logger(
            stream_level=log_stream_level, tofile=args.log_filename)

        # Load all dictionaries
        self._load_dicts()

        # Make sure repositories are cloned
        self._clone_repos()
        return

    def _clone_repos(self):
        all_repos = [self.repos_dict[x] for x in self.repos_dict]
        all_repos = [i for x in all_repos for i in x]
        for repo in pbar(all_repos):
            if not os.path.isdir(PATH.ROOT + "/" + repo):
                try:
                    self.log.warning('Cloning "' + repo +
                                     '" (only needs to be done ' +
                                     'once, may take few minutes per repo).')
                    Repo.clone_from("git@github.com:astrocatalogs/" +
                                    repo + ".git",
                                    PATH.ROOT + "/" + repo)
                except:
                    self.log.error("CLONING '{}' INTERRUPTED".format(repo))
                    raise
                    sys.exit()
        return

    def _load_dicts(self):
        # Create empty `events` collection
        self.events = OrderedDict()
        # Create/Load auxiliary dictionaries
        self.nedd_dict = OrderedDict()
        self.repos_dict = get_repos_dict()
        self.bibauthor_dict = get_bibauthor_dict()
        self.biberror_dict = get_biberror_dict()
        self.extinctions_dict = get_extinctions_dict()
        return

    def add_event(self, name, load=True, delete=True):
        """Find an existing event in, or add a new one to, the `events` dict.

        FIX: rename to `create_event`???

        Returns
        -------
        events : OrderedDict of EVENT objects
        newname : str
            Name of matching event found in `events`, or new event added to
            `events`
        """
        self.log.debug("Events.add_event()")
        newname = name_clean(name)
        # If event already exists, return
        if newname in self.events:
            self.log.debug(
                "`newname`: '{}' (name: '{}') already exists.".
                format(newname, name))
            return

        # If event is alias of another event *in `events`*, find and return that
        match_name = self.find_event_name_of_alias(self.events, newname)
        if match_name is not None:
            self.log.debug("`newname`: '{}' (name: '{}') already exist as alias for "
                           "'{}'.".format(newname, name, match_name))
            return match_name

        # Load Event from file
        if load:
            loaded_event = EVENT.init_from_file(name=newname)
            if loaded_event is not None:
                self.events[newname] = loaded_event
                self.log.debug("Added '{}', from '{}', to `self.events`".format(
                    newname, loaded_event.filename))
                # Delete source file, if desired
                if delete:
                    self._delete_event_file(event=loaded_event)
                return

        # Create new event
        new_event = Events.EVENT(newname)
        new_event['schema'] = SCHEMA.URL
        self.log.log(self.log._LOADED, "Created new event for '{}'".format(newname))
        # Add event to dictionary
        self.events[newname] = new_event
        return newname

    def delete_old_event_files(self):
        if len(self.events):
            err_str = "`delete_old_event_files` with `events` not empty!"
            self.log.error(err_str)
            raise RuntimeError(err_str)
        # Delete all old event JSON files
        repo_files = repo_file_list()
        for rfil in pbar(repo_files, desc='Deleting old events'):
            os.remove(rfil)
            self.log.debug("Deleted '{}'".format(os.path.split(rfil)[-1]))
        return

    def find_event_name_of_alias(self, events, alias):
        """Return the first event name with the given 'alias' included in its
        list of aliases.

        Returns
        -------
        name of matching event (str) or 'None' if no matches

        """
        for name, event in events.items():
            aliases = event.get_aliases()
            if alias in aliases:
                if ((KEYS.DISTINCTS not in event.keys()) or
                        (alias not in event[KEYS.DISTINCTS])):
                    return name

        return None

    def copy_to_event(self, fromname, destname):
        """

        Used by `merge_duplicates`
        """
        self.log.debug("Events.copy_to_event()")
        self.log.info("Copy '{}' to '{}'".format(fromname, destname))
        newsourcealiases = {}
        keys = list(sorted(self.events[fromname].keys(),
                           key=lambda xx: event_attr_priority(xx)))

        if 'sources' in self.events[fromname]:
            for source in self.events[fromname]['sources']:
                newsourcealiases[source['alias']] = self.events[destname].add_source(
                    bibcode=source['bibcode'] if 'bibcode' in source else '',
                    srcname=source['name'] if 'name' in source else '',
                    reference=source['reference'] if 'reference' in source else '',
                    url=source['url'] if 'url' in source else '')

        if 'errors' in self.events[fromname]:
            for err in self.events[fromname]['errors']:
                self.events[destname].setdefault('errors', []).append(err)

        for key in keys:
            if key not in ['schema', 'name', 'sources', 'errors']:
                for item in self.events[fromname][key]:
                    # isd = False
                    sources = []
                    if 'source' not in item:
                        ValueError("Item has no source!")
                    for sid in item['source'].split(','):
                        if sid == 'D':
                            sources.append('D')
                        elif sid in newsourcealiases:
                            sources.append(newsourcealiases[sid])
                        else:
                            ValueError("Couldn't find source alias!")
                    sources = uniq_cdl(sources)

                    if key == 'photometry':
                        add_photometry(
                            self.events, destname, u_time=null_field(item, "u_time"), time=null_field(item, "time"),
                            e_time=null_field(item, "e_time"), telescope=null_field(item, "telescope"),
                            instrument=null_field(item, "instrument"), band=null_field(item, "band"),
                            magnitude=null_field(item, "magnitude"), e_magnitude=null_field(item, "e_magnitude"),
                            source=sources, upperlimit=null_field(item, "upperlimit"), system=null_field(item, "system"),
                            observatory=null_field(item, "observatory"), observer=null_field(item, "observer"),
                            host=null_field(item, "host"), survey=null_field(item, "survey"))
                    elif key == 'spectra':
                        add_spectrum(
                            self.events, destname, null_field(item, "waveunit"), null_field(item, "fluxunit"),
                            data=null_field(item, "data"),
                            u_time=null_field(item, "u_time"), time=null_field(item, "time"),
                            instrument=null_field(item, "instrument"), deredshifted=null_field(item, "deredshifted"),
                            dereddened=null_field(item, "dereddened"), errorunit=null_field(item, "errorunit"),
                            source=sources, snr=null_field(item, "snr"),
                            telescope=null_field(item, "telescope"), observer=null_field(item, "observer"),
                            reducer=null_field(item, "reducer"), filename=null_field(item, "filename"),
                            observatory=null_field(item, "observatory"))
                    elif key == 'errors':
                        self.events[destname].add_quantity(
                            key, item['value'], sources,
                            kind=null_field(item, "kind"), extra=null_field(item, "extra"))
                    else:
                        self.events[destname].add_quantity(
                            key, item['value'], sources, error=null_field(item, "error"),
                            unit=null_field(item, "unit"), probability=null_field(item, "probability"),
                            kind=null_field(item, "kind"))

        return

    def new_event(self, name, load=True, delete=True,
                  loadifempty=True, srcname='', reference='', url='',
                  bibcode='', secondary='', acknowledgment=''):
        oldname = name
        self.add_event(name, load=load, delete=delete)
        source = self.events[name].add_source(
            bibcode=bibcode, srcname=srcname, reference=reference, url=url,
            secondary=secondary, acknowledgment=acknowledgment)
        self.events[name].add_quantity('alias', oldname, source)
        return name, source

    def merge_duplicates(self):
        """Merge and remove duplicate events.

        Compares each entry ('name') in `stubs` to all later entries to check for
        duplicates in name or alias.  If a duplicate is found, they are merged and
        written to file.
        """
        self.log.debug("Events.merge_duplicates()")
        if len(self.events) == 0:
            self.log.error("WARNING: `events` is empty, loading stubs")
            if self.args.update:
                self.log.warning("No sources changed, event files unchanged in update."
                                 "  Skipping merge.")
                return
            events = self.load_stubs()

        currenttask = 'Merging duplicate events'

        keys = list(sorted(events.keys()))
        for n1, name1 in enumerate(pbar(keys, currenttask)):
            allnames1 = set(events[name1].get_aliases())
            if name1.startswith('SN') and is_number(name1[2:6]):
                allnames1 = allnames1.union(['AT' + name1[2:]])

            # Search all later names
            for name2 in keys[n1+1:]:
                allnames2 = set(events[name2].get_aliases())
                if name2.startswith('SN') and is_number(name2[2:6]):
                    allnames2.union(['AT' + name2[2:]])

                # If there are any common names or aliases, merge
                if len(allnames1 & allnames2):
                    self.log.warning("Found single event with multiple entries "
                                     "('{}' and '{}'), merging.".format(name1, name2))

                    load1 = EVENT.init_from_file(name=name1, delete=True)
                    load2 = EVENT.init_from_file(name=name2, delete=True)
                    if load1 is not None and load2 is not None:
                        # Delete old files
                        self._delete_event_file(event=load1)
                        self._delete_event_file(event=load2)
                        priority1 = 0
                        priority2 = 0
                        for an in allnames1:
                            if an.startswith(('SN', 'AT')):
                                priority1 += 1
                        for an in allnames2:
                            if an.startswith(('SN', 'AT')):
                                priority2 += 1

                        if priority1 > priority2:
                            self.copy_to_event(name2, name1)
                            keys.append(name1)
                            del events[name2]
                        else:
                            self.copy_to_event(name1, name2)
                            keys.append(name2)
                            del events[name1]
                    else:
                        self.log.warning('Duplicate already deleted')

                    if len(events) != 1:
                        self.log.error("WARNING: len(events) = {}, expected 1.  "
                                       "Still journaling...".format(len(events)))
                    events = self.journal_events()

            if self.args.travis and n1 > TRAVIS_QUERY_LIMIT:
                break

        return events

    def set_preferred_names(self):
        """Choose between each events given name and its possible aliases for
        the best one.

        Highest preference goes to names of the form 'SN####AA'.
        Otherwise base the name on whichever survey is the 'discoverer'.

        FIX: create function to match SN####AA type names.
        """
        self.log.debug("Events.set_preferred_names()")

        if len(self.events) == 0:
            self.log.error("WARNING: `events` is empty, loading stubs")
            self.load_stubs()

        currenttask = 'Setting preferred names'
        for ni, name in pbar(list(sorted(self.events.keys())), currenttask):
            newname = ''
            aliases = self.events[name].get_aliases()
            # if there are no other options to choose from, skip
            if len(aliases) <= 1:
                continue
            # If the name is already in the form 'SN####AA' then keep using that
            if (name.startswith('SN') and
                ((is_number(name[2:6]) and not is_number(name[6:])) or
                 (is_number(name[2:5]) and not is_number(name[5:])))):
                continue
            # If one of the aliases is in the form 'SN####AA' then use that
            for alias in aliases:
                if (alias[:2] == 'SN' and
                    ((is_number(alias[2:6]) and not is_number(alias[6:])) or
                     (is_number(alias[2:5]) and not is_number(alias[5:])))):
                    newname = alias
                    break
            # Otherwise, name based on the 'discoverer' survey
            if not newname and 'discoverer' in self.events[name]:
                discoverer = ','.join([x['value'].upper()
                                       for x in self.events[name]['discoverer']])
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
                        if True in [x in alias.upper()
                                    for x in ['CSS', 'MLS', 'SSS', 'SNHUNT']]:
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
            # Always prefer another alias over PSN
            if not newname and name.startswith('PSN'):
                newname = aliases[0]
            if newname and name != newname:
                # Make sure new name doesn't already exist
                if EVENT.init_from_file(name=newname):
                    self.log.error("WARNING: `newname` already exists... "
                                   "should do something about that...")
                    continue

                new_event = EVENT.init_from_file(name=name)
                if new_event is None:
                    self.log.error("Could not load `new_event` with name '{}'".format(
                        name))
                else:
                    self.log.info("Changing event from name '{}' to preferred"
                                  " name '{}'".format(name, newname))
                    self._delete_event_file(event=new_event)
                    self.events[newname] = new_event
                    self.events[newname][KEYS.NAME] = newname
                    if name in self.events:
                        self.log.error("WARNING: `name` = '{}' is in `events` "
                                       "shouldnt happen?".format(name))
                        del self.events[name]
                    self.journal_events()

            if self.args.travis and ni > TRAVIS_QUERY_LIMIT:
                break

        return

    def load_stubs(self):
        """
        """
        currenttask = 'Loading event stubs'
        files = repo_file_list()
        for fi in pbar(files, currenttask):
            fname = fi
            # FIX: should this be ``fi.endswith(``.gz')`` ?
            if '.gz' in fi:
                fname = _uncompress_gz(fi)
            name = os.path.basename(
                os.path.splitext(fname)[0]).replace('.json', '')
            new_event = EVENT.init_from_file(path=fname, delete=False)
            # Make sure a non-stub event doesnt already exist with this name
            if name in self.events and not self.events[name]._stub:
                err_str = ("ERROR: non-stub event already exists with name '{}'"
                           "".format(name))
                self.log.error(err_str)
                raise RuntimeError(err_str)

            self.events[name] = new_event.get_stub()
            self.log.debug("Added stub for '{}'".format(name))

        return

    def _delete_event_file(self, event_name=None, event=None):
        """Delete the file associated with the given event.
        """
        if event_name is None and event is None:
            raise RuntimeError("Either `event_name` or `event` must be given.")
        elif event_name is not None and event is not None:
            raise RuntimeError("Cannot use both `event_name` and `event`.")

        if event_name is not None:
            event_filename = self.events[event_name]
        else:
            event_name = event[KEYS.NAME]
            event_filename = event.filename

        if self.args.write_events:
            os.remove(event_filename)
            self.log.info("Deleted event '{}' file '{}'".format(
                event_name, event_filename))
        else:
            self.log.debug("Not deleting '{}' because `write_events`"
                           " is Failse".format(event_filename))

        return

    def journal_events(self, clear=True, gz=False, bury=False, write_stubs=False):
        """Write all events in `events` to files, and clear.  Depending on
        arguments and `tasks`.

        Iterates over all elements of `events`, saving (possibly 'burying') and
        deleting.
        -   If ``clear == True``, then each element of `events` is deleted, and
            a `stubs` entry is added
        """
        self.log.debug("Events.journal_events()")
        # FIX: store this somewhere instead of re-loading each time
        with open(FILENAME.NON_SNE_TYPES, 'r') as f:
            non_sne_types = json.loads(f.read(), object_pairs_hook=OrderedDict)
            non_sne_types = [x.upper() for x in non_sne_types]

        # Write it all out!
        # NOTE: this needs to use a `list` wrapper to allow modification of dict
        for name in list(self.events.keys()):
            if self.args.write_events:
                # If this is a stub and we aren't writing stubs, skip
                if self.events[name]._stub and not write_stubs:
                    continue

                # Bury non-SN events here if only claimed type is non-SN type,
                # or if primary name starts with a non-SN prefix.
                buryevent = False
                save_event = True
                ct_val = None
                if bury:
                    if name.startswith(NON_SNE_PREFIXES):
                        self.log.debug("Killing '{}', non-SNe prefix.".format(name))
                        save_event = False
                    else:
                        if KEYS.CLAIMED_TYPE in self.events[name]:
                            for ct in self.events[name][KEYS.CLAIMED_TYPE]:
                                up_val = ct['value'].upper()
                                if up_val not in non_sne_types and \
                                        up_val != 'CANDIDATE':
                                    buryevent = False
                                    break
                                if up_val in non_sne_types:
                                    buryevent = True
                                    ct_val = ct['value']

                        if buryevent:
                            self.log.debug("Burying '{}', {}.".format(name, ct_val))

                if save_event:
                    save_name = self.events[name].save(bury=buryevent)
                    self.log.info("Saved {} to '{}'.".format(name.ljust(20), save_name))
                    if gz and os.path.getsize(save_name) > COMPRESS_ABOVE_FILESIZE:
                        save_name = _compress_gz(save_name)
                        self.log.debug("Compressed '{}' to '{}'".format(
                            name, save_name))
                        # FIX: use subprocess
                        outdir, filename = os.path.split(save_name)
                        filename = filename.split('.')[:-1]
                        os.system('cd ' + outdir + '; git rm ' + filename +
                                  '.json; git add -f ' + filename +
                                  '.json.gz; cd ' + '../scripts')

            if clear:
                self.events[name] = self.events[name].get_stub()
                self.log.debug("Entry for '{}' converted to stub".format(name))

        return

    def count(self):
        full = 0
        stub = 0
        for ev in self.events:
            if self.events[ev]._stub:
                stub += 1
            else:
                full += 1
        return full, stub


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
