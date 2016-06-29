import os
import sys
from collections import OrderedDict

from git import Repo

from scripts import FILENAME, PATH, SCHEMA

from .Entry import KEYS
from .importer.constants import COMPRESS_ABOVE_FILESIZE, TRAVIS_QUERY_LIMIT
from .importer.funcs import (event_attr_priority, name_clean, read_json_arr,
                             read_json_dict, uniq_cdl)
from .utils import is_number, logger, pbar, repo_file_list


class Catalog():
    """
    Object to hold the main catalog dictionary and other catalog globals.
    """

    def __init__(self, proto, args):
        # Store runtime arguments
        self.args = args

        # Set the catalog prototype class
        self.proto = proto

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

        # Create empty `entries` collection
        self.entries = OrderedDict()

        # Load repositories dictionary
        self.repos_dict = read_json_dict(FILENAME.REPOS)

        # Load auxiliary data
        self._load_aux()

        # Make sure repositories are cloned
        self._clone_repos()
        return

    def _clone_repos(self):
        all_repos = [self.repos_dict[x] for x in self.repos_dict]
        all_repos = [i for x in all_repos for i in x]
        for repo in pbar(all_repos):
            if not os.path.isdir(PATH.ROOT + "/" + repo):
                try:
                    self.log.warning(
                        'Cloning "' + repo + '" (only needs to be done ' +
                        'once, may take few minutes per repo).')
                    Repo.clone_from("git@github.com:astrocatalogs/" +
                                    repo + ".git",
                                    PATH.ROOT + "/" + repo)
                except:
                    self.log.error("CLONING '{}' INTERRUPTED".format(repo))
                    raise
                    sys.exit()
        return

    def _load_aux(self):
        # Create/Load auxiliary dictionaries
        self.nedd_dict = OrderedDict()
        self.bibauthor_dict = read_json_dict(FILENAME.BIBAUTHORS)
        self.biberror_dict = read_json_dict(FILENAME.BIBERRORS)
        self.extinctions_dict = read_json_dict(FILENAME.EXTINCT)
        # Create/Load auxiliary arrays
        self.nonsneprefixes_dict = read_json_arr(FILENAME.NON_SNE_PREFIXES)
        self.nonsnetypes = read_json_arr(FILENAME.NON_SNE_TYPES)
        return

    def add_entry(self, name, load=True, delete=True):
        """Find an existing entry in, or add a new one to, the `entries` dict.

        FIX: rename to `create_entry`???

        Returns
        -------
        entries : OrderedDict of Entry objects
        newname : str
            Name of matching entry found in `entries`, or new entry added to
            `entries`
        """
        self.log.debug("catalog.add_entry()")
        newname = name_clean(name)
        # If entry already exists, return
        if newname in self.entries:
            self.log.debug(
                "`newname`: '{}' (name: '{}') already exists.".
                format(newname, name))
            return newname

        # If entry is alias of another entry *in `entries`*, find and return
        # that
        match_name = self.find_entry_name_of_alias(self.entries, newname)
        if match_name is not None:
            self.log.debug(
                "`newname`: '{}' (name: '{}') already exist as alias for "
                "'{}'.".format(newname, name, match_name))
            return match_name

        # Load Event from file
        if load:
            loaded_entry = self.proto.init_from_file(name=newname)
            if loaded_entry is not None:
                self.entries[newname] = loaded_entry
                self.log.debug(
                    "Added '{}', from '{}', to `self.entries`"
                    .format(newname, loaded_entry.filename))
                # Delete source file, if desired
                if delete:
                    self._delete_entry_file(entry=loaded_entry)
                return newname

        # Create new entry
        new_entry = self.proto(self, newname)
        new_entry['schema'] = SCHEMA.URL
        self.log.log(self.log._LOADED,
                     "Created new entry for '{}'".format(newname))
        # Add entry to dictionary
        self.entries[newname] = new_entry
        return newname

    def delete_old_entry_files(self):
        if len(self.entries):
            err_str = "`delete_old_entry_files` with `entries` not empty!"
            self.log.error(err_str)
            raise RuntimeError(err_str)
        # Delete all old entry JSON files
        repo_files = repo_file_list()
        for rfil in pbar(repo_files, desc='Deleting old entries'):
            os.remove(rfil)
            self.log.debug("Deleted '{}'".format(os.path.split(rfil)[-1]))
        return

    def find_entry_name_of_alias(self, entries, alias):
        """Return the first entry name with the given 'alias' included in its
        list of aliases.

        Returns
        -------
        name of matching entry (str) or 'None' if no matches

        """
        for name, entry in entries.items():
            aliases = entry.get_aliases()
            if alias in aliases:
                if ((KEYS.DISTINCTS not in entry.keys()) or
                        (alias not in entry[KEYS.DISTINCTS])):
                    return name

        return None

    def copy_to_entry(self, fromname, destname):
        """

        Used by `merge_duplicates`
        """
        self.log.debug("Events.copy_to_entry()")
        self.log.info("Copy '{}' to '{}'".format(fromname, destname))
        newsourcealiases = {}
        keys = list(sorted(self.entries[fromname].keys(),
                           key=lambda xx: entry_attr_priority(xx)))

        if 'sources' in self.entries[fromname]:
            for source in self.entries[fromname]['sources']:
                newsourcealiases[source['alias']] = (self.entries[destname]
                                                     .add_source(
                    bibcode=source['bibcode'] if 'bibcode' in source else '',
                    srcname=source['name'] if 'name' in source else '',
                    reference=source['reference'] if
                    'reference' in source else '',
                    url=source['url'] if 'url' in source else ''))

        if 'errors' in self.entries[fromname]:
            for err in self.entries[fromname]['errors']:
                self.entries[destname].setdefault('errors', []).append(err)

        for key in keys:
            if key not in ['schema', 'name', 'sources', 'errors']:
                for item in self.entries[fromname][key]:
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
                        self.entries[destname].add_photometry(
                            u_time=item.get("u_time", ""),
                            time=item.get("time", ""),
                            e_time=item.get("e_time", ""),
                            telescope=item.get("telescope", ""),
                            instrument=item.get("instrument", ""),
                            band=item.get("band", ""),
                            magnitude=item.get("magnitude", ""),
                            e_magnitude=item.get("e_magnitude", ""),
                            source=sources,
                            upperlimit=item.get("upperlimit", ""),
                            system=item.get("system", ""),
                            observatory=item.get("observatory", ""),
                            observer=item.get("observer", ""),
                            host=item.get("host", ""),
                            survey=item.get("survey"), "")
                    elif key == 'spectra':
                        self.entries[destname].add_spectrum(
                            item.get("waveunit", ""),
                            item.get("fluxunit", ""),
                            data=item.get("data", ""),
                            u_time=item.get("u_time", ""),
                            time=item.get("time", ""),
                            instrument=item.get("instrument", ""),
                            deredshifted=item.get("deredshifted", ""),
                            dereddened=item.get("dereddened", ""),
                            errorunit=item.get("errorunit", ""),
                            source=sources, snr=item.get("snr", ""),
                            telescope=item.get("telescope", ""),
                            observer=item.get("observer", ""),
                            reducer=item.get("reducer", ""),
                            filename=item.get("filename", ""),
                            observatory=item.get("observatory"), "")
                    elif key == 'errors':
                        self.entries[destname].add_quantity(
                            key, item['value'], sources,
                            kind=item.get("kind", ""),
                            extra=item.get("extra"), "")
                    else:
                        self.entries[destname].add_quantity(
                            key, item['value'], sources,
                            error=item.get("error", ""),
                            unit=item.get("unit", ""),
                            probability=item.get("probability", ""),
                            kind=item.get("kind"), "")

        return

    def new_entry(self, name, load=True, delete=True,
                  loadifempty=True, srcname='', reference='', url='',
                  bibcode='', secondary='', acknowledgment=''):
        newname = self.add_entry(name, load=load, delete=delete)
        source = self.entries[newname].add_source(
            bibcode=bibcode, srcname=srcname, reference=reference, url=url,
            secondary=secondary, acknowledgment=acknowledgment)
        self.entries[newname].add_quantity('alias', name, source)
        return newname, source

    def merge_duplicates(self):
        """Merge and remove duplicate entries.

        Compares each entry ('name') in `stubs` to all later entries to check
        for duplicates in name or alias.  If a duplicate is found, they are
        merged and written to file.
        """
        self.log.debug("Events.merge_duplicates()")
        if len(self.entries) == 0:
            self.log.error("WARNING: `entries` is empty, loading stubs")
            if self.args.update:
                self.log.warning(
                    "No sources changed, entry files unchanged in update."
                    "  Skipping merge.")
                return
            entries = self.load_stubs()

        currenttask = 'Merging duplicate entries'

        keys = list(sorted(entries.keys()))
        for n1, name1 in enumerate(pbar(keys, currenttask)):
            allnames1 = set(entries[name1].get_aliases())
            if name1.startswith('SN') and is_number(name1[2:6]):
                allnames1 = allnames1.union(['AT' + name1[2:]])

            # Search all later names
            for name2 in keys[n1 + 1:]:
                allnames2 = set(entries[name2].get_aliases())
                if name2.startswith('SN') and is_number(name2[2:6]):
                    allnames2.union(['AT' + name2[2:]])

                # If there are any common names or aliases, merge
                if len(allnames1 & allnames2):
                    self.log.warning(
                        "Found single entry with multiple entries "
                        "('{}' and '{}'), merging.".format(name1, name2))

                    load1 = self.proto.init_from_file(name=name1, delete=True)
                    load2 = self.proto.init_from_file(name=name2, delete=True)
                    if load1 is not None and load2 is not None:
                        # Delete old files
                        self._delete_entry_file(entry=load1)
                        self._delete_entry_file(entry=load2)
                        priority1 = 0
                        priority2 = 0
                        for an in allnames1:
                            if an.startswith(('SN', 'AT')):
                                priority1 += 1
                        for an in allnames2:
                            if an.startswith(('SN', 'AT')):
                                priority2 += 1

                        if priority1 > priority2:
                            self.copy_to_entry(name2, name1)
                            keys.append(name1)
                            del entries[name2]
                        else:
                            self.copy_to_entry(name1, name2)
                            keys.append(name2)
                            del entries[name1]
                    else:
                        self.log.warning('Duplicate already deleted')

                    if len(entries) != 1:
                        self.log.error(
                            "WARNING: len(entries) = {}, expected 1.  "
                            "Still journaling...".format(len(entries)))
                    entries = self.journal_entries()

            if self.args.travis and n1 > TRAVIS_QUERY_LIMIT:
                break

        return entries

    def set_preferred_names(self):
        """Choose between each entries given name and its possible aliases for
        the best one.

        Highest preference goes to names of the form 'SN####AA'.
        Otherwise base the name on whichever survey is the 'discoverer'.

        FIX: create function to match SN####AA type names.
        """
        self.log.debug("Events.set_preferred_names()")

        if len(self.entries) == 0:
            self.log.error("WARNING: `entries` is empty, loading stubs")
            self.load_stubs()

        currenttask = 'Setting preferred names'
        for ni, name in pbar(list(sorted(self.entries.keys())), currenttask):
            newname = ''
            aliases = self.entries[name].get_aliases()
            # if there are no other options to choose from, skip
            if len(aliases) <= 1:
                continue
            # If the name is already in the form 'SN####AA' then keep using
            # that
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
            if not newname and 'discoverer' in self.entries[name]:
                discoverer = ','.join(
                    [x['value'].upper() for x in
                     self.entries[name]['discoverer']])
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
                if self.proto.init_from_file(name=newname):
                    self.log.error("WARNING: `newname` already exists... "
                                   "should do something about that...")
                    continue

                new_entry = self.proto.init_from_file(name=name)
                if new_entry is None:
                    self.log.error(
                        "Could not load `new_entry` with name '{}'"
                        .format(name))
                else:
                    self.log.info("Changing entry from name '{}' to preferred"
                                  " name '{}'".format(name, newname))
                    self._delete_entry_file(entry=new_entry)
                    self.entries[newname] = new_entry
                    self.entries[newname][KEYS.NAME] = newname
                    if name in self.entries:
                        self.log.error("WARNING: `name` = '{}' is in `entries` "
                                       "shouldnt happen?".format(name))
                        del self.entries[name]
                    self.journal_entries()

            if self.args.travis and ni > TRAVIS_QUERY_LIMIT:
                break

        return

    def load_stubs(self):
        """
        """
        currenttask = 'Loading entry stubs'
        files = repo_file_list()
        for fi in pbar(files, currenttask):
            fname = fi
            # FIX: should this be ``fi.endswith(``.gz')`` ?
            if '.gz' in fi:
                fname = _uncompress_gz(fi)
            name = os.path.basename(
                os.path.splitext(fname)[0]).replace('.json', '')
            new_entry = self.proto.init_from_file(path=fname, delete=False)
            # Make sure a non-stub entry doesnt already exist with this name
            if name in self.entries and not self.entries[name]._stub:
                err_str = (
                    "ERROR: non-stub entry already exists with name '{}'"
                    .format(name))
                self.log.error(err_str)
                raise RuntimeError(err_str)

            self.entries[name] = new_entry.get_stub()
            self.log.debug("Added stub for '{}'".format(name))

        return

    def _delete_entry_file(self, entry_name=None, entry=None):
        """Delete the file associated with the given entry.
        """
        if entry_name is None and entry is None:
            raise RuntimeError("Either `entry_name` or `entry` must be given.")
        elif entry_name is not None and entry is not None:
            raise RuntimeError("Cannot use both `entry_name` and `entry`.")

        if entry_name is not None:
            entry_filename = self.entries[entry_name]
        else:
            entry_name = entry[KEYS.NAME]
            entry_filename = entry.filename

        if self.args.write_entries:
            os.remove(entry_filename)
            self.log.info("Deleted entry '{}' file '{}'".format(
                entry_name, entry_filename))
        else:
            self.log.debug("Not deleting '{}' because `write_entries`"
                           " is Failse".format(entry_filename))

        return

    def journal_entries(self, clear=True, gz=False, bury=False,
                        write_stubs=False):
        """Write all entries in `entries` to files, and clear.  Depending on
        arguments and `tasks`.

        Iterates over all elements of `entries`, saving (possibly 'burying') and
        deleting.
        -   If ``clear == True``, then each element of `entries` is deleted, and
            a `stubs` entry is added
        """
        self.log.debug("Events.journal_entries()")
        # FIX: store this somewhere instead of re-loading each time

        # Write it all out!
        # NOTE: this needs to use a `list` wrapper to allow modification of
        # dict
        for name in list(self.entries.keys()):
            if self.args.write_entries:
                # If this is a stub and we aren't writing stubs, skip
                if self.entries[name]._stub and not write_stubs:
                    continue

                # Bury non-SN entries here if only claimed type is non-SN type,
                # or if primary name starts with a non-SN prefix.
                buryentry = False
                save_entry = True
                ct_val = None
                if bury:
                    if name.startswith(self.nonsneprefixes_dict):
                        self.log.debug(
                            "Killing '{}', non-SNe prefix.".format(name))
                        save_entry = False
                    else:
                        if KEYS.CLAIMED_TYPE in self.entries[name]:
                            for ct in self.entries[name][KEYS.CLAIMED_TYPE]:
                                up_val = ct['value'].upper()
                                if up_val not in self.non_sne_types and \
                                        up_val != 'CANDIDATE':
                                    buryentry = False
                                    break
                                if up_val in self.non_sne_types:
                                    buryentry = True
                                    ct_val = ct['value']

                        if buryentry:
                            self.log.debug(
                                "Burying '{}', {}.".format(name, ct_val))

                if save_entry:
                    save_name = self.entries[name].save(bury=buryentry)
                    self.log.info(
                        "Saved {} to '{}'.".format(name.ljust(20), save_name))
                    if (gz and os.path.getsize(save_name) >
                            COMPRESS_ABOVE_FILESIZE):
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
                self.entries[name] = self.entries[name].get_stub()
                self.log.debug("Entry for '{}' converted to stub".format(name))

        return

    def entry_exists(self, name):
        if name in self.entries:
            return True
        for ev in self.entries:
            if name in self.entries[ev].get_aliases():
                return True
        return False

    def count(self):
        full = 0
        stub = 0
        for ev in self.entries:
            if self.entries[ev]._stub:
                stub += 1
            else:
                full += 1
        return full, stub

    def get_current_task_str(self):
        return self.current_task.current_task(self.args)


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
