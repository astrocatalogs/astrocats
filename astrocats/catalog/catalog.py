"""Overarching catalog object for all open catalogs.
"""
import codecs
import importlib
import json
import os
import psutil
import sys
import warnings
from collections import OrderedDict
from glob import glob

from git import Repo

from astrocats.catalog.entry import KEYS
from astrocats.catalog.task import Task
from astrocats.catalog.source import SOURCE
from astrocats.catalog.utils import (compress_gz, is_integer, is_number,
                                     logger, pbar, read_json_dict,
                                     uncompress_gz, uniq_cdl)
from astrocats.supernovae.utils import name_clean


class Catalog:
    """Object to hold the main catalog dictionary and other catalog globals.
    """

    OSC_BIBCODE = '2016arXiv160501054G'
    OSC_NAME = 'The Open Supernova Catalog'
    OSC_URL = 'https://sne.space'

    ADS_BIB_URL = ("http://adsabs.harvard.edu/cgi-bin/nph-abs_connect?"
                   "db_key=ALL&version=1&bibcode=")

    TRAVIS_QUERY_LIMIT = 10
    COMPRESS_ABOVE_FILESIZE = 90000000   # bytes

    class PATHS:
        """Store and control catalog file-structure information.

        Individual catalogs must provide the below file structure.
        -   `repos.json`
        -   `tasks.json`

        """
        PATH_BASE = os.path.abspath(os.path.dirname(__file__))

        def __init__(self, catalog):
            self.catalog = catalog
            self.PATH_INPUT = os.path.join(self.PATH_BASE, 'input', '')
            self.PATH_OUTPUT = os.path.join(self.PATH_BASE, 'output', '')
            # critical datafiles
            self.REPOS = os.path.join(self.PATH_INPUT, 'repos.json')
            self.TASK_LIST = os.path.join(self.PATH_INPUT, 'tasks.json')
            self.repos_dict = read_json_dict(self.REPOS)

        def _get_repo_file_list(self, repo_folders, normal=True, bones=True):
            """Get filenames for files in each repository, `boneyard` optional.
            """
            # repo_folders = get_repo_output_folders()
            files = []
            for rep in repo_folders:
                if 'boneyard' not in rep and not normal:
                    continue
                if not bones and 'boneyard' in rep:
                    continue
                these_files = glob(rep + "/*.json") + glob(rep + "/*.json.gz")
                self.catalog.log.debug("Found {} files in '{}'".format(
                    len(these_files), rep))
                files += these_files

            return files

        def get_repo_boneyard(self):
            bone_path = self.repos_dict['boneyard']
            try:
                bone_path = bone_path[0]
            except TypeError:
                pass
            bone_path = os.path.join(self.PATH_OUTPUT, bone_path, '')
            return bone_path

    class SCHEMA:
        HASH = ''
        URL = ''

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

        # Instantiate PATHS
        self.PATHS = self.PATHS(self)

        # Load repos dictionary (required)
        self.repos_dict = read_json_dict(self.PATHS.REPOS)

        # Create empty `entries` collection
        self.entries = OrderedDict()
        return

    def import_data(self):
        """Run all of the import tasks.

        This is executed by the 'scripts.main.py' when the module is run as an
        executable. This can also be run as a method, in which case default
        arguments are loaded, but can be overriden using `**kwargs`.
        """

        tasks_list = self.load_task_list()
        warnings.filterwarnings(
            'ignore', r'Warning: converting a masked element to nan.')
        # FIX
        warnings.filterwarnings(
            'ignore', category=DeprecationWarning)

        # Delete all old (previously constructored) output files
        if self.args.delete_old:
            self.log.warning("Deleting all old entry files.")
            self.delete_old_entry_files()

        prev_priority = 0
        prev_task_name = ''
        # for task, task_obj in tasks_list.items():
        for task_name, task_obj in tasks_list.items():
            if not task_obj.active:
                continue
            self.log.warning("Task: '{}'".format(task_name))

            nice_name = task_obj.nice_name
            mod_name = task_obj.module
            func_name = task_obj.function
            priority = task_obj.priority

            # Make sure things are running in the correct order
            if priority < prev_priority and priority > 0:
                raise RuntimeError(
                    "Priority for '{}': '{}', less than prev,"
                    "'{}': '{}'.\n{}"
                    .format(task_name, priority, prev_task_name, prev_priority,
                            task_obj))

            self.log.debug("\t{}, {}, {}, {}".format(
                nice_name, priority, mod_name, func_name))
            mod = importlib.import_module('.' + mod_name, package='astrocats')
            self.current_task = task_obj
            getattr(mod, func_name)(self)

            num_events, num_stubs = self.count()
            self.log.warning("Task finished.  Events: {},  Stubs: {}".format(
                num_events, num_stubs))
            self.journal_entries()
            num_events, num_stubs = self.count()
            self.log.warning("Journal finished.  Events: {}, Stubs: {}".format(
                num_events, num_stubs))

            prev_priority = priority
            prev_task_name = task_name

        process = psutil.Process(os.getpid())
        memory = process.memory_info().rss
        print('Memory used (MBs): '
              '{:,}'.format(memory / 1024. / 1024.))
        return

    def load_task_list(self):
        """Load the list of tasks in this catalog's 'input/tasks.json' file.

        A `Task` object is created for each entry, with the parameters filled
        in. These are placed in an OrderedDict, sorted by the `priority`
        parameter, with positive values and then negative values,
            e.g. [0, 2, 10, -10, -1].
        """

        # Make sure appropriate command-line arguments are used

        # Dont allow both a 'min' and 'max' task priority
        if ((self.args.min_task_priority is not None and
             self.args.max_task_priority is not None)):
            raise ValueError("Can only use *either* 'min' *or* 'max' priority")

        # Load tasks data from input json file
        def_task_list_filename = self.PATHS.TASK_LIST
        self.log.debug(
            "Loading task-list from '{}'".format(def_task_list_filename))
        data = json.load(open(def_task_list_filename, 'r'))

        # Make sure 'active' modification lists are all valid
        args_lists = [self.args.args_task_list,
                      self.args.yes_task_list, self.args.no_task_list]
        args_names = ['--tasks', '--yes', '--no']
        for arglist, lname in zip(args_lists, args_names):
            if arglist is not None:
                for tname in arglist:
                    if tname not in data.keys():
                        raise ValueError(
                            "Value '{}' in '{}' list does not match"
                            " any tasks".format(tname, lname))

        # Create `Task` objects for each element in the tasks data file
        tasks = {}
        for key, val in data.items():
            tasks[key] = Task(name=key, **val)

        # Process min/max priority specification ('None' if none given)
        min_priority = _get_task_priority(tasks, self.args.min_task_priority)
        max_priority = _get_task_priority(tasks, self.args.max_task_priority)
        task_groups = self.args.task_groups
        if task_groups is not None:
            if not isinstance(task_groups, list):
                task_groups = [task_groups]

        # Iterate over all tasks to determine which should be (in)active
        # --------------------------------------------------------------
        for key in tasks:
            # If specific list of tasks is given, make only those active
            if self.args.args_task_list is not None:
                if key in self.args.args_task_list:
                    tasks[key].active = True
                else:
                    tasks[key].active = False

            # Only run tasks above minimum priority
            # (doesn't modify negtive priority tasks)
            if min_priority is not None and tasks[key].priority >= 0:
                tasks[key].active = False
                if tasks[key].priority >= min_priority:
                    tasks[key].active = True

            # Only run tasks below maximum priority
            # (doesnt modify negative priority tasks)
            if max_priority is not None and tasks[key].priority >= 0:
                tasks[key].active = False
                if tasks[key].priority <= max_priority:
                    tasks[key].active = True

            # Set 'yes' tasks to *active*
            if self.args.yes_task_list is not None:
                if key in self.args.yes_task_list:
                    tasks[key].active = True
            # Set 'no' tasks to *inactive*
            if self.args.no_task_list is not None:
                if key in self.args.no_task_list:
                    tasks[key].active = False
            # Set tasks in target 'groups' to *active*
            if task_groups is not None and tasks[key].groups is not None:
                # Go through each group defined in the command line
                for given_group in task_groups:
                    # If this task is a member of any of those groups
                    if given_group in tasks[key].groups:
                        tasks[key].active = True
                        break

        # Sort entries as positive values, then negative values
        #    [0, 1, 2, 2, 10, -100, -10, -1]
        # Tuples are sorted by first element (here: '0' if positive), then
        # second (here normal order)
        tasks = OrderedDict(sorted(tasks.items(), key=lambda t: (
            t[1].priority < 0, t[1].priority, t[1].name)))

        names_act = []
        names_inact = []
        for key, val in tasks.items():
            if val.active:
                names_act.append(key)
            else:
                names_inact.append(key)

        self.log.info("Active Tasks:\n\t" + ", ".join(nn for nn in names_act))
        self.log.debug("Inactive Tasks:\n\t" +
                       ", ".join(nn for nn in names_inact))
        return tasks

    def _clone_repos(self, all_repos):
        """Given a list of repositories, make sure they're all cloned.

        Should be called from the subclassed `Catalog` objects, passed a list
        of specific repository names.

        Arguments
        ---------
        all_repos : list of str
            *Absolute* path specification of each target repository.

        """
        for repo in pbar(all_repos):
            if not os.path.isdir(repo):
                try:
                    repo_name = os.path.split(repo)[-1]
                    self.log.warning(
                        'Cloning "' + repo + '" (only needs to be done ' +
                        'once, may take few minutes per repo).')
                    Repo.clone_from("git@github.com:astrocatalogs/" +
                                    repo_name + ".git", repo)
                except:
                    self.log.error("CLONING '{}' INTERRUPTED".format(repo))
                    raise
                    sys.exit()
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
        newname = name_clean(name)
        # If entry already exists, return
        if newname in self.entries:
            self.log.debug(
                "`newname`: '{}' (name: '{}') already exists.".
                format(newname, name))
            # If this is a stub, we need to continue, possibly load file
            if self.entries[newname]._stub:
                self.log.debug("'{}' is a stub".format(newname))
            # If a full (non-stub) event exists, return its name
            else:
                self.log.debug("'{}' is not a stub, returning".format(newname))
                return newname

        # If entry is alias of another entry in `entries`, find and return that
        match_name = self.find_entry_name_of_alias(self.entries, newname)
        if match_name is not None:
            self.log.debug(
                "`newname`: '{}' (name: '{}') already exists as alias for "
                "'{}'.".format(newname, name, match_name))
            newname = match_name

        # Load Event from file
        if load:
            loaded_entry = self.proto.init_from_file(self, name=newname)
            if loaded_entry is not None:
                self.entries[newname] = loaded_entry
                self.log.debug(
                    "Added '{}', from '{}', to `self.entries`".format(
                        newname, loaded_entry.filename))
                # Delete source file, if desired
                if delete:
                    self._delete_entry_file(entry=loaded_entry)
                return newname

        # If we match an existing event, return that
        if match_name is not None:
            return match_name

        # Create new entry
        new_entry = self.proto(self, newname)
        new_entry[self.proto._KEYS.SCHEMA] = self.SCHEMA.URL
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
        repo_files = self.PATHS.get_repo_output_file_list()
        for rfil in pbar(repo_files, desc='Deleting old entries'):
            os.remove(rfil)
            self.log.debug("Deleted '{}'".format(os.path.split(rfil)[-1]))
        return

    def get_preferred_name(self, name):
        if name not in self.entries:
            # matches = []
            for entry in self.entries:
                aliases = self.entries[entry].get_aliases()
                if len(aliases) > 1 and name in aliases:
                    return entry
            return name
        else:
            return name

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
                if ((KEYS.DISTINCT_FROM not in entry.keys()) or
                        (alias not in entry[KEYS.DISTINCT_FROM])):
                    return name

        return None

    def copy_to_entry(self, fromname, destname):
        """

        Used by `merge_duplicates`
        """
        self.log.info("Copy '{}' to '{}'".format(fromname, destname))
        newsourcealiases = {}

        if self.proto._KEYS.SOURCES in self.entries[fromname]:
            for source in self.entries[fromname][self.proto._KEYS.SOURCES]:
                alias = source.pop(SOURCE.ALIAS)
                newsourcealiases[alias] = source

        if self.proto._KEYS.ERRORS in self.entries[fromname]:
            for err in self.entries[fromname][self.proto._KEYS.ERRORS]:
                self.entries[destname].setdefault(
                    self.proto._KEYS.ERRORS, []).append(err)

        for key in self.entries[fromname].keys():
            if key in ['schema', 'name', 'sources', 'errors']:
                continue
            for item in self.entries[fromname][key]:
                # isd = False
                if 'source' not in item:
                    ValueError("Item has no source!")

                nsid = []
                for sid in item['source'].split(','):
                    if sid in newsourcealiases:
                        source = newsourcealiases[sid]
                        nsid.append(self.entries[destname]
                                    .add_source(**source))
                    else:
                        ValueError("Couldn't find source alias!")

                item['source'] = uniq_cdl(nsid)

                if key == 'photometry':
                    self.entries[destname].add_photometry(**item)
                elif key == 'spectra':
                    self.entries[destname].add_spectrum(**item)
                elif key == 'errors':
                    self.entries[destname].add_error(**item)
                else:
                    self.entries[destname].add_quantity(quantity=key, **item)

        return

    def new_entry(self, name, load=True, delete=True,
                  loadifempty=True, srcname='', reference='', url='',
                  bibcode='', secondary=False, acknowledgment=''):
        newname = self.add_entry(name, load=load, delete=delete)
        source = self.entries[newname].add_source(
            bibcode=bibcode, name=srcname, reference=reference, url=url,
            secondary=secondary, acknowledgment=acknowledgment)
        self.entries[newname].add_quantity('alias', name, source)
        return newname, source

    def merge_duplicates(self):
        """Merge and remove duplicate entries.

        Compares each entry ('name') in `stubs` to all later entries to check
        for duplicates in name or alias.  If a duplicate is found, they are
        merged and written to file.
        """
        if len(self.entries) == 0:
            self.log.error("WARNING: `entries` is empty, loading stubs")
            if self.args.update:
                self.log.warning(
                    "No sources changed, entry files unchanged in update."
                    "  Skipping merge.")
                return
            self.entries = self.load_stubs()

        task_str = self.get_current_task_str()

        keys = list(sorted(self.entries.keys()))
        for n1, name1 in enumerate(pbar(keys, task_str)):
            allnames1 = set(self.entries[name1].get_aliases())

            # Search all later names
            for name2 in keys[n1 + 1:]:
                allnames2 = set(self.entries[name2].get_aliases())

                # If there are any common names or aliases, merge
                if len(allnames1 & allnames2):
                    self.log.warning(
                        "Found single entry with multiple entries "
                        "('{}' and '{}'), merging.".format(name1, name2))

                    load1 = self.proto.init_from_file(
                        self, name=name1)
                    load2 = self.proto.init_from_file(
                        self, name=name2)
                    if load1 is not None and load2 is not None:
                        # Delete old files
                        self._delete_entry_file(entry=load1)
                        self._delete_entry_file(entry=load2)
                        priority1 = 0
                        priority2 = 0

                        if priority1 > priority2:
                            self.entries[name1] = load1
                            self.copy_to_entry(name2, name1)
                            keys.append(name1)
                            del self.entries[name2]
                        else:
                            self.entries[name1] = load2
                            self.copy_to_entry(name1, name2)
                            keys.append(name2)
                            del self.entries[name1]
                    else:
                        self.log.warning('Duplicate already deleted')

                    if len(self.entries) != 1:
                        self.log.error(
                            "WARNING: len(entries) = {}, expected 1.  "
                            "Still journaling...".format(len(self.entries)))
                    self.journal_entries()

            if self.args.travis and n1 > self.TRAVIS_QUERY_LIMIT:
                break

    def sanitize(self):
        task_str = self.get_current_task_str()
        for name in pbar(list(sorted(self.entries.keys())), task_str):
            self.add_entry(name)
            self.journal_entries(bury=True, final=True)

    def set_preferred_names(self):
        """Choose between each entries given name and its possible aliases for
        the best one.

        Highest preference goes to names of the form 'SN####AA'.
        Otherwise base the name on whichever survey is the 'discoverer'.

        FIX: create function to match SN####AA type names.
        """
        pass

    def load_stubs(self):
        """
        """
        currenttask = 'Loading entry stubs'
        files = self.PATHS.get_repo_output_file_list()
        for fi in pbar(files, currenttask):
            fname = fi
            # FIX: should this be ``fi.endswith(``.gz')`` ?
            if '.gz' in fi:
                fname = uncompress_gz(fi)
            name = os.path.basename(
                os.path.splitext(fname)[0]).replace('.json', '')
            new_entry = self.proto.init_from_file(
                self, path=fname, delete=False)
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
            entry = self.entries[entry_name]
        else:
            entry_name = entry[KEYS.NAME]

        outdir, filename = entry._get_save_path()
        # FIX: do we also need to check for gzipped files??
        entry_filename = os.path.join(outdir, filename + '.json')
        # entry_filename = entry.filename

        if self.args.write_entries:
            self.log.info("Deleting entry file '{}' of entry '{}'".format(
                entry_filename, entry_name))
            if not os.path.exists(entry_filename):
                self.log.error("Filename '{}' does not exist".format(
                    entry_filename))
            os.remove(entry_filename)
        else:
            self.log.debug("Not deleting '{}' because `write_entries`"
                           " is False".format(entry_filename))

        return

    def should_bury(self, name):
        return (False, True)

    def journal_entries(self, clear=True, gz=False, bury=False,
                        write_stubs=False, final=False):
        """Write all entries in `entries` to files, and clear.  Depending on
        arguments and `tasks`.

        Iterates over all elements of `entries`, saving (possibly 'burying')
        and deleting.
        -   If ``clear == True``, then each element of `entries` is deleted,
            and a `stubs` entry is added
        """
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
                bury_entry = False
                save_entry = True
                if bury:
                    (bury_entry, save_entry) = self.should_bury(name)

                if save_entry:
                    save_name = self.entries[name].save(bury=bury_entry,
                                                        final=final)
                    self.log.info(
                        "Saved {} to '{}'.".format(name.ljust(20), save_name))
                    if (gz and os.path.getsize(save_name) >
                            self.COMPRESS_ABOVE_FILESIZE):
                        save_name = compress_gz(save_name)
                        self.log.debug("Compressed '{}' to '{}'".format(
                            name, save_name))
                        # FIX: use subprocess
                        outdir, filename = os.path.split(save_name)
                        filename = filename.split('.')[:-1]
                        os.system('cd ' + outdir + '; git rm ' + filename +
                                  '.json; git add -f ' + filename +
                                  '.json.gz; cd ' + self.PATHS.PATH_BASE)

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

    def get_current_task_repo(self):
        return self.current_task._get_repo_path(self.PATHS.PATH_BASE)

    def has_task(self, task):
        return task in self.tasks and (not self.args.update or
                                       self.tasks[task]['update'])

    def load_cached_url(self, url, filepath, timeout=120, write=True,
                        failhard=False):
        from hashlib import md5
        filemd5 = ''
        filetxt = ''
        if not self.args.refresh and os.path.isfile(filepath):
            with codecs.open(filepath, 'r', encoding='utf8') as f:
                filetxt = f.read()
                if self.args.update:
                    filemd5 = md5(filetxt.encode('utf-8')).hexdigest()

        try:
            import requests
            session = requests.Session()
            response = session.get(url, timeout=timeout)
            response.raise_for_status()
            for x in response.history:
                x.raise_for_status()
                if (x.status_code == 500 or x.status_code == 307 or
                        x.status_code == 404):
                    raise
            txt = response.text
            newmd5 = md5(txt.encode('utf-8')).hexdigest()
            # tprint(filemd5 + ": " + newmd5)
            if self.args.update and newmd5 == filemd5:
                self.log.debug(
                    'Skipping file in "' + self.current_task +
                    '," local and remote copies identical [' + newmd5 + '].')
                return False
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            if failhard:
                return ''
            return filetxt
        else:
            if write:
                with codecs.open(filepath, 'w', encoding='utf8') as f:
                    f.write(txt if txt else filetxt)
        return txt


def _get_task_priority(tasks, task_priority):
    """Get the task `priority` corresponding to the given `task_priority`.

    If `task_priority` is an integer or 'None', return it.
    If `task_priority` is a str, return the priority of the task it matches.
    Otherwise, raise `ValueError`.
    """
    if task_priority is None:
        return None
    if is_integer(task_priority):
        return task_priority
    if isinstance(task_priority, str):
        if task_priority in tasks:
            return tasks[task_priority].priority

    raise ValueError("Unrecognized task priority '{}'".format(task_priority))
