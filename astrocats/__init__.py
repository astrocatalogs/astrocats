"""Astrocats.

Scripts for creating and analyzing catalogs of astronomical data.
"""

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import os
import glob

_PATH_VERSION = os.path.join(os.path.dirname(__file__), os.path.pardir, "VERSION")
with open(_PATH_VERSION, "r") as inn:
    version = inn.read().strip()

__version__ = version
__author__ = 'James Guillochon & Luke Zoltan Kelley'
__license__ = 'MIT'


class Paths(object):
    """
    """

    ASTROCATS = os.path.join(os.path.dirname(__file__), "")

    # These need to be overridden in subclasses
    ROOT = os.path.join(os.path.dirname(__file__), "")
    NAME = __name__
    FILE = __file__

    def __init__(self):
        # Determine if this instance is from a derived class (as apposed to this class itself)
        is_astrocats = (type(self) == Paths)
        is_derived = (not is_astrocats)

        # Make sure these attributes have been set (use strings for error messages)
        check_not_none = []

        if self.ASTROCATS != os.path.join(os.path.dirname(__file__), ""):
            raise ValueError("`ASTROCATS` attribute must not be overridden!")

        # Make sure that derived classes look okay
        if is_derived:
            if os.path.realpath(self.ASTROCATS).lower() == os.path.realpath(self.ROOT).lower():
                raise ValueError("`ROOT` attribute must be overridden!")

            # Make sure required parameters have been changed to legit values (instead of `None`)
            for check in check_not_none:
                check_val = getattr(self, check, None)
                if (check_val is None):
                    err = "Attribute `{}` is 'None' but must be overridden!".format(check)
                    raise ValueError(err)

        check_dirs_astrocats = []
        check_dirs_derived = []

        required_files_derived = []

        self.STRUCTURES = os.path.join(self.ROOT, "structures", "")

        # Schema files
        # -----------------------------
        self.SCHEMA = os.path.join(self.STRUCTURES, "schema", "")
        self.SCHEMA_INPUT = os.path.join(self.SCHEMA, "input", "")
        self.SCHEMA_OUTPUT = os.path.join(self.SCHEMA, "output", "")
        # check_dirs_derived.extend([self.SCHEMA_OUTPUT])
        check_dirs_astrocats.extend([self.SCHEMA_OUTPUT])

        self.SCHEMA_INPUT_ASTROSCHEMA = os.path.join(self.SCHEMA_INPUT, "astroschema", "")
        self.SCHEMA_INPUT_ASTROCATS = os.path.join(self.SCHEMA_INPUT, "astrocats", "")
        check_dirs_astrocats.extend([self.SCHEMA_INPUT_ASTROSCHEMA, self.SCHEMA_INPUT_ASTROCATS])

        self.INPUT = os.path.join(self.ROOT, 'input', '')
        # check_dirs.append(self.INPUT)

        # Output Paths and Files
        # ----------------------
        self.OUTPUT = os.path.join(self.ROOT, 'output', '')
        check_dirs_derived.append(self.OUTPUT)

        self.BIBLIO_FILE = os.path.join(self.OUTPUT, 'biblio.json')
        self.BIBLIO_MIN_FILE = os.path.join(self.OUTPUT, 'biblio.min.json')

        self.WEB_TABLE_FILE = os.path.join(self.OUTPUT, 'catalog.json')
        self.WEB_TABLE_MIN_FILE = os.path.join(self.OUTPUT, 'catalog.min.json')

        self.NAMES_FILE = os.path.join(self.OUTPUT, 'names.json')
        self.NAMES_MIN_FILE = os.path.join(self.OUTPUT, 'names.min.json')

        # Cache path and files
        self.CACHE = os.path.join(self.OUTPUT, 'cache', '')
        check_dirs_derived.append(self.CACHE)

        self.MD5_FILE = os.path.join(self.CACHE, 'md5s.json')
        self.AUTHORS_FILE = os.path.join(self.CACHE, 'bibauthors.json')
        self.ALL_AUTHORS_FILE = os.path.join(self.CACHE, 'biballauthors.json')
        self.HOST_IMAGES_FILE = os.path.join(self.CACHE, 'host_images.json')

        # json path and files
        self.JSON = os.path.join(self.OUTPUT, 'json', '')
        check_dirs_derived.append(self.JSON)
        # html path and files
        self.HTML = os.path.join(self.OUTPUT, 'html', '')
        check_dirs_derived.append(self.HTML)

        # critical datafiles
        # ------------------
        self.TASKS = os.path.join(self.ROOT, 'tasks', "")
        self.REPOS_FILE = os.path.join(self.INPUT, 'repos.json')
        self.TASKS_FILE = os.path.join(self.INPUT, 'tasks.json')
        required_files_derived.extend([self.REPOS_FILE, self.TASKS_FILE])

        self._derived = is_derived

        if is_astrocats:
            # Make sure astrocats required directories exist
            for cd in check_dirs_astrocats:
                self._check_create_dir(cd)

        if is_derived:
            # Make sure derived required directories exist
            for cd in check_dirs_derived:
                self._check_create_dir(cd)

            # Make sure derived required files exist
            for fname in required_files_derived:
                if not os.path.exists(fname):
                    raise RuntimeError("Required file '{}' does not exist!".format(fname))

                if not os.path.isfile(fname):
                    raise RuntimeError("Required file '{}' is not a file!".format(fname))

        # self.catalog = catalog
        # this_file = os.path.abspath(sys.modules[self.__module__].__file__)
        # Path of the `catalog`
        # self.catalog_dir = os.path.dirname(this_file)
        # Path in which `catalog` resides
        # self.root_dir = os.path.realpath(os.path.join(self.catalog_dir, os.path.pardir))
        # self.tasks_dir = os.path.join(self.catalog_dir, 'tasks')

        return

    def _check_create_dir(self, direct):
        if not os.path.isdir(direct):
            print("Creating path '{}'".format(direct))
            os.makedirs(direct)

        return

    def __str__(self):
        rstr = "`ROOT` = '{}'".format(self.ROOT)
        rstr += "\n`INPUT` = '{}'".format(self.INPUT)
        rstr += "\n`OUTPUT` = '{}'".format(self.OUTPUT)
        rstr += "\n\t`BIBLIO_FILE` = '{}'".format(self.BIBLIO_FILE)
        rstr += "\n\t`NAMES_MIN_FILE` = '{}'".format(self.NAMES_MIN_FILE)
        rstr += "\n\t`WEB_TABLE_FILE` = '{}'".format(self.WEB_TABLE_FILE)
        rstr += "\n\t`CACHE` = '{}'".format(self.CACHE)
        rstr += "\n\t\t`MD5_FILE` = '{}'".format(self.MD5_FILE)
        rstr += "\n\t\t`AUTHORS_FILE` = '{}'".format(self.AUTHORS_FILE)
        rstr += "\n\t\t`ALL_AUTHORS_FILE` = '{}'".format(self.ALL_AUTHORS_FILE)
        rstr += "\n\t`JSON` = '{}'".format(self.JSON)
        rstr += "\n\t`HTML` = '{}'".format(self.HTML)
        rstr += "\n`REPOS_FILE` = '{}'".format(self.REPOS_FILE)
        rstr += "\n`TASKS_FILE` = '{}'".format(self.TASKS_FILE)
        return rstr

    def _get_repo_file_list(self, repo_folders, normal=True, bones=True):
        """Get filenames for files in each repository.

        `boneyard` optionally include with `bones=True`.
        """
        # repo_folders = get_repo_output_folders()
        files = []
        for rep in repo_folders:
            if 'boneyard' not in rep and not normal:
                continue
            if not bones and 'boneyard' in rep:
                continue
            these_files = glob.glob(rep + "/*.json") + glob.glob(rep + "/*.json.gz")
            # self.log.debug("Found {} files in '{}'".format(len(these_files), rep))
            files += these_files

        return files

    @property
    def repos_dict(self):
        if hasattr(self, '_repos_dict'):
            return self._repos_dict
        from astrocats import utils
        repos = utils.read_json_dict(self.REPOS_FILE)
        self._repos_dict = repos
        return repos

    def get_all_repo_folders(self, boneyard=True, private=False):
        """Get the full paths of all data repositories."""
        all_repos = self.get_repo_input_folders(private=private)
        all_repos.extend(self.get_repo_output_folders(bones=boneyard))
        return all_repos

    def get_repo_boneyard(self):
        bone_path = self.repos_dict.get('boneyard', [])
        try:
            bone_path = bone_path[0]
        except TypeError:
            pass
        bone_path = os.path.join(self.OUTPUT, bone_path, '')
        return bone_path

    def get_filename_for_internal_event(self, event_name, suffix='.json'):
        """Get the full path to the target input folder.
        """
        folder = self.repos_dict.get('internal')
        # 'internal' must be a unique folder
        if (folder is None) or (len(folder) != 1):
            return None

        fname = os.path.join(self.INPUT, folder[0], event_name + suffix)
        return fname

    def get_repo_input_folders(self, private=False):
        """Get the full paths of the input data repositories."""
        repo_folders = []
        repo_folders += self.repos_dict.get('external', [])
        repo_folders += self.repos_dict.get('internal', [])
        if private:
            repo_folders += self.repos_dict.get('private', [])
        repo_folders = list(sorted(set(repo_folders)))
        repo_folders = [
            os.path.join(self.INPUT, rf) for rf in repo_folders
            if len(rf)
        ]
        return repo_folders

    def get_md5_filename(self, fname='md5s.json'):
        return os.path.join(self.CACHE, fname)

    def get_repo_output_file_list(self, normal=True, bones=True):
        """Get a list of all existing output files.

        These are the files deleted in the `delete_old_entry_files` task.
        """
        repo_folders = self.get_repo_output_folders(bones=bones)
        repo_files = self._get_repo_file_list(repo_folders, normal=normal, bones=bones)
        return repo_files

    def get_repo_output_folders(self, bones=True):
        """Get the full paths of the output data repositories."""
        repo_folders = []
        repo_folders += self.repos_dict.get('output', [])
        if bones:
            repo_folders += self.repos_dict.get('boneyard', [])
        # repo_folders = list(
        #     sorted(
        #         list(set(repo_folders)),
        #         key=lambda key: utils.repo_priority(key)))
        repo_folders = [
            os.path.join(self.OUTPUT, rf) for rf in repo_folders
            if len(rf)
        ]
        return repo_folders

    def is_internal_event(self, event_name):
        """Check if the given event corresponds to an 'internal' file.
        """
        fname = self.get_filename_for_internal_event(event_name)
        if os.path.isfile(fname):
            return True
        return False


PATHS = Paths()


# =============================================

'''
warnings.warn("Adding `pyastroschema` to `sys.path` for easy access... fix this")
_PAS_PATH = "/Users/lzkelley/Research/catalogs/astroschema"
if _PAS_PATH not in sys.path:
    sys.path.append(_PAS_PATH)
'''

# =====================================

'''
_RM_PATH = "/Users/lzkelley/Research/catalogs/redesign/astrocats/astrocats/catalog/schema/output/"
warnings.warn("Removing schema files from {} for testing!".format(_RM_PATH))
pattern = os.path.join(_RM_PATH, "", "*.json")
old_files = sorted(glob.glob(pattern))
old_files = [os.path.join(_RM_PATH, fname) for fname in old_files]
for fname in old_files:
    os.remove(fname)
'''

# ======================================

from .utils import gitter

__git_version__ = gitter.get_git_sha(os.path.dirname(__file__))
