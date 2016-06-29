import os
import warnings
from glob import glob

# from astrocats import FILENAME, PATH

from .digits import is_number

__all__ = []
# __all__ = ['get_repo_file_list', 'get_repo_folder_for_year', 'get_output_repo_folders',
#            'get_repo_paths', 'get_repo_years']


def get_repo_file_list(normal=True, bones=True):
    """Get filenames for all files in each repository, with `boneyard` files
    optional.
    """
    repo_folders = get_output_repo_folders()
    files = []
    for rep in repo_folders:
        rep_path = os.path.join(PATH.ROOT, rep)
        if 'boneyard' not in rep and not normal:
            continue
        if not bones and 'boneyard' in rep:
            continue
        files += glob(rep_path + "/*.json") + glob(rep_path + "/*.json.gz")

    return files


def get_repo_folder_for_year(entry):
    """Determine the appropriate SN repository based on the `discoverdate` year
    of this `entry`.
    """
    repo_folders = get_output_repo_folders()
    if 'discoverdate' not in entry:
        return repo_folders[0]
    if not is_number(entry['discoverdate'][0]['value'].split('/')[0]):
        warnings.warn('Discovery year is not a number')
        return repo_folders[0]

    repo_years = get_repo_years(repo_folders)
    for r, repoyear in enumerate(repo_years):
        if int(entry['discoverdate'][0]['value'].split('/')[0]) <= repoyear:
            return repo_folders[r]
    return repo_folders[0]


def get_output_repo_folders():
    """Get the names of all repositories given in the 'rep-folders.txt' file.
    """
    with open(FILENAME.REPOS_LIST, 'r') as f:
        repo_folders = f.read().splitlines()
    return repo_folders


def get_repo_paths():
    repo_folders = get_output_repo_folders()
    repo_paths = [os.path.join(PATH.ROOT, rf, '') for rf in repo_folders]
    return repo_paths


def get_repo_years(repo_folders):
    """Get the years section of all repository names given in the
    'rep-folders.txt' file.
    """
    repo_years = [int(repo_folders[x][-4:])
                  for x in range(len(repo_folders) - 1)]
    repo_years[0] -= 1
    return repo_years
