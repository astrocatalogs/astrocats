import sys
from glob import glob
from . digits import is_number

__all__ = ['repo_file_list', 'get_rep_folder']


def repo_file_list(bones=True):
    """
    """
    repo_folders = _get_repo_folders()
    files = []
    for rep in repo_folders:
        if not bones and 'boneyard' in rep:
            continue
        files += glob('../' + rep + "/*.json") + glob('../' + rep + "/*.json.gz")
    return files


def get_rep_folder(entry):
    repo_folders = _get_repo_folders()
    if 'discoverdate' not in entry:
        return repo_folders[0]
    if not is_number(entry['discoverdate'][0]['value'].split('/')[0]):
        raise(ValueError('Discovery year is not a number!'))
        sys.exit()

    repo_years = _get_repo_years(repo_folders)
    for r, repoyear in enumerate(repo_years):
        if int(entry['discoverdate'][0]['value'].split('/')[0]) <= repoyear:
            return repo_folders[r]
    return repo_folders[0]


def _get_repo_folders():
    """Get the names of all repositories given in the 'rep-folders.txt' file.
    """
    _REPO_FILENAME = '../rep-folders.txt'
    with open(_REPO_FILENAME, 'r') as f:
        repo_folders = f.read().splitlines()
    return repo_folders


def _get_repo_years(repo_folders):
    """Get the years section of all repository names given in the 'rep-folders.txt' file.
    """
    repo_years = [int(repo_folders[x][-4:]) for x in range(len(repo_folders)-1)]
    repo_years[0] -= 1
    return repo_years
