import warnings
from glob import glob

from astrocats.catalog.utils import is_number



def get_rep_folders(module):
    with open('astrocats/' + module + '/input/rep-folders.txt', 'r') as f:
        return f.read().splitlines()

def get_rep_years(repofolders):
    repoyears = [int(repofolders[x][-4:]) for x in range(len(repofolders) - 1)]
    repoyears[0] -= 1
    return repoyears

def repo_file_list(module, repofolders, normal=True, bones=True):
    outdir = 'astrocats/' + module + '/output/'
    files = []
    for rep in repofolders:
        if 'boneyard' not in rep and not normal:
            continue
        if not bones and 'boneyard' in rep:
            continue
        files += glob(outdir + rep + "/*.json") + \
            glob(outdir + rep + "/*.json.gz")
    return files

def get_rep_folder(entry, repofolders):
    if 'discoverdate' not in entry:
        return repofolders[0]
    if not is_number(entry['discoverdate'][0]['value'].split('/')[0]):
        warnings.warn('Discovery year is not a number')
        return repofolders[0]
    repoyears = get_rep_years(repofolders)
    for r, repoyear in enumerate(repoyears):
        if int(entry['discoverdate'][0]['value'].split('/')[0]) <= repoyear:
            return repofolders[r]
    return repofolders[0]
