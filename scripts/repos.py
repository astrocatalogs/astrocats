import glob
import sys
from digits import *

with open('rep-folders.txt', 'r') as f:
    repofolders = f.read().splitlines()

repoyears = [int(repofolders[x][-4:]) for x in range(len(repofolders)-1)]
repoyears[0] -= 1

def repo_file_list(bones = True):
    files = []
    for rep in repofolders:
        if not bones and 'boneyard' in rep:
            continue
        files += glob.glob('../' + rep + "/*.json") + glob.glob('../' + rep + "/*.json.gz")
    return files

def get_rep_folder(entry):
    if 'discoverdate' not in entry:
        return repofolders[0]
    if not is_number(entry['discoverdate'][0]['value'].split('/')[0]):
        raise(ValueError('Discovery year is not a number!'))
        sys.exit()
    for r, repoyear in enumerate(repoyears):
        if int(entry['discoverdate'][0]['value'].split('/')[0]) <= repoyear:
            return repofolders[r]
    return repofolders[0]
