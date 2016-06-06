"""Open Supernova Catalog (OCS) Scripts for Downloading and Processing Supernova data.
"""

print("scripts/__init__.py")

import os
# from enum import Enum

_PATH_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


class PATH:
    ROOT = str(_PATH_ROOT)
    DATA = os.path.join(ROOT, 'scripts/aux_data')
    REPO_BONEYARD = os.path.join(ROOT, 'sne-boneyard')
    REPO_INTERNAL = os.path.join(ROOT, 'sne-internal')
    REPO_EXTERNAL_RADIO = os.path.join(ROOT, 'sne-external-radio')
    REPO_EXTERNAL_XRAY = os.path.join(ROOT, 'sne-external-xray')


class FILENAME:
    TYPE_SYNONYMS = os.path.join(PATH.DATA, 'type-synonyms.json')
    SOURCE_SYNONYMS = os.path.join(PATH.DATA, 'source-synonyms.json')
    NON_SNE_TYPES = os.path.join(PATH.DATA, 'non-sne-types.json')
    REPOS_LIST = os.path.join(PATH.DATA, 'rep-folders.txt')
    ATELS = os.path.join(PATH.DATA, 'atels.json')
    CBETS = os.path.join(PATH.DATA, 'cbets.json')
    IAUCS = os.path.join(PATH.DATA, 'iaucs.json')
    EXTINCT = os.path.join(PATH.DATA, 'extinctions.json')
    BIBAUTHORS = os.path.join(PATH.DATA, 'bibauthors.json')

# _PATH_DATA = os.path.join(_PATH_ROOT, 'scripts/aux_data')
# _PATH_REPO_BONEYARD = os.path.join(_PATH_ROOT, 'sne-boneyard')
# _PATH_REPO_INTERNAL = os.path.join(_PATH_ROOT, 'sne-internal')

# _FILENAME_TYPE_SYNONYMS = os.path.join(_PATH_DATA, 'type-synonyms.json')
# _FILENAME_SOURCE_SYNONYMS = os.path.join(_PATH_DATA, 'source-synonyms.json')
# _FILENAME_NON_SNE_TYPES = os.path.join(_PATH_DATA, 'non-sne-types.json')
# _FILENAME_REPOS = os.path.join(_PATH_DATA, 'rep-folders.txt')
# _FILENAME_ATELS = os.path.join(_PATH_DATA, 'atels.json')
# _FILENAME_CBETS = os.path.join(_PATH_DATA, 'cbets.json')
# _FILENAME_IAUCS = os.path.join(_PATH_DATA, 'iaucs.json')
# _FILENAME_EXTINCT = os.path.join(_PATH_DATA, 'extinctions.json')
# _FILENAME_BIBAUTHORS = os.path.join(_PATH_DATA, 'bibauthors.json')
