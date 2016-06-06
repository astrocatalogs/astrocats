"""Open Supernova Catalog (OCS) Scripts for Downloading and Processing Supernova data.
"""

print("scripts/__init__.py")

import os

_PATH_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
_PATH_DATA = os.path.join(_PATH_ROOT, 'scripts/aux_data')

_FILENAME_TYPE_SYNONYMS = os.path.join(_PATH_DATA, 'type-synonyms.json')
_FILENAME_SOURCE_SYNONYMS = os.path.join(_PATH_DATA, 'source-synonyms.json')
_FILENAME_NON_SNE_TYPES = os.path.join(_PATH_DATA, 'non-sne-types.json')
_FILENAME_REPOS = os.path.join(_PATH_DATA, 'rep-folders.txt')
_FILENAME_ATELS = os.path.join(_PATH_DATA, 'atels.json')
_FILENAME_CBETS = os.path.join(_PATH_DATA, 'cbets.json')
_FILENAME_IAUCS = os.path.join(_PATH_DATA, 'iaucs.json')
