"""Open Supernova Catalog (OCS) Scripts for Downloading and Processing Supernova data.
"""

print("scripts/__init__.py")

import os

_PATH_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
_PATH_DATA = os.path.join(_PATH_ROOT, 'aux_data')

_FILENAME_TYPE_SYNONYMS = os.path.join(_PATH_DATA, 'type-synonyms.json')
_FILENAME_SOURCE_SYNONYMS = os.path.join(_PATH_DATA, 'source-synonyms.json')
_FILENAME_NON_SNE_TYPES = os.path.join(_PATH_DATA, 'non-sne-types.json')
