"""Open Supernova Catalog (OCS) Scripts for Downloading and Processing Supernova data.
"""

import os
from subprocess import check_output

_PATH_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


class PATH:
    ROOT = str(_PATH_ROOT)
    DATA = os.path.join(ROOT, 'scripts/aux_data')
    IMPORT = os.path.join(ROOT, 'scripts/importer')
    REPO_BONEYARD = os.path.join(ROOT, 'sne-boneyard')
    REPO_INTERNAL = os.path.join(ROOT, 'sne-internal')
    REPO_EXTERNAL = os.path.join(ROOT, 'sne-external')
    REPO_EXTERNAL_RADIO = os.path.join(ROOT, 'sne-external-radio')
    REPO_EXTERNAL_XRAY = os.path.join(ROOT, 'sne-external-xray')
    REPO_EXTERNAL_SPECTRA = os.path.join(ROOT, 'sne-external-spectra')
    REPO_EXTERNAL_WISEREP = os.path.join(ROOT, 'sne-external-WISEREP')


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
    TASK_LIST = os.path.join(PATH.IMPORT, 'tasks.json')


class SCHEMA:
    HASH = check_output(['git', 'log', '-n', '1', '--format="%H"',
        '--', 'OSC-JSON-format.md']).decode('ascii').strip().strip('"').strip()
    URL = 'https://github.com/astrocatalogs/sne/blob/' + HASH + '/OSC-JSON-format.md'
