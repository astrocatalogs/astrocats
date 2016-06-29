"""Supernovae specific Constant variables.
"""
import os
from astropy import constants as const
from astropy import units as un

_PATH_SUPERNOVAE = os.path.abspath(os.path.dirname(__file__))
_PATH_INPUT = os.path.join(_PATH_SUPERNOVAE, 'input', '')
_PATH_OUTPUT = os.path.join(_PATH_SUPERNOVAE, 'output', '')


class FILENAME:
    # critical datafiles
    REPOS = os.path.join(_PATH_INPUT, 'repos.json')
    TASK_LIST = os.path.join(_PATH_INPUT, 'tasks.json')
    # auxiliary datafiles
    TYPE_SYNONYMS = os.path.join(_PATH_INPUT, 'type-synonyms.json')
    SOURCE_SYNONYMS = os.path.join(_PATH_INPUT, 'source-synonyms.json')
    NON_SNE_TYPES = os.path.join(_PATH_INPUT, 'non-sne-types.json')
    NON_SNE_PREFIXES = os.path.join(_PATH_INPUT, 'non-sne-prefixes.json')
    BIBERRORS = os.path.join(_PATH_INPUT, 'biberrors.json')

    BIBAUTHORS = os.path.join(_PATH_OUTPUT, 'cache', 'bibauthors.json')
    EXTINCT = os.path.join(_PATH_OUTPUT, 'cache', 'extinctions.json')


CLIGHT = const.c.cgs.value
KM = (1.0 * un.km).cgs.value

PREF_KINDS = ['heliocentric', 'cmb', 'spectroscopic',
              'photometric', 'host', 'cluster', '']

REPR_BETTER_QUANTITY = {
    'redshift',
    'ebv',
    'velocity',
    'lumdist',
    'discoverdate',
    'maxdate'
}

MAX_BANDS = [
    ['B', 'b', 'g'],  # B-like bands first
    ['V', 'G'],       # if not, V-like bands
    ['R', 'r']        # if not, R-like bands
]
