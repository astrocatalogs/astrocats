"""Constant variables used by OSC import methods.
"""

OSC_BIBCODE = '2016arXiv160501054G'
OSC_NAME = 'The Open Supernova Catalog'
OSC_URL = 'https://sne.space'

ACKN_CFA = ("This research has made use of the CfA Supernova Archive, "
          "which is funded in part by the National Science Foundation "
          "through grant AST 0907903.")

CLIGHT = const.c.cgs.value
KM = (1.0 * un.KM).cgs.value
TRAVIS_QUERY_LIMIT = 10


REPR_BETTER_QUANTITY = {
    'redshift',
    'ebv',
    'velocity',
    'lumdist',
    'discoverdate',
    'maxdate'
}

MAX_BANDS = [
    ['B', 'b', 'g'], # B-like bands first
    ['V', 'G'],      # if not, V-like bands
    ['R', 'r']       # if not, R-like bands
]


PREF_KINDS = ['heliocentric', 'cmb', 'spectroscopic', 'photometric', 'host', 'cluster', '']

