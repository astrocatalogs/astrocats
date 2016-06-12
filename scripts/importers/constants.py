"""Constant variables used by OSC import methods.
"""
from astropy import constants as const
from astropy import units as un

OSC_BIBCODE = '2016arXiv160501054G'
OSC_NAME = 'The Open Supernova Catalog'
OSC_URL = 'https://sne.space'

ACKN_CFA = ("This research has made use of the CfA Supernova Archive, "
            "which is funded in part by the National Science Foundation "
            "through grant AST 0907903.")

CLIGHT = const.c.cgs.value
KM = (1.0 * un.km).cgs.value
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
    ['B', 'b', 'g'],  # B-like bands first
    ['V', 'G'],       # if not, V-like bands
    ['R', 'r']        # if not, R-like bands
]

PREF_KINDS = ['heliocentric', 'cmb', 'spectroscopic', 'photometric', 'host', 'cluster', '']


class TASK():
    name = ''
    nice_name = ''   # Name for pretty printing
    update = False   # Perform task during update
    archived = None  # Use archived data ???
    active = True    # Whether this task should be performed or not

    module = None    # Module in which to find `function` for carrying out this task
    function = ''    # Function to execute when carrying out this task
    priority = None  # Order in which tasks should be executed

    def __init__(self, **kwargs):
        for key, val in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, val)
            else:
                raise ValueError("No attribute '{}'".format(key))

    def __repr__(self):
        retval = ("TASK(name='{}', nice_name='{}', update='{}', archived='{}', "
                  "module='{}', function='{}', priority='{}'")
        retval = retval.format(self.name, self.nice_name, self.update, self.archived,
                               self.module, self.function, self.priority)
        return retval
