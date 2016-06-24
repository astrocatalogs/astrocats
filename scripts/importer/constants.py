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

COMPRESS_ABOVE_FILESIZE = 90000000   # FIX: units?

MAX_BANDS = [
    ['B', 'b', 'g'],  # B-like bands first
    ['V', 'G'],       # if not, V-like bands
    ['R', 'r']        # if not, R-like bands
]

PREF_KINDS = ['heliocentric', 'cmb', 'spectroscopic', 'photometric', 'host', 'cluster', '']

NON_SNE_PREFIXES = ['PNVJ', 'PNV J', 'OGLE-2013-NOVA', 'EV*', 'V*', 'Nova']


class TASK:
    name = ''
    nice_name = ''   # Name for pretty printing
    update = False   # Perform task during update
    archived = False  # Use archived data ???
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
        retval = ("TASK(name='{}', nice_name='{}', active='{}', update='{}', archived='{}', "
                  "module='{}', function='{}', priority='{}'")
        retval = retval.format(self.name, self.nice_name, self.active, self.update, self.archived,
                               self.module, self.function, self.priority)
        return retval

    def current_task(self, args):
        """Name of current action for progress-bar output, depends on run configuration.
        """
        ctask = self.nice_name if self.nice_name else self.name
        if args is not None:
            if args.update:
                ctask = ctask.replace('%pre', 'Updating')
            else:
                ctask = ctask.replace('%pre', 'Loading')
        return ctask

    def load_archive(self, args):
        """Depending on run configuration, whether previously archived data should be loaded.
        """
        if not self.archived:
            return False
        # If we're running in 'archived' mode, and only loading 'archived' things, then True
        if args.archived not in args.refresh_list and not args.full_refresh:
            return True
        # For normal running, if we are not sepcifically refreshing this task, then True
        if self.name not in args.refresh_list and not args.full_refresh:
            return True

        return False
