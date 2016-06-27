from collections import OrderedDict

from .importer.funcs import (get_bibauthor_dict, get_biberror_dict,
                             get_extinctions_dict)


class Catalog():
    """
    Object to hold the main catalog dictionary and other catalog globals.
    """

    nedd_dict = OrderedDict()
    bibauthor_dict = get_bibauthor_dict()
    biberror_dict = get_biberror_dict()
    extinctions_dict = get_extinctions_dict()

    def __init__(self):
        return
