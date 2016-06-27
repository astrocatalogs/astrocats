from collections import OrderedDict

from .funcs import get_bibauthor_dict, get_extinctions_dict


class Catalog():
    """
    Object to hold the main catalog dictionary and other catalog globals.
    """

    biberrordict = {
        "2012Sci..337..942D": "2012Sci...337..942D",
        "2012MNRAS.420.1135": "2012MNRAS.420.1135S",
        "2014MNRAS.438,368": "2014MNRAS.438..368T",
        "2006ApJ...636...400Q": "2006ApJ...636..400Q",
        "0609268": "2007AJ....133...58K",
        "2004MNRAS.tmp..131P": "2004MNRAS.352..457P",
        "2013MNRAS.tmp.1499F": "2013MNRAS.433.1312F",
        "1991MNRAS.247P.410B": "1991A&A...247..410B",
        "2011Sci.333..856S": "2011Sci...333..856S"
    }

    nedd_dict = OrderedDict()
    bibauthor_dict = get_bibauthor_dict()
    extinctions_dict = get_extinctions_dict()

    def __init__(self):
        return