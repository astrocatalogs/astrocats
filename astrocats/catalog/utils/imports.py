'''Utility functions related to importing data.
'''
import json
import os
from collections import OrderedDict

from .digits import is_number

__all__ = ['compress_gz', 'convert_aq_output', 'read_json_dict',
           'read_json_arr', 'uncompress_gz']


def convert_aq_output(row):
    return OrderedDict([(x, str(row[x]) if is_number(row[x]) else row[x])
                        for x in row.colnames])


def read_json_dict(filename):
    # path = '../atels.json'
    if os.path.isfile(filename):
        with open(filename, 'r') as f:
            mydict = json.loads(f.read(), object_pairs_hook=OrderedDict)
    else:
        mydict = OrderedDict()
    return mydict


def read_json_arr(filename):
    if os.path.isfile(filename):
        with open(filename, 'r') as f:
            myarr = json.loads(f.read())
    else:
        myarr = []
    return myarr


def compress_gz(fname):
    """Compress the file with the given name and delete the uncompressed file.

    The compressed filename is simply the input filename with '.gz' appended.

    Arguments
    ---------
    fname : str
        Name of the file to compress and delete.

    Returns
    -------
    comp_fname : str
        Name of the compressed file produced.  Equal to `fname + '.gz'`.
    """
    import shutil
    import gzip
    comp_fname = fname + '.gz'
    with open(fname, 'rb') as f_in, gzip.open(comp_fname, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(fname)
    return comp_fname


def uncompress_gz(fname):
    import shutil
    import gzip
    uncomp_name = fname.replace('.gz', '')
    with gzip.open(fname, 'rb') as f_in, open(uncomp_name, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(fname)
    return uncomp_name
