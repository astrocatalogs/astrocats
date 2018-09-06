"""Utility functions related to importing data."""
import codecs
import json
import os
from collections import OrderedDict

from .digits import is_number

__all__ = ['ADD_FAIL_ACTION', 'compress_gz', 'convert_aq_output', 'read_json_dict',
           'read_json_arr', 'uncompress_gz', 'import_ads']


class ADD_FAIL_ACTION:
    IGNORE = "ignore"
    WARN = "warn"
    RAISE = "raise"


def convert_aq_output(row):
    return OrderedDict([(x, str(row[x]) if is_number(row[x]) else row[x]) for x in row.colnames])


def read_json_dict(filename):
    # path = '../atels.json'
    if os.path.isfile(filename):
        with codecs.open(filename, 'r') as f:
            mydict = json.loads(f.read(), object_pairs_hook=OrderedDict)
    else:
        mydict = OrderedDict()
    return mydict


def read_json_arr(filename):
    if os.path.isfile(filename):
        with codecs.open(filename, 'r') as f:
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
    with codecs.open(fname, 'rb') as f_in, gzip.open(
            comp_fname, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(fname)
    return comp_fname


def uncompress_gz(fname):
    import shutil
    import gzip
    uncomp_name = fname.replace('.gz', '')
    with gzip.open(fname, 'rb') as f_in, codecs.open(
            uncomp_name, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(fname)
    return uncomp_name


def import_ads():
    """Load and return the ads package checking that a token file exists.
    """
    import ads

    ads_key_path = os.path.join(os.path.expanduser('~'), '.ads/dev_key')
    if not os.path.exists(ads_key_path):
        local_path = 'ads.key'
        if os.path.isfile(local_path):
            with open(local_path, 'r') as ff:
                ads.config.token = ff.read().splitlines()[0]
        else:
            token_url = "https://ui.adsabs.harvard.edu/#user/settings/token"
            err = "Cannot find '{}' or '{}'.".format(ads_key_path, local_path)
            err += "Generate one at '{}', and place it in one of these files.".format(token_url)
            raise IOError(err)

    return ads
