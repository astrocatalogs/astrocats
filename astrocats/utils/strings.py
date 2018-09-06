"""
"""
import json
import urllib

from .digits import round_sig

__all__ = [
    'decode_url', 'dict_to_pretty_string', 'rep_chars', 'get_filename',
    'single_spaces', 'trim_str_arr', 'uniq_cdl', 'utf8'
]


def rep_chars(string, chars, rep=''):
    for c in chars:
        if c in string:
            string = string.replace(c, rep)
    return string


def single_spaces(string):
    return ' '.join(list(filter(None, string.split())))


def trim_str_arr(arr, length=10, max_rows=10):
    do_full = False
    for i, x in enumerate(arr):
        # Check the first max_rows rows, if no changes needed don't bother
        # converting values.
        if not do_full and i >= max_rows:
            return arr
        lenx = len(x)
        if lenx <= length:
            continue
        round_str = str(round_sig(float(x), length))
        if len(round_str) < lenx:
            arr[i] = round_str
            do_full = True
    return arr


def uniq_cdl(values):
    return ','.join(sorted(list(set(values))))


def utf8(x):
    return str(x, 'utf-8')


def dict_to_pretty_string(odict):
    jsonstring = json.dumps(odict, indent=4, separators=(',', ':'), ensure_ascii=False)
    return jsonstring


def get_filename(name):
    return name.replace('/', '_')


def decode_url(url):
    decode = urllib.parse.unquote(url)
    return decode
