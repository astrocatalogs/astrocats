"""
"""
import json
from .digits import round_sig

__all__ = ['dict_to_pretty_string', 'get_event_filename', 'rep_chars',
           'single_spaces', 'trim_str_arr', 'uniq_cdl', 'utf8']


def get_event_filename(name):
    return name.replace('/', '_')


def rep_chars(string, chars, rep=''):
    for c in chars:
        if c in string:
            string = string.replace(c, rep)
    return string


def single_spaces(string):
    return ' '.join(list(filter(None, string.split())))


def trim_str_arr(arr, length=10):
    return [str(round_sig(float(x), length)) if
            (len(x) > length and
             len(str(round_sig(float(x), length))) < len(x))
            else x for x in arr]


def uniq_cdl(values):
    return ','.join(sorted(list(set(values))))


def utf8(x):
    return str(x, 'utf-8')


def dict_to_pretty_string(odict):
    jsonstring = json.dumps(odict,
                            indent='\t', separators=(',', ':'),
                            ensure_ascii=False)
    return jsonstring
