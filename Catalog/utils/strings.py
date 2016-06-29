"""
"""

__all__ = ['rep_chars', 'single_spaces', 'get_event_filename']


def rep_chars(string, chars, rep=''):
    for c in chars:
        if c in string:
            string = string.replace(c, rep)
    return string


def single_spaces(string):
    return ' '.join(list(filter(None, string.split())))


def get_event_filename(name):
    return name.replace('/', '_')


def uniq_cdl(values):
    return ','.join(sorted(list(set(values))))


def trim_str_arr(arr, length=10):
    return [str(round_sig(float(x), length)) if
            (len(x) > length and
             len(str(round_sig(float(x), length))) < len(x))
            else x for x in arr]


def utf8(x):
    return str(x, 'utf-8')
