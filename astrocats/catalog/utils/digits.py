"""
"""

from math import floor, log10

from past.builtins import basestring

__all__ = [
    'get_sig_digits', 'is_integer', 'is_number', 'pretty_num', 'round_sig',
    'zpad'
]


def get_sig_digits(x):
    return len((''.join(x.split('.'))).strip('0'))


def is_integer(s):
    if isinstance(s, list) and not isinstance(s, basestring):
        try:
            [int(x) for x in s]
            return True
        except ValueError:
            return False
    else:
        try:
            int(s)
            return True
        except ValueError:
            return False


def is_number(s):
    if isinstance(s, list) and not isinstance(s, basestring):
        try:
            for x in s:
                if isinstance(x, basestring) and ' ' in x:
                    raise ValueError
            [float(x) for x in s]
            return True
        except ValueError:
            return False
    else:
        try:
            if isinstance(s, basestring) and ' ' in s:
                raise ValueError
            float(s)
            return True
        except ValueError:
            return False


def pretty_num(x, sig=4):
    return str('%g' % (round_sig(x, sig)))


def round_sig(x, sig=4):
    if x == 0.0:
        return 0.0
    return round(x, sig - int(floor(log10(abs(x)))) - 1)


def zpad(val, n=2):
    bits = val.split('.')
    if len(bits) != 2:
        return val.zfill(n)
    return "%s.%s" % (bits[0].zfill(n), bits[1])
