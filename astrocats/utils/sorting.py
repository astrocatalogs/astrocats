'''Key sorting functions
'''

from collections import OrderedDict

from .digits import is_integer

__all__ = ['alias_priority', 'bib_priority', 'repo_priority', 'sortOD']


def alias_priority(name, attr):
    if name == attr:
        return 0
    return 1


def bib_priority(attr):
    if attr.get('secondary', False):
        if 'bibcode' in attr:
            if is_integer(attr['bibcode'][:4]):
                return (3000 - int(attr['bibcode'][:4]), '')
        if 'name' in attr:
            return (0, attr['name'])
        return (0, '')
    if 'bibcode' in attr:
        if is_integer(attr['bibcode'][:4]):
            return (-int(attr['bibcode'][:4]), '')
        return (0, '')
    return (0, '')


def repo_priority(attr):
    if is_integer(attr[-4:]):
        return int(attr[-4:])
    return 1000000000


def sortOD(od):
    res = OrderedDict()
    for k, v in sorted(od.items()):
        if isinstance(v, dict):
            res[k] = sortOD(v)
        else:
            res[k] = v
    return res
