'''Key sorting functions
'''

from .digits import is_integer

__all__ = ['alias_priority', 'bib_priority', 'repo_priority']


def alias_priority(name, attr):
    if name == attr:
        return 0
    return 1


def bib_priority(attr):
    if attr.get('secondary', False):
        if 'bibcode' in attr:
            if is_integer(attr['bibcode'][:4]):
                return '0%d' % (-3000+int(attr['bibcode'][:4]))
        if 'name' in attr:
            return attr['name']
        return ''
    if 'bibcode' in attr:
        if is_integer(attr['bibcode'][:4]):
            return '0%d' % -int(attr['bibcode'][:4])
        return '0%d' % 0
    return '0%d' % 0


def repo_priority(attr):
    if is_integer(attr[-4:]):
        return int(attr[-4:])
    return 1000000000
