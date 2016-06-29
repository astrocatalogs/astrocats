'''Key sorting functions
'''

__all__ = ['alias_priority']


def alias_priority(name, attr):
    if name == attr:
        return 0
    return 1
