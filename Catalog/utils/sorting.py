'''Key sorting functions
'''

__all__ = ['alias_priority', 'event_attr_priority', 'frame_priority']


PREF_KINDS = ['heliocentric', 'cmb', 'spectroscopic',
              'photometric', 'host', 'cluster', '']


def alias_priority(name, attr):
    if name == attr:
        return 0
    return 1


def event_attr_priority(attr):
    if attr == 'photometry':
        return 'zzy'
    if attr == 'spectra':
        return 'zzz'
    if attr == 'schema':
        return 'aaa'
    if attr == 'name':
        return 'aab'
    if attr == 'sources':
        return 'aac'
    if attr == 'alias':
        return 'aad'
    return attr


def frame_priority(attr):
    if 'kind' in attr:
        if attr['kind'] in PREF_KINDS:
            return PREF_KINDS.index(attr['kind'])
        else:
            return len(PREF_KINDS)
    return len(PREF_KINDS)
