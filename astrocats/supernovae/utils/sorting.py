"""
"""
from ..constants import PREF_KINDS

__all__ = ['event_attr_priority', 'frame_priority']


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
