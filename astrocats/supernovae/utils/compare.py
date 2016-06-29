'''Utility functions related to comparing quanta to one another to determine if
they are unique.
'''
from cdecimal import Decimal


def same_tag_num(photo, val, tag, canbelist=False):
    issame = (
        (tag not in photo and not val) or
        (tag in photo and not val) or
        (tag in photo and
         ((not canbelist and Decimal(photo[tag]) == Decimal(val)) or
          (canbelist and
           ((isinstance(photo[tag], str) and isinstance(val, str) and
             Decimal(photo[tag]) == Decimal(val)) or
            (isinstance(photo[tag], list) and isinstance(val, list) and
             photo[tag] == val))))))
    return issame


def same_tag_str(photo, val, tag):
    issame = ((tag not in photo and not val) or (
        tag in photo and not val) or (tag in photo and photo[tag] == val))
    return issame
