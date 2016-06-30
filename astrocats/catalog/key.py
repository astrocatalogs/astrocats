"""
"""
from astrocats.catalog.utils import is_bool, is_number
import json


class KEY_TYPES:
    NUMERIC = 'numeric'
    STRING = 'string'
    BOOL = 'bool'
    ANY = None

    _keys = sorted([kk for kk in dir() if not kk.startswith('_')])


class Key(str):
    def __new__(cls, string, type=None, canbelist=False):
        return str.__new__(cls, string)

    def __init__(self, string, type=None, canbelist=False):
        # Make sure type is allowed
        if type is not None and type not in KEY_TYPES._keys:
            raise ValueError(
                "Key `type` ('{}') must be 'None' or one of '{}'".format(
                    type, KEY_TYPES._keys))
        self.string = str(string)
        self.type = type
        self.canbelist = canbelist

    def check(self, val):
        """Make sure given value is consistent with this `Key` specification.
        """
        is_list = isinstance(val, list)
        # If lists are not allowed, and this is a list --> false
        if not self.canbelist and is_list:
            return False

        # If there is a type requirement, check that
        if self.type is not None:
            if self.type == KEY_TYPES.NUMERIC and not is_number(val):
                return False
            elif self.type == KEY_TYPES.STRING:
                # If its a list, check first element
                if is_list and not isinstance(val[0], str):
                    return False
                # Otherwise, check it
                elif not isinstance(val, str):
                    return False
            elif self.type == KEY_TYPES.BOOL:
                if is_list and not isinstance(json.loads(val[0]), bool):
                    return False

        return True
