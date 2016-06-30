"""
"""
from astrocats.catalog.utils import is_number


class KeyCollection:
    """
    """
    @classmethod
    def keys(cls):
        _keys = [kk for kk in vars(cls).keys() if not kk.startswith('_')]
        return _keys

    @classmethod
    def vals(cls):
        _vals = [vv for kk, vv in vars(cls).items() if not kk.startswith('_')]
        return _vals


class KEY_TYPES(KeyCollection):
    NUMERIC = 'numeric'
    STRING = 'string'
    BOOL = 'bool'
    ANY = None


class Key(str):
    def __new__(cls, name, type=None, listable=False, **kwargs):
        return str.__new__(cls, name)

    def __init__(self, name, type=None, listable=False, **kwargs):
        # Make sure type is allowed
        if type is not None and type not in KEY_TYPES.vals():
            raise ValueError(
                "Key `type` ('{}') must be 'None' or one of '{}'".format(
                    type, KEY_TYPES.keys()))
        self.name = str(name)
        self.type = type
        self.listable = listable
        for key, val in kwargs.items():
            setattr(self, key, val)

    def __repr__(self):
        retval = "Key(name={}, type={}, cambelist={})".format(
            self.name, self.key, self.listable)
        return retval

    def check(self, val):
        """Make sure given value is consistent with this `Key` specification.
        """
        is_list = isinstance(val, list)
        # If lists are not allowed, and this is a list --> false
        if not self.listable and is_list:
            return False

        # If there is a type requirement, check that
        if self.type is not None:
            # `is_number` already checks for either list or single value
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
                if is_list and not isinstance(val[0], bool):
                    return False
                elif not isinstance(val, bool):
                    return False

        return True
