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
    """

    NOTE: if `type` is 'None', then `listable` also is *not* checked.
    """
    def __new__(cls, name, type=None, listable=False, compare=True, **kwargs):
        return str.__new__(cls, name)

    def __init__(self, name, type=None, listable=False, compare=True, **kwargs):
        # Make sure type is allowed
        if type is not None and type not in KEY_TYPES.vals():
            raise ValueError(
                "Key `type` ('{}') must be 'None' or one of '{}'".format(
                    type, KEY_TYPES.keys()))
        self.name = str(name)
        self.type = type
        self.listable = listable
        self.compare = compare
        for key, val in kwargs.items():
            setattr(self, key, val)

    def __repr__(self):
        retval = "Key(name={}, type={}, listable={}, compare={})".format(
            self.name, self.key, self.listable, self.compare)
        return retval

    def check(self, val):
        """Make sure given value is consistent with this `Key` specification.

        NOTE: if `type` is 'None', then `listable` also is *not* checked.
        """
        # If there is no `type` requirement, everything is okay
        if self.type is None:
            return True

        is_list = isinstance(val, list)
        # If lists are not allowed, and this is a list --> false
        if not self.listable and is_list:
            return False

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
