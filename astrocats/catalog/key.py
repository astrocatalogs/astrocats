"""
"""
from astrocats.catalog.utils import is_number


class KeyCollection:
    """General container class with methods to list attribute names and values.

    Used mostly by different `CatDict` subclasses to contain the 'keys' to
    their internal dictionaries.
    """
    _keys = []
    _vals = []

    @classmethod
    def keys(cls):
        """Return this class's attribute names (those not stating with '_').

        Returns
        -------
        _keys : list of str
            List of names of internal attributes.  Order is effectiely random.
        """
        if cls._keys:
            return cls._keys

        cls._keys = [kk for kk in vars(cls).keys() if not kk.startswith('_')]
        return cls._keys

    @classmethod
    def vals(cls):
        """Return this class's attribute values (those not stating with '_').

        Returns
        -------
        _vals : list of objects
            List of values of internal attributes.  Order is effectiely random,
            but should match that returned by `keys()`.
        """
        if cls._vals:
            return cls._vals

        cls._vals = [vv for kk, vv in
                     vars(cls).items() if not kk.startswith('_')]
        return cls._vals


class KEY_TYPES(KeyCollection):
    NUMERIC = 'numeric'
    STRING = 'string'
    BOOL = 'bool'
    ANY = None


class Key(str):
    """Class to act as a 'key' (with metadata) to a `CatDict` dictionary.

    Used in `KeyCollection` subclasses, for example `PHOTOMETRY`. The most
    basic function is to contain the str 'key' for a `CatDict` (or other)
    dictionary, but can also contain (and manage) additional specifictions and
    metadata. For example, if a `type` parameter is specified (from
    `KEY_TYPES`), then the builtin `check()` method will ensure that a value is
    consistent with that type.  Arbitrary additional attributes can also be
    stored.

    Attributes
    ----------
    name : str
        The str used as a key for the corresponding dictionary.
    type : str or None (one of `KEY_TYPES`)
        Specification of variable type for the value corresponding to this key.
        'None' means any type is valid.
        NOTE: if ``type == None``, then the `listable` check is disabled.
        note: if `type` is not one of `KEY_TYPES`, a `ValueError` is raised.
    listable : bool
        Whether the value corresponding to this key is allowed to be a list.
        NOTE: this is only enforced if `type` is not 'None'.
    compare : bool
        Whether the value corresponding to this key should be compared to that
        of other dictionaries when looking for duplicate entries.
        i.e. is equality of the value corresponding to this key requisite for
        a different dictionary to be considered a duplicate.

    """
    def __new__(cls, name, type=None, listable=False, compare=True, **kwargs):
        """Construct the underlying str object with `name`.
        """
        return str.__new__(cls, name)

    def __init__(self, name, type=None,
                 listable=False, compare=True, **kwargs):
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

    def pretty(self):
        """Return a 'pretty' string representation of this `Key`.

        note: do not override the builtin `__str__` or `__repr__` methods!
        """
        retval = "Key(name={}, type={}, listable={}, compare={})".format(
            self.name, self.type, self.listable, self.compare)
        return retval

    def check(self, val):
        """Make sure given value is consistent with this `Key` specification.

        NOTE: if `type` is 'None', then `listable` also is *not* checked.
        """
        # If there is no `type` requirement, everything is allowed
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
