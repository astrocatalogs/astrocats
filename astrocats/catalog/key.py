"""
"""
from astrocats.catalog.utils import is_number

from past.builtins import basestring


class KeyCollection(object):
    """General container class with methods to list attribute names and values.

    Used mostly by different `CatDict` subclasses to contain the 'keys' to
    their internal dictionaries.
    """

    _keys = []
    _vals = []
    _compare_vals = []

    @classmethod
    def keys(cls):
        """Return this class's attribute names (those not stating with '_').

        Also retrieves the attributes from base classes, e.g.
        For: ``ENTRY(KeyCollection)``, ``ENTRY.keys()`` gives just the
             attributes of `ENTRY` (`KeyCollection` has no keys).
        For: ``SUPERNOVA(ENTRY)``, ``SUPERNOVA.keys()`` gives both the
             attributes of `SUPERNOVAE` itself, and of `ENTRY`.

        Returns
        -------
        _keys : list of str
            List of names of internal attributes.  Order is effectiely random.
        """
        if cls._keys:
            return cls._keys

        # If `_keys` is not yet defined, create it
        # ----------------------------------------
        _keys = []
        # get the keys from all base-classes aswell (when this is subclasses)
        for mro in cls.__bases__:
            # base classes below `KeyCollection` (e.g. `object`) wont work
            if issubclass(mro, KeyCollection):
                _keys.extend(mro.keys())

        # Get the keys from this particular subclass
        # Only non-hidden (no '_') and variables (non-callable)
        _keys.extend([
            kk for kk in vars(cls).keys()
            if not kk.startswith('_') and not callable(getattr(cls, kk))
        ])
        # Store for future retrieval
        cls._keys = _keys
        return cls._keys

    @classmethod
    def vals(cls):
        """Return this class's attribute values (those not stating with '_').

        Returns
        -------
        _vals : list of objects
            List of values of internal attributes.  Order is effectiely random.
        """
        if cls._vals:
            return cls._vals

        # If `_vals` is not yet defined, create it
        # ----------------------------------------
        _vals = []
        # get the keys from all base-classes aswell (when this is subclasses)
        for mro in cls.__bases__:
            # base classes below `KeyCollection` (e.g. `object`) wont work
            if issubclass(mro, KeyCollection):
                _vals.extend(mro.vals())

        # Get the keys from this particular subclass
        # Only non-hidden (no '_') and variables (non-callable)
        _vals.extend([
            vv for kk, vv in vars(cls).items()
            if not kk.startswith('_') and not callable(getattr(cls, kk))
        ])
        # Store for future retrieval
        cls._vals = _vals
        return cls._vals

    @classmethod
    def compare_vals(cls, sort=True):
        """Return this class's attribute values (those not stating with '_'),
        but only for attributes with `compare` set to `True`.

        Returns
        -------
        _compare_vals : list of objects
            List of values of internal attributes to use when comparing
            `CatDict` objects. Order sorted by `Key` priority, followed by
            alphabetical.
        """
        if cls._compare_vals:
            return cls._compare_vals

        # If `_compare_vals` is not yet defined, create it
        # ----------------------------------------
        _compare_vals = []
        # get the keys from all base-classes aswell (when this is subclasses)
        for mro in cls.__bases__:
            # base classes below `KeyCollection` (e.g. `object`) wont work
            if issubclass(mro, KeyCollection):
                _compare_vals.extend(mro.compare_vals(sort=False))

        # Get the keys from this particular subclass
        # Only non-hidden (no '_') and variables (non-callable)
        _compare_vals.extend([
            vv for kk, vv in vars(cls).items()
            if (not kk.startswith('_') and not callable(getattr(cls, kk)) and
                vv.compare)
        ])

        # Sort keys based on priority, high priority values first
        if sort:
            _compare_vals = sorted(
                _compare_vals,
                reverse=True,
                key=lambda key: (key.priority, key.name))

        # Store for future retrieval
        cls._compare_vals = _compare_vals
        return cls._compare_vals

    @classmethod
    def get_key_by_name(cls, name):
        for val in cls.vals():
            if name == val.name:
                return val
        return Key(name)


class KEY_TYPES(KeyCollection):
    NUMERIC = 'numeric'
    TIME = 'time'
    STRING = 'string'
    BOOL = 'bool'
    DICT = 'dict'
    LIST = 'list'
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

    def __new__(cls,
                name,
                type=None,
                no_source=False,
                listable=False,
                compare=True,
                priority=0,
                kind_preference=[],
                replace_better=False,
                **kwargs):
        """Construct the underlying str object with `name`.
        """
        return str.__new__(cls, name)

    def __init__(self,
                 name,
                 type=None,
                 no_source=False,
                 listable=False,
                 compare=True,
                 priority=0,
                 kind_preference=[],
                 replace_better=False,
                 **kwargs):
        super(Key, self).__init__()
        # Make sure type is allowed
        if type is not None and type not in KEY_TYPES.vals():
            raise ValueError("Key `type` ('{}') must be 'None' or one of '{}'".
                             format(type, KEY_TYPES.keys()))
        self.name = str(name)
        self.type = type
        self.listable = listable
        self.compare = compare
        self.no_source = no_source
        self.priority = priority
        self.kind_preference = kind_preference
        self.replace_better = replace_better
        for key, val in kwargs.items():
            setattr(self, key, val)

    def pretty(self):
        """Return a 'pretty' string representation of this `Key`.

        note: do not override the builtin `__str__` or `__repr__` methods!
        """
        retval = ("Key(name={}, type={}, listable={}, compare={}, "
                  "priority={}, kind_preference={}, "
                  "replace_better={})").format(
                      self.name, self.type, self.listable, self.compare,
                      self.priority, self.kind_preference, self.replace_better)
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
        elif (self.type == KEY_TYPES.TIME and
              not is_number(val) and '-' not in val and '/' not in val):
            return False
        elif self.type == KEY_TYPES.STRING:
            # If its a list, check first element
            if is_list:
                if not isinstance(val[0], basestring):
                    return False
            # Otherwise, check it
            elif not isinstance(val, basestring):
                return False
        elif self.type == KEY_TYPES.BOOL:
            if is_list and not isinstance(val[0], bool):
                return False
            elif not isinstance(val, bool):
                return False

        return True
