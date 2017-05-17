# -*- coding: utf-8 -*-
"""Class defitions for `CatDict` and `CatDictError`, data storage classes."""
from collections import OrderedDict
from copy import deepcopy

from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection
from astrocats.catalog.utils import uniq_cdl

try:
    basestring
except NameError:
    basestring = str


class CatDictError(Exception):
    """Special Error class for non-fatal errors raised in CatDict."""

    def __init__(self, *args, **kwargs):
        """Initialize `CatDictError`."""
        # If `warn` is True, then a warning should be issues.  Otherwise ignore
        # completely
        self.warn = True
        if 'warn' in kwargs:
            self.warn = kwargs.pop('warn')
        Exception.__init__(self, *args, **kwargs)


class CatDict(OrderedDict):
    """General data storage super-class used throughout catalogs.

    In general, `CatDict` subclasses are used to manage/store data which ends
    up in the lowest-level dictionaries saved to output json files.  For
    example, and individual `Source` or `Photometry` entry.

    Attributes
    ----------
    These attributes should be overridden in subclasses as needed:
    _KEYS : `KeyCollection` object
        Contains the `Key` strings and associated specifications which makeup
        the keys for data in this object.  For example: `PHOTOMETRY` for the
        `Photometry` `CatDict` subclass.
    _ALLOW_UNKNOWN_KEYS : bool
        Whether data with keys *not* in the `_KEYS` object are allowed.
        'True': additional keyword arguments passed to the constructor are
            add without validation or sanitization.
        'False': keyword arguments passed to the constructor which are not in
            `ENTRY` will produce errors.
    _REQ_KEY_SETS : list of lists
        Which elements of the associated `_KEYS` are required for each instance
        of this class.
        The structure of this variable is a list of lists, where each inner
        list contains a set of 'Key's, *at least one* of which are required for
        validity.  For example, if
        ``_REQ_KEY_SETS = [[_KEYS.ONE, _KEYS.TWO], [_KEYS.THREE]]``
        then either `_KEYS.ONE` *or* `_KEYS.TWO` is required, and `_KEYS.THREE`
        if also required.

    Notes
    -----
    -   Invalid data and Errors
        If, for any reason, the `CatDict` being constructed looks invalid in a
        not unexpected way (e.g. required arguments are missing), then a
        `CatDictError` exception should be raised.  The method in the parent
        class (e.g. `Supernova`) which tries to construct the `CatDict` (e.g.
        `add_photometry`) should catch those errors specifically and deal
        with them appropriately.

    """

    # The `KeyCollection` object associated with this dictionary For example,
    # `PHOTOMETRY` for the `Photometry(CatDict)` subclass
    _KEYS = KeyCollection

    # If `_ALLOW_UNKNOWN_KEYS` is 'False', then only parameters with names
    # included in `_KEYS` are allowed.  Others will raise an error. If this
    # parameter is 'True', then parameters corresponding to those in `_KEYS`
    # are still checked (for type etc), but additional parameters are just
    # tacked onto the `CatDict` object without any checks or errors.
    _ALLOW_UNKNOWN_KEYS = True

    _REQ_KEY_SETS = []

    def __init__(self, parent, key=None, **kwargs):
        """Initialize `CatDict`."""
        super(CatDict, self).__init__()
        # Store the parent object (an `Entry` subclass) to which this instance
        # will belong.  e.g. a `Supernova` entry.
        self._parent = parent
        self._key = key
        self._log = parent.catalog.log

        # Store any individual keys which are required
        self._req_keys = []
        for rks in self._REQ_KEY_SETS:
            # If this set is only 1 long, that key is individually required
            if len(rks) == 1:
                self._req_keys.append(rks[0])

        # Iterate over all `_KEYS` parameters, load each if given note that the
        # stored 'values' are the `Key` objects, referred to here with the name
        # 'key'.
        vals = self._KEYS.vals()
        for key in kwargs.copy():
            # If we allow unknown keys, or key is in list of known keys,
            # process and store it.
            kiv = key in vals
            if self._ALLOW_UNKNOWN_KEYS or kiv:
                # Load associated Key object if it exists, otherwise construct
                # a default Key object.
                if kiv:
                    key_obj = vals[vals.index(key)]
                else:
                    self._log.info('[{}] `{}` not in list of keys for `{}`, '
                                   'adding anyway as allow unknown keys is '
                                   '`{}`.'.format(parent[
                                       parent._KEYS.NAME], key,
                                       type(self).__name__,
                                       self._ALLOW_UNKNOWN_KEYS))
                    key_obj = Key(key)

                # Handle Special Cases
                # --------------------
                # Only keep booleans and strings if they evaluate true.
                if ((key_obj.type == KEY_TYPES.BOOL or
                     key_obj.type == KEY_TYPES.STRING) and not kwargs[key]):
                    del kwargs[key]
                    continue

                # Make sure value is compatible with the 'Key' specification.
                check_fail = False
                if not key_obj.check(kwargs[key]):
                    check_fail = True
                    self._log.info("Value for '{}' is invalid "
                                   "'{}':'{}'".format(key_obj.pretty(), key,
                                                      kwargs[key]))
                    # Have the parent log a warning if this is a required key
                    if key in self._req_keys:
                        raise CatDictError(
                            "Value for required key '{}' is invalid "
                            "'{}:{}'".format(key_obj.pretty(), key,
                                             kwargs[key]),
                            warn=True)

                # Check and store values
                # ----------------------
                # Remove key-value pair from `kwargs` dictionary.
                value = kwargs.pop(key)
                value = self._clean_value_for_key(key_obj, value)
                # only store values that are not empty
                if value and not check_fail:
                    self[key] = value

        # If we require all parameters to be a key in `PHOTOMETRY`, then all
        # elements should have been removed from `kwargs`.
        if not self._ALLOW_UNKNOWN_KEYS and len(kwargs):
            raise CatDictError(
                "All permitted keys stored, remaining: '{}'".format(kwargs))

        # Make sure that currently stored values are valid
        self._check()

        return

    def __deepcopy__(self, memo):
        dict_copy = OrderedDict()
        for key in self:
            if not key.startswith('__'):
                dict_copy[key] = deepcopy(self[key])
        return self.__class__(self._parent, key=self._key,
                              **dict_copy)

    def sort_func(self, key):
        return key

    def pretty(self):
        retval = "{}({}, Parent:{})".format(
            type(self), self._key, self._parent)
        return retval

    def is_duplicate_of(self, other):
        # If these are not the same type, return False
        if type(other) is not type(self):
            return False

        # Go over all expected parameters and check equality of each
        for key in self._KEYS.compare_vals():
            kis = key in self
            kio = key in other
            # If only one object has this parameter, not the same
            if kis != kio:
                return False
            # If self doesnt have this parameter (and thus neither does), skip
            if not kis:
                continue

            # Now, both objects have the same parameter, compare them
            if self[key] != other[key]:
                return False

        return True

    def append_sources_from(self, other):
        """Merge the source alias lists of two CatDicts."""
        # Get aliases lists from this `CatDict` and other
        self_aliases = self[self._KEYS.SOURCE].split(',')
        other_aliases = other[self._KEYS.SOURCE].split(',')

        # Store alias to `self`
        self[self._KEYS.SOURCE] = uniq_cdl(self_aliases + other_aliases)

        return

    def _check(self):
        """
        """
        for req_any in self._REQ_KEY_SETS:
            if not any([req_key in self for req_key in req_any]):
                err_str = ("'{}' Requires one or more of: ".format(self._key) +
                           ",".join("'{}'".format(rk) for rk in req_any))
                self._log.info(err_str)
                raise CatDictError(err_str)

        return

    def _clean_value_for_key(self, key, value):
        """

        FIX: should this be in `key`??  like 'check()'??
        """
        if key.type is None:
            return value

        # Store whether given value started as a list or not
        single = True
        # Make everything a list for conversions below
        if isinstance(value, list):
            # But if lists arent allowed, and this is, raise error
            if not key.listable:
                raise CatDictError("`value` '{}' for '{}' shouldnt be a list.".
                                   format(value, key.pretty()))

            single = False
        else:
            single = True
            value = [value]

        # Store booleans as booleans, make sure each element of list is bool
        if key.type == KEY_TYPES.BOOL:
            if not all(isinstance(val, bool) for val in value):
                raise CatDictError("`value` '{}' for '{}' should be boolean".
                                   format(value, key.pretty()))
        # Strings and numeric types should be stored as strings
        elif key.type in [KEY_TYPES.STRING, KEY_TYPES.NUMERIC, KEY_TYPES.TIME]:
            # Clean leading/trailing whitespace
            value = [
                val.strip() if isinstance(val, (str, basestring)) else str(val)
                for val in value
            ]
            # Only keep values that are not empty
            value = [val for val in value if len(val)]

        # Convert back to single value, if thats how it started
        if single and len(value):
            value = value[0]

        return value
