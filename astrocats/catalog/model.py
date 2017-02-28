"""Class for representing models.
"""

from astrocats.catalog.catdict import CatDict, CatDictError
from astrocats.catalog.error import Error
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection
from astrocats.catalog.realization import Realization


class MODEL(KeyCollection):
    # Strings
    SOURCE = Key('source', KEY_TYPES.STRING, compare=False)
    CODE = Key('code', KEY_TYPES.STRING)
    NAME = Key('name', KEY_TYPES.STRING)
    VERSION = Key('version', KEY_TYPES.STRING)
    DATE = Key('date', KEY_TYPES.STRING)
    DESCRIPTION = Key('description', KEY_TYPES.STRING, compare=False)
    # Numbers
    ALIAS = Key('alias', KEY_TYPES.NUMERIC, compare=False)
    STEPS = Key('steps', KEY_TYPES.NUMERIC, compare=False)
    # Arrays/dictionaries
    CONVERGENCE = Key('convergence', compare=False)
    REALIZATIONS = Key('realizations', compare=False)
    SCORE = Key('score', compare=False)
    SETUP = Key('setup', compare=False)


class Model(CatDict):
    """Container for a model with associated metadata.

    `Source` citation required.
    """

    _ALLOW_UNKNOWN_KEYS = True
    _KEYS = MODEL

    def __init__(self, parent, **kwargs):
        self._REQ_KEY_SETS = [[MODEL.SOURCE], [MODEL.ALIAS],
                              [MODEL.NAME, MODEL.CODE]]
        # Note: `_check()` is called at end of `super().__init__`
        super(Model, self).__init__(parent, **kwargs)
        self.catalog = parent.catalog

        return

    def _check(self):
        """

        """
        # Run the super method
        super(Model, self)._check()

        err_str = None

        # Future checks go here

        if err_str is not None:
            raise ValueError(err_str)

        return

    def _clean_value_for_key(self, key, value):
        value = super(Model, self)._clean_value_for_key(key, value)

        return value

    def _init_cat_dict(self, cat_dict_class, key_in_self, **kwargs):
        """Initialize a CatDict object, checking for errors.
        """
        # Catch errors associated with crappy, but not unexpected data
        try:
            new_entry = cat_dict_class(self, key=key_in_self, **kwargs)
        except CatDictError as err:
            if err.warn:
                self._log.info("'{}' Not adding '{}': '{}'".format(self[
                    self._KEYS.NAME], key_in_self, str(err)))
            return None
        return new_entry

    def _add_cat_dict(self,
                      cat_dict_class,
                      key_in_self,
                      check_for_dupes=True,
                      **kwargs):
        """Add a CatDict to this Entry if initialization succeeds and it
        doesn't already exist within the Entry.
        """
        # Try to create a new instance of this subclass of `CatDict`
        new_entry = self._init_cat_dict(cat_dict_class, key_in_self, **kwargs)
        if new_entry is None:
            return False

        # Compare this new entry with all previous entries to make sure is new
        if cat_dict_class != Error:
            for item in self.get(key_in_self, []):
                if new_entry.is_duplicate_of(item):
                    item.append_sources_from(new_entry)
                    # Return the entry in case we want to use any additional
                    # tags to augment the old entry
                    return new_entry

        self.setdefault(key_in_self, []).append(new_entry)

        return True

    def add_realization(self, **kwargs):
        self._add_cat_dict(Realization, self._KEYS.REALIZATIONS, **kwargs)

    def sort_func(self, key):
        if key == self._KEYS.SOURCE:
            return 'zzy'
        if key == self._KEYS.ALIAS:
            return 'zzz'
        return key
