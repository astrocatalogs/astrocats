"""Class for representing sources of data.
"""
from .key import Key, KEY_TYPES


class SOURCES:
    # Strings
    NAME = Key('name', KEY_TYPES.STRING)
    ALIAS = Key('alias', KEY_TYPES.STRING)
    BIBCODE = Key('bibcode', KEY_TYPES.STRING)
    URL = Key('url', KEY_TYPES.STRING)
    ACKNOWLEDGMENT = Key('acknowledgment', KEY_TYPES.STRING)
    # Booleans
    SECONDARY = Key('secondary', KEY_TYPES.BOOL)
