"""Class for representing sources of data."""


from . import struct

Source = struct._Source
SOURCE = Source.get_keychain()
Source._KEYS = SOURCE
