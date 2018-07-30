"""Class for representing model realizations."""

from . import struct

Realization = struct._Realization
REALIZATION = Realization.get_keychain()
Realization._KEYS = REALIZATION
