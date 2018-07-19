"""Class for representing quantities."""

from . import struct

Quantity = struct._Quantity
QUANTITY = Quantity.get_keychain()
Quantity._KEYS = QUANTITY
