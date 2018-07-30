"""Handle known errors in catalog data."""

from . import struct

Error = struct._Error
ERROR = Error.get_keychain()
Error._KEYS = ERROR
