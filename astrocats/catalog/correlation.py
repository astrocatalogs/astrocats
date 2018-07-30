"""Class for representing correlations."""

from . import struct

Correlation = struct._Correlation
CORRELATION = Correlation.get_keychain(extendable=True)
Correlation._KEYS = CORRELATION
