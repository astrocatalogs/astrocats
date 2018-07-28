"""Class for representing photometric data."""

from . import struct

Photometry = struct._Photometry
PHOTOMETRY = Photometry.get_keychain()
Photometry._KEYS = PHOTOMETRY
