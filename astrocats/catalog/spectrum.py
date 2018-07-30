"""Class for representing spectra."""

from . import struct

Spectrum = struct._Spectrum
SPECTRUM = Spectrum.get_keychain(extendable=True)
Spectrum._KEYS = SPECTRUM
