"""Class for representing models."""

from . import struct

Model = struct._Model
MODEL = Model.get_keychain()
Model._KEYS = MODEL
