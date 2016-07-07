"""Astrocats: Scripts for creating and analyzing catalogs of astronomical data.
"""

import os

__version__ = '0.1.7'
__author__ = 'James Guillochon'
__license__ = 'MIT'

_CONFIG_PATH = os.path.join(os.path.expanduser('~'),
                            '.config', 'astrocats', 'astrocatsrc')

if not os.path.isfile(_CONFIG_PATH):
    raise RuntimeError("'{}' does not exist.  "
                       "Run `astrocats setup` to configure."
                       "".format(_CONFIG_PATH))
