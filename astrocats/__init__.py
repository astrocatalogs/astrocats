"""Astrocats: Scripts for creating and analyzing catalogs of astronomical data.
"""

import os
import sys

__version__ = '0.2.6'
__author__ = 'James Guillochon'
__license__ = 'MIT'

# Set the path for the user's configuration file
_CONFIG_PATH = os.path.join(os.path.expanduser('~'),
                            '.config', 'astrocats', 'astrocatsrc')

# Make sure configuration file exists, or that's what we're doing
# (with the 'setup' subcommand)
if not os.path.isfile(_CONFIG_PATH) and 'setup' not in sys.argv:
    raise RuntimeError("'{}' does not exist.  "
                       "Run `astrocats setup` to configure."
                       "".format(_CONFIG_PATH))
