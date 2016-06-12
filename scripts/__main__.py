"""Entry point for OCS import scripts.
"""

print("scripts/__main__.py")

from . import importers
from . import utils


if __name__ == "__main__":
    importers.importer.import_main()
