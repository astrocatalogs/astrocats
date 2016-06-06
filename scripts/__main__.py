"""Entry point for OCS import scripts.
"""

print("scripts/__main__.py")

from . import import_sn
from . import utils


if __name__ == "__main__":
    import_sn.importer.import_main()
