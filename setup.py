"""WARNING : THIS SCRIPT IS NOT CURRENTLY OPERATIONAL.
"""

import os
import re
import uuid

from pip.req import parse_requirements
from setuptools import find_packages, setup

raise(SystemExit('Setup not setup yet.'))

VERSIONFILE = "astrocats/__init__.py"
ver_file = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, ver_file, re.M)

if mo:
    version = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE))

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

install_reqs = parse_requirements('requirements.txt', session=uuid.uuid1())
reqs = [str(req.req) for req in install_reqs]

setup(
    name="astrocats",
    version=version,
    author="James Guillochon",
    author_email="guillochon@gmail.com",
    description=("Package for downloading, analyzing, and constructing open "
                 "astronomy catalogs."),
    license="MIT",
    keywords="astronomy",
    url="https://github.com/astrocatalogs/astrocats",
    packages=find_packages(),
    # scripts=['astrocats/main.py'],
    entry_points={'console_scripts': 'astrocats = astrocats.main:main'},
    long_description=read('README.md'),
    install_requires=reqs,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.5"
    ],
    zip_safe=True
)
