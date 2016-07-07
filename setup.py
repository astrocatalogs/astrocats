import os

from setuptools import setup


# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="astrocats",
    version="0.1.0",
    author="James Guillochon",
    author_email="guillochon@gmail.com",
    description=("Package for downloading, analyzing, and constructing open "
                 "astronomy catalogs."),
    license="MIT",
    keywords="astronomy",
    url="http://packages.python.org/astrocats",
    packages=['astrocats'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha"
    ],
)
