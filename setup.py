import os
import uuid

from setuptools import setup, find_packages
from pip.req import parse_requirements


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
    version="0.1.6",
    author="James Guillochon",
    author_email="guillochon@gmail.com",
    description=("Package for downloading, analyzing, and constructing open "
                 "astronomy catalogs."),
    license="MIT",
    keywords="astronomy",
    url="https://github.com/astrocatalogs/astrocats",
    packages=find_packages(),
    long_description=read('README.md'),
    install_requires=reqs,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.5"
    ],
    zip_safe=True
)
