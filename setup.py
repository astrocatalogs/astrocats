"""WARNING : THIS SCRIPT IS NOT CURRENTLY OPERATIONAL.
"""

import logging
import os

from setuptools import find_packages, setup
from setuptools.command.develop import develop
from setuptools.command.install import install

with open('requirements.txt') as f:
    required = f.read().splitlines()

dir_path = os.path.dirname(os.path.realpath(__file__))
exec(open(os.path.join(dir_path, 'astrocats', '__init__.py')).read())


def read(fname):
    return open(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), fname)).read()


def setup_uc():
    from astrocats.main import setup_user_config

    setup_user_config(logging.getLogger())


class PostDevelopCommand(develop):
    """Post-develop command."""

    def run(self):
        setup_uc()
        develop.run(self)


class PostInstallCommand(install):
    """Post-installation command."""

    def run(self):
        setup_uc()
        install.run(self)


setup(
    name="astrocats",
    packages=find_packages(exclude=('*supernovae*', '*tidaldisruptions*',
                                    '*novae*')),
    include_package_data=True,
    version=__version__,  # noqa
    description=("Package for downloading, analyzing, and constructing open "
                 "astronomy catalogs."),
    license=__license__,  # noqa
    author=__author__,  # noqa
    author_email="guillochon@gmail.com",
    install_requires=required,
    setup_requires=required,
    url="https://github.com/astrocatalogs/astrocats",
    download_url=(
        'https://github.com/astrocatalogs/astrocats/tarball/' +
        __version__  # noqa
    ),
    cmdclass={
        'develop': PostDevelopCommand,
        'install': PostInstallCommand,
    },
    keywords="astronomy",
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.5"
    ])
