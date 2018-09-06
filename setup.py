"""WARNING : THIS SCRIPT IS NOT CURRENTLY OPERATIONAL.
"""

# import logging
# import os

from setuptools import find_packages, setup
# from setuptools.command.develop import develop
# from setuptools.command.install import install

with open('requirements.txt') as f:
    required = f.read().splitlines()

with open("README.md", "r") as inn:
    long_description = inn.read().strip()

with open("VERSION", "r") as inn:
    version = inn.read().strip()

'''
dir_path = os.path.dirname(os.path.realpath(__file__))
exec(open(os.path.join(dir_path, 'astrocats', '__init__.py')).read())

def read(fname):
    return open(os.path.join(os.path.dirname(os.path.abspath(__file__)), fname)).read()

def setup_uc():
    from astrocats.main import setup_user_config
    setup_user_config(logging.getLogger())


class PostDevelopCommand(develop):
    """Post-develop command."""

    def run(self):
        # setup_uc()
        develop.run(self)


class PostInstallCommand(install):
    """Post-installation command."""

    def run(self):
        # setup_uc()
        install.run(self)
'''

setup(
    name="astrocats",
    packages=find_packages(exclude=('*tidaldisruptions*', '*novae*', '*faststars*')),
    include_package_data=True,
    version=version,
    description="Package for downloading, analyzing, and constructing open astronomy catalogs.",
    license="MIT",
    author="James Guillochon & Luke Zoltan Kelley",
    author_email="guillochon@gmail.com",
    install_requires=required,
    setup_requires=required,
    url="https://github.com/astrocatalogs/astrocats",
    download_url=('https://github.com/astrocatalogs/astrocats/tarball/' + version),
    entry_points={
        'console_scripts': [
            'astrocats = astrocats.__main__:main'
        ]
    },
    # cmdclass={
    #     'develop': PostDevelopCommand,
    #     'install': PostInstallCommand,
    # },
    keywords="astronomy",
    long_description=long_description,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.5"
    ]
)
