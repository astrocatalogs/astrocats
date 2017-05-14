"""Methods relating to production scripts (formerly 'webcat.py').

Copied in from 'webcat.py', 'scripts/events.py'.

"""

import gzip
import os
import hashlib


def get_event_text(eventfile):
    if eventfile.split('.')[-1] == 'gz':
        with gzip.open(eventfile, 'rt') as f:
            filetext = f.read()
    else:
        with open(eventfile, 'r') as f:
            filetext = f.read()
    return filetext


def get_event_filename(name):
    return name.replace('/', '_')


def touch(fname, times=None):
    with open(fname, 'a'):
        os.utime(fname, times)


def label_format(label):
    newlabel = label.replace('Angstrom', 'Å')
    newlabel = newlabel.replace('^2', '²')
    return newlabel


def get_first_value(catalog, name, field):
    return catalog[name][field][0]['value'] if field in catalog[
        name] and catalog[name][field] else ''


def get_first_kind(catalog, name, field):
    return (catalog[name][field][0]['kind'] if field in catalog[name] and
            catalog[name][field] and 'kind' in catalog[name][field][0] else '')


def md5file(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()
