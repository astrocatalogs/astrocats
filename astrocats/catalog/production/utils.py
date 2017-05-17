"""Methods relating to production scripts (formerly 'webcat.py').

Copied in from 'webcat.py', 'scripts/events.py'.

"""

import gzip
import os
import hashlib
import json
from collections import OrderedDict


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


def load_md5_file(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def get_event_name_from_filename(fname):
    event_name = os.path.basename(fname)
    # Remove the terminal extension (either 'json' or 'gz')
    event_name = os.path.splitext(event_name)[0]
    # if there is still a json extension, remove that.
    event_name = event_name.replace('.json', '')
    return event_name


def load_event_from_filename(event_fname, log):
    file_text = get_event_text(event_fname)
    event_data = json.loads(file_text, object_pairs_hook=OrderedDict)
    try:
        entry, event_data = list(event_data.items())[0]
    except Exception as err:
        log.error("Error unpacking event from '{}'".format(event_fname))
        log.error("Error: {}".format(str(err)))
        log.error("Contents: {}".format(event_data))
        raise

    return entry, event_data
