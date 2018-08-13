"""
"""

import os
# import sys
# import warnings
import shutil
import glob
# from collections import OrderedDict

# import astrocats
from astrocats import PATHS as AC_PATHS
# from astrocats.structures import struct

import pyastroschema as pas  # noqa


def setup(log):

    # Clear output directory
    pattern = os.path.join(AC_PATHS.SCHEMA_OUTPUT, "", "*.json")
    old_files = sorted(glob.glob(pattern))
    old_files = [os.path.join(AC_PATHS.SCHEMA_OUTPUT, fname) for fname in old_files]
    log.info("Deleting {} old schema files.".format(len(old_files)))
    for fname in old_files:
        os.remove(fname)

    # Copy all `astroschema` schema files to the 'input' directory
    log.info("Copying `astroschema` schema files to input")
    astroschema_fnames = pas.copy_schema_files(AC_PATHS.SCHEMA_INPUT_ASTROSCHEMA, verbose=False)
    for fname in astroschema_fnames:
        base = os.path.basename(fname)
        log.debug("Copied '{}' to '{}'".format(base, fname))

    log.info("Copied {} files".format(len(astroschema_fnames)))

    # Load names of astrocats files (already in 'input')
    log.info("Loading `astrocats` schema file names from input")
    pattern = os.path.join(AC_PATHS.SCHEMA_INPUT_ASTROCATS, "*.json")
    astrocats_fnames = list(sorted(glob.glob(pattern)))
    # Add directory path
    pattern_path = os.path.dirname(pattern)
    astrocats_fnames = [os.path.join(pattern_path, fname) for fname in astrocats_fnames]

    log.info("Loaded {} file names".format(len(astrocats_fnames)))

    def get_output_fname(src, group):
        base = os.path.basename(src)
        if not base.lower().startswith(group.lower()):
            base = group.lower() + "_" + base
        dst = os.path.join(AC_PATHS.SCHEMA_OUTPUT, base)
        return dst

    # Store all source filenames from input and construct the destination names for output
    src_list = []
    dst_list = []
    fnames = [astroschema_fnames, astrocats_fnames]
    groups = ["astroschema", "astrocats"]
    log.info("Copying `astrocats` schema files to input")
    for src_files, grp in zip(fnames, groups):
        for src in src_files:
            src_list.append(src)
            dst = get_output_fname(src, grp)
            # Make sure there are no duplicates
            if dst in dst_list:
                idx = dst_list.index(dst)
                err = "Destination '{}' already exists for '{}', cannot be used for '{}'!".format(
                    dst, src_list[idx], src)
                log.raise_error(err, ValueError)

            dst_list.append(dst)

    # Copy all input schema files to output with modified names
    log.info("Copying {} schema files to output".format(len(src_list)))
    for src, dst in zip(src_list, dst_list):
        shutil.copy(src, dst)
        base = os.path.basename(dst)
        log.debug("Copied '{}' : '{}' ==> '{}'".format(base, src, dst))

    return


'''
def output_schema(path_out, struct_obj, verbose=True):
    # Order of dictionary keys in output file
    sort_order = ['$schema', 'id', 'title', 'description', 'version',
                  'type', 'definitions', 'properties', 'required']

    def sort_item_func(item):
        """Method to generate key for sorting each item of a dictionary
        """
        kk = item[0]

        # Find the index of this key in the `sort_order` list
        idx = None
        for ii, sval in enumerate(sort_order):
            # Ignore text-case
            if kk.lower().startswith(sval.lower()):
                idx = ii
                break

        # Sort things first if found in `sort_order`
        if idx is not None:
            rv = "a_" + str(idx)
        # Then just sort alphabetically by the key
        else:
            rv = "b_" + str(kk)

        return rv

    def sort_dict_func(odict):
        """Method to sort an entire dictionary
        """
        out = OrderedDict(sorted(odict.items(), key=sort_item_func))
        return out

    if hasattr(struct_obj, "_SCHEMA"):
        schema = struct_obj._SCHEMA
    elif isinstance(struct_obj, str):
        if not os.path.exists(struct_obj):
            struct_obj = os.path.join(PATH_SCHEMA_INPUT, struct_obj)
        schema = pas.schema.SchemaDict(struct_obj)
    else:
        err = "Unexptected type '{}' for `struct_obj` '{}'".format(
            type(struct_obj), struct_obj)
        raise RuntimeError(err)

    # schema = struct_obj._SCHEMA
    _title = schema['title'].lower()
    fname = os.path.join(path_out, _title + ".json")
    schema.dump(fname, sort_func=sort_dict_func)
    # if verbose:
    #     print("saved structure '{}' schema to '{}'".format(_title, fname))
    return fname
'''
