"""
"""

import os
import sys
import warnings
import shutil
from collections import OrderedDict

import astrocats
from astrocats.catalog import struct

_PAS_PATH = "/Users/lzkelley/Research/catalogs/astroschema"
if _PAS_PATH not in sys.path:
    sys.path.append(_PAS_PATH)

import pyastroschema as pas  # noqa

PATH_SCHEMA_INPUT = os.path.join(astrocats._PATH_SCHEMA, "input", "")
PATH_SCHEMA_OUTPUT = os.path.join(astrocats._PATH_SCHEMA, "output", "")

# This describes replacements/renames that should be applied to the file corresponding to the
# given schema name
_REPLACE = {
    'quantity': [
        ["error_value", "e_value"],
        ["error_lower", "e_lower_value"],
        ["error_upper", "e_upper_value"],
        ["units_value", "u_value"],
        ["units_error", "u_e_value"],
    ]
}


def main(add_files=[], add_structures=[]):

    # Copy all `astroschema` schema files to the 'input' directory
    schema_fnames = pas.copy_schema_files(PATH_SCHEMA_INPUT, verbose=True)
    for src in add_files:
        dst = os.path.basename(src)
        dst = os.path.join(PATH_SCHEMA_INPUT, dst)
        shutil.copy(src, dst)
        # print("Copied '{}' ==> '{}'".format(src, dst))

    mod_files = schema_fnames + add_files

    # Make any schema-file modifications
    warnings.warn("WARNING: `modify_schema_file` might not take effect until a second run...")
    for sch_file in mod_files:
        modify_schema_file(sch_file)

    # Save all schema to output directory
    output_structures = struct.STRUCTURES + add_structures
    for obj in output_structures:
        output_schema(PATH_SCHEMA_OUTPUT, obj)

    return


def modify_schema_file(fname):

    if not os.path.exists(fname):
        raise RuntimeError("Schema file '{}' does not exist!".format(fname))

    def do_replace(fil, mapping):
        with open(fil, 'r') as inn:
            schema = inn.read()

        for old, new in mapping:
            schema = schema.replace(old, new)

        with open(fil, 'w') as out:
            out.write(schema)

        return

    # Get the name of the schema itself
    schema_name = os.path.basename(fname).lower().split('.json')[0]
    # If there is any entry for this schema, do the replacement
    if schema_name in _REPLACE:
        # print("Running replace in '{}'".format(schema_name))
        # Load the mapping (src-->dst) for replacements
        mapping = _REPLACE[schema_name]
        # Do the replacing
        do_replace(fname, mapping)

    return


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
