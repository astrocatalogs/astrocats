"""
"""
import os
import shutil
# import sys
import glob

SRC_DIR = "/Users/lzkelley/Research/catalogs/astroschema/schema/"
DST_DIR = "/Users/lzkelley/Research/catalogs/astrocats/astrocats/schema/"

# FILES = ["source.json", "quantity.json", "meta-schema.json"]

# AS_SUBDIR = "astroschema"
AS_SUBDIR = ""

REPLACE = {
    "quantity.json": [
        ["error_value", "e_value"],
        ["error_lower", "e_lower_value"],
        ["error_upper", "e_upper_value"],
        ["units_value", "u_value"],
        ["units_error", "u_e_value"],
    ],
}


def main():
    as_subdir = os.path.join(DST_DIR, AS_SUBDIR, "")
    if not os.path.exists(as_subdir):
        os.makedirs(as_subdir)

    # Get list of source files
    schema_pattern = os.path.join(SRC_DIR, "*.json")
    files = sorted(glob.glob(schema_pattern))
    files = [fs for fs in files if not os.path.split(fs)[-1].startswith('_')]
    # print(files)

    # Copy source files to local subdirectory
    for fs in files:
        base_name = os.path.basename(fs)
        target_name = os.path.join(as_subdir, base_name)
        shutil.copy(fs, target_name)
        print("'{}' ==>\n\t'{}'".format(fs, target_name))

    # Make `REPLACE` changes to local files
    for rep_file, mappings in REPLACE.items():
        target_name = os.path.join(as_subdir, rep_file)
        print(target_name)
        replace_in_file(target_name, mappings)

    return


def replace_in_file(fname, mappings):
    # Read in the file
    filedata = None
    with open(fname, 'r') as fil:
        filedata = fil.read()

    # Replace the target string
    for src, dst in mappings:
        print("\t\t'{}' --> '{}'".format(src, dst))
        filedata = filedata.replace(src, dst)

    # Write the file out again
    with open(fname, 'w') as fil:
        fil.write(filedata)

    return


if __name__ == "__main__":
    main()
