"""
"""

import os
# import warnings
import shutil
import re
import glob

from astrocats import PATHS as AC_PATHS

import pyastroschema as pas
import pyastroschema.utils  # noqa

'''
class SchemaCompiler(object):

    # name = None
    # path = None
    # base = None
    # extensions = None
    # updates = None

    def __init__(self, name, base=None, extensions=None, updates=None):

        # if path is None:
        #     path = os.path.realpath(os.curdir)
        # path = os.path.join(path, '')

        if extensions is None:
            extensions = []
        elif not isinstance(extensions, list):
            extensions = [extensions]

        if updates is None:
            updates = []
        elif not isinstance(updates, list):
            updates = [updates]

        # Convert from filenames to dictionaries (as needed)
        #    note: `load_schema_dict` returns (schema, path, title)
        base, _, _ = pas.utils.load_schema_dict(base)[0]
        extensions = [ext if isinstance(ext, dict) else pas.utils.load_schema_dict(ext)[0]
                      for ext in extensions]
        updates = [upd if isinstance(upd, dict) else pas.utils.load_schema_dict(upd)[0]
                   for upd in updates]

        self.name = name
        self.base = base
        # self.path = path

        self.extensions = extensions
        self.updates = updates

        self._schema = None
        self._keychain = None

        return

    def compile(self, check_conflict=True, mutable=False, extendable=True):
        schema_dict = pas.schema.SchemaDict(self.base)
        for ext in self.extensions:
            schema_dict.extend(ext, check_conflict=check_conflict)
        for upd in self.updates:
            schema_dict.update(upd)
        schema = schema_dict
        keychain = pas.keys.Keychain(schema_dict, mutable=mutable, extendable=extendable)
        self._schema = schema
        self._keychain = keychain
        return schema, keychain

    def save(self, path):
        name = self.name
        name = name if name.lower().endswith('.json') else name + '.json'
        if not os.path.exists(path) or not os.path.isdir(path):
            raise ValueError("Path '{}' for '{}' does not exist or is not a directory!".format(
                path, self.name))

        output = os.path.join(path, name)


    @property
    def schema(self):
        if self._schema is None:
            schema, keychain = self.compile()
            return schema

        return self._schema

    @property
    def keychain(self):
        if self._keychain is None:
            schema, keychain = self.compile()
            return keychain

        return self._keychain
'''


def setup(log, catalog_info=None):

    # Path to catalog-specific schema-files
    if catalog_info is not None:
        catalog_name = catalog_info['catalog_name']
        catalog_schema_path = catalog_info.get('schema_path', None)
    else:
        catalog_name = None
        catalog_schema_path = None

    # Clear output directories
    clear_dirs = [AC_PATHS.SCHEMA_OUTPUT]
    for clear in clear_dirs:
        pattern = os.path.join(clear, "", "*.json")
        old_files = sorted(glob.glob(pattern))
        old_files = [os.path.join(clear, fname) for fname in old_files]
        log.info("Deleting {} old schema files in '{}'.".format(len(old_files), clear))
        for fname in old_files:
            os.remove(fname)

    # Copy all `astroschema` schema files to the 'input' directory
    # astroschema_fnames = pas.copy_schema_files(AC_PATHS.SCHEMA_INPUT_ASTROSCHEMA, verbose=False)
    # log.info("Copied {} 'astroschema' files to input".format(len(astroschema_fnames)))
    astroschema_fnames = pas.get_schema_file_names()
    log.info("Loaded {} 'astroschema' file names".format(len(astroschema_fnames)))

    # Load names of astrocats schema files (already in 'input')
    # ------------------------------------------------------------
    pattern = os.path.join(AC_PATHS.SCHEMA_INPUT_ASTROCATS, "*.json")
    astrocats_fnames = list(sorted(glob.glob(pattern)))
    # Add directory component of path to filenames
    pattern_path = os.path.dirname(pattern)
    astrocats_fnames = [os.path.join(pattern_path, fname) for fname in astrocats_fnames]
    log.info("Loaded {} 'astrocats' file names from input".format(len(astrocats_fnames)))

    # Load names of catalog-specific schema files
    # ----------------------------------------------------------
    if catalog_schema_path is None:
        catalog_fnames = []
    else:
        pattern = os.path.join(catalog_schema_path, "*.json")
        catalog_fnames = sorted(glob.glob(pattern))

    log.info("Loaded {} file names from '{}' catalog input".format(
        len(catalog_fnames), catalog_name))

    # Compile Schema
    # ===============================
    schema = []
    match_asc = [False for fn in astrocats_fnames]
    match_cat = [False for fn in catalog_fnames]
    # Go through each `astroschema` schema-files, and create schema
    # --------------------------------------------------------------------------
    for fname_sch in astroschema_fnames:
        base_sch = os.path.basename(fname_sch).lower()
        exts = []
        # Add matching `astrocats` schema-files to extensions list
        for jj, fname_asc in enumerate(astrocats_fnames):
            base_asc = os.path.basename(fname_asc).lower()
            if base_sch == base_asc:
                exts.append(fname_asc)
                match_asc[jj] = True
                break

        # Add matching `catalog` schema-files to extensions list
        for kk, fname_cat in enumerate(catalog_fnames):
            base_cat = os.path.basename(fname_cat).lower()
            if base_sch == base_cat:
                exts.append(fname_cat)
                match_cat[kk] = True
                break

        # Load schema
        sch = pas.schema.SchemaDict(schema=fname_sch)
        log.debug("Created schema for '{}'".format(fname_sch))
        # Apply extensions
        for ex in exts:
            sch.extend(ex)
            log.debug("    Extended with '{}'".format(ex))

        schema.append([sch, base_sch])

    exts = []
    # Go through each `astrocats` schema-file not already matched, and create schema
    # ------------------------------------------------------------------------------------
    for jj, fname_asc in enumerate(astrocats_fnames):
        # Skip if already matched, and used in a schema
        if match_asc[jj]:
            continue

        base_asc = os.path.basename(fname_asc).lower()

        exts = []
        # Add matching `catalog` schema-files to extensions list
        for kk, fname_cat in enumerate(catalog_fnames):
            base_cat = os.path.basename(fname_cat).lower()
            if base_asc == base_cat:
                exts.append(fname_cat)
                match_cat[kk] = True
                break

        # Load schema
        sch = pas.schema.SchemaDict(schema=fname_asc)
        log.debug("Created schema for '{}'".format(fname_asc))
        # Apply extensions
        for ex in exts:
            sch.extend(ex)
            log.debug("    Extended with '{}'".format(ex))

        schema.append([sch, base_asc])

    # Go through each `catalog` schema-file not already matched, and create schema
    # ------------------------------------------------------------------------------------
    for kk, fname_cat in enumerate(catalog_fnames):
        # Skip if already matched, and used in a schema
        if match_cat[kk]:
            continue

        base_cat = os.path.basename(fname_cat).lower()

        # Load schema
        sch = pas.schema.SchemaDict(schema=fname_cat)
        log.debug("Created schema for '{}'".format(fname_cat))
        schema.append([sch, base_cat])

    # Save compiled schema to output
    # ----------------------------------------
    path_output = AC_PATHS.SCHEMA_OUTPUT
    if not os.path.exists(path_output) or not os.path.isdir(path_output):
        log.raise_error("Path '{}' does not exist or is not a directory!".format(path_output))

    errors = []
    for sch, name in schema:
        base = name.split('.json')[0]
        fname_out = os.path.join(path_output, name)
        log.info("Saving schema '{}' to '{}'".format(base, fname_out))

        sch.dump(fname_out)
        if not os.path.exists(fname_out):
            log.raise_error("Output of schema to '{}' failed!".format(fname_out))

        try:
            sch.validate()
        except pas.ValidationError as err:
            log.error("Validation failure on '{}' ({})!".format(name, fname_out))
            log.error(str(err))
            errors.append(err)

    if len(errors) > 0:
        log.raise_error("Validation of schema failed!")

    return


'''
def setup(log, catalog_info=None):

    # Path to catalog-specific schema-files
    if catalog_info is not None:
        catalog_name = catalog_info['catalog_name']
        catalog_schema_path = catalog_info.get('schema_path', None)
    else:
        catalog_name = None
        catalog_schema_path = None

    # Clear 'metaput' and 'output' directories
    clear_dirs = [AC_PATHS.SCHEMA_METAPUT, AC_PATHS.SCHEMA_OUTPUT]
    for clear in clear_dirs:
        pattern = os.path.join(clear, "", "*.json")
        old_files = sorted(glob.glob(pattern))
        old_files = [os.path.join(clear, fname) for fname in old_files]
        log.info("Deleting {} old schema files in '{}'.".format(len(old_files), clear))
        for fname in old_files:
            os.remove(fname)

    # Copy all `astroschema` schema files to the 'input' directory
    astroschema_fnames = pas.copy_schema_files(AC_PATHS.SCHEMA_INPUT_ASTROSCHEMA, verbose=False)
    log.info("Copied {} 'astroschema' files to input".format(len(astroschema_fnames)))

    # Load names of astrocats files (already in 'input')
    # ------------------------------------------------------------
    pattern = os.path.join(AC_PATHS.SCHEMA_INPUT_ASTROCATS, "*.json")
    astrocats_fnames = list(sorted(glob.glob(pattern)))
    # Add directory component of path to filenames
    pattern_path = os.path.dirname(pattern)
    astrocats_fnames = [os.path.join(pattern_path, fname) for fname in astrocats_fnames]
    log.info("Loaded {} 'astrocats' file names from input".format(len(astrocats_fnames)))

    # Copy catalog-specific files to 'input' directory
    # ----------------------------------------------------------
    if catalog_schema_path is None:
        catalog_fnames = []
    else:
        pattern = os.path.join(catalog_schema_path, "*.json")
        catalog_fnames = sorted(glob.glob(pattern))

    log.info("Loaded {} file names from '{}' catalog input".format(
        len(catalog_fnames), catalog_name))

    # Store all source filenames from input and construct the destination names for 'metaput' dir
    # -------------------------------------------------------------------------------------------

    def get_metaput_fname(src, group):
        base = os.path.basename(src)
        if not base.lower().startswith(group.lower()):
            base = group.lower() + "_" + base
        dst = os.path.join(AC_PATHS.SCHEMA_METAPUT, base)
        return dst

    src_list = []
    dst_list = []
    rename_file_list = []
    fnames = [astroschema_fnames, astrocats_fnames, catalog_fnames]
    groups = ["astroschema", "astrocats", catalog_name]
    log.info("Loading and constructing schema file names")
    for src_files, grp in zip(fnames, groups):
        for src in src_files:
            src_list.append(src)
            dst = get_metaput_fname(src, grp)
            # Make sure there are no duplicates
            if dst in dst_list:
                idx = dst_list.index(dst)
                err = "Destination '{}' already exists for '{}', cannot be used for '{}'!".format(
                    dst, src_list[idx], src)
                log.raise_error(err, ValueError)

            dst_list.append(dst)
            rename_file_list.append([src, dst])

    # Copy all input schema files to output with modified names
    # ----------------------------------------------------------------------
    log.info("Copying {} schema files to metaput".format(len(src_list)))
    for src, dst in zip(src_list, dst_list):
        shutil.copy(src, dst)
        # base = os.path.basename(dst)
        # log.debug("Copied '{}' : '{}' ==> '{}'".format(base, src, dst))
        # Make sure `id` parameter exists; check that it matches filename later (after rename)
        check_schema_id(dst, log, id_equals_name=False)

    # Modify astroschema files to use (e.g.) "astroschema_source.json" instead of "source.json"
    # This is also applied to 'id' attributes which is required for correct reference resolution
    log.info("Replacing astroschema internal references")
    rename_internal_astroschema_refs(rename_file_list, log)

    for dst in dst_list:
        # Make sure that `id` matches the filename
        check_schema_id(dst, log, id_equals_name=True)

    # Try loading and validating all schema files
    for dst in dst_list:
        base = os.path.basename(dst)
        log.debug("Loading and validating '{}'".format(base))
        schema = pas.SchemaDict(dst)
        schema.validate()

    return
'''


def check_schema_id(fname, log, id_equals_name=False):
    """Make sure schema has `id` attribute matching filename.  Required for reference resolution.
    """
    data, *args = pas.utils.load_schema_dict(fname)
    base = os.path.basename(fname)
    # Make sure `id` is present
    if "id" not in data:
        err = "`id` is missing, but required in all schema (must match filename)!"
        log.raise_error(err, ValueError)
    # Also check that the `id` matches the filename
    elif id_equals_name and (data['id'] != base):
        err = "`id` '{}' must match filename '{}'!".format(data['id'], fname)
        log.raise_error(err, ValueError)

    return


def rename_internal_astroschema_refs(astroschema_list, log):
    """Rename internal schema references from old to new filenames.

    This is applied to both the 'id' and 'title' values and properties.

    """
    # Construct a mapping from source to destination filenames (base-name with extension)
    mapping = [[os.path.basename(fn) for fn in fnames] for fnames in astroschema_list]

    for _, fname in astroschema_list:
        if not os.path.exists(fname):
            raise RuntimeError("Schema file '{}' does not exist!".format(fname))

        # base = os.path.basename(fname)
        with open(fname, 'r') as inn:
            schema = inn.read()

        for old, new in mapping:
            # log.debug("In '{}' renaming '{}' ==> '{}'".format(base, old, new))
            # schema = schema.replace(old, new)
            # Use `re` to avoid substrings e.g. 'foo_source.json' ==> 'foo_foo_source.json'
            schema = re.sub(r"\b{}\b".format(old), new, schema)

        with open(fname, 'w') as out:
            out.write(schema)

    return
