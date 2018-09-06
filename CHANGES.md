# AstroCats Release Notes #

[AstroCats on GitHub](https://github.com/astrocatalogs/astrocats)

**AstroCats** maintains a single, primary, active branch: `master`.  We try to follow [Semantic Versioning](http://semver.org/) (`MAJOR.MINOR.PATCH`), where changes to the `master` branch should generally correspond to changes to (at least) the `MINOR` version number.

The ['change-log' (below)](#changelog) in this file should summarize **all API changes** that have effects on individual catalogs.

**Contents**
* [Versions - Change Log](#changelog)
    * [v0.2.0 - 2016/07/18](#v0.2.0)  
      Code restructured into **AstroCats** package, with template classes for subcatalogs.  **Supernovae** catalog divorced into its own repository.

---

### Future / To-Do ###

* FIX THE DAMN ARGPARSE
* Create a settings/arguments/parameters class which handles storing parameters (instead of just an `argparse.Namespace`) and has methods for running `argparse`. 
* Create specific 'tasks' for cloning/pulling input/output repositories.
* Add to the `test` task to include coverage of the basic types (e.g. `spectrum`, etc).
* Add fuzzy-string logic to command-line arguments

---

<a name='changelog'>
## Versions ##


### Current ###

- Remove `utils.em_band.py` dependencies, should be moved to `supernovae`, causes circular import errors
- Created `Error`, `Correlation`, `Model`, `Realization` subclasses of `pyastroschema.struct`
- Move `CatDictError` to `struct.py`
- Added existing schema ('.json') files to repository
- Moved most struct classes into `struct.py` (not `Photometry` and `Entry` which still have independent logic.


<a name='v0.3.42'>
### v0.3.42 - 2018/07/29 ###

- `utils.strings.get_entry_filename` ==> `utils.strings.get_filename`
- `entry.get_filename` ==> `utils.strings.get_filename`
- All usage of *astrocats* (specifically) 'schema' removed for the moment.
- `Entry` changed to be a subclass of astroschema struct instead of `OrderedDict`, using `pyastroschema.Keychain` instead of `KeyCollection`.
- Remove `check()` method, use `validate()` (inherited from pyastroschema) instead
- `Photometry` and `Spectra` (and associated keys) updated to use astroschema
    - Random methods that used to be in `photometry.py` moved to new utilities submodule `em_bands.py`
- `utils.dates` new `astrotime` method to replace separate usaged of `astropy.time.Time`; this needs to be fully integrated still.
- Updated keys for `Testnova` as temporary fixes for key-value incompatibility issues of whether or not underscores are included...

<a name='v0.3.41'>
### v0.3.41 - 2018/07/27 ###


<a name='v0.3.40'>
### v0.3.40 - 2018/07/16 ###


<a name='v0.3.39'>
### v0.3.38 - 2018/06/27 ###
- Introduced the `testcat` submodule which is a stripped-down version of the `supernovae` catalog for testing purposes.  Only the 'internal', 'radio', 'xray', 'cfa_photo', and 'cfa_spectra' tasks have been preserved.  A dedicated output repo was created, 'astrocatalogs/testcat-output' so that comparisons can be madel; currently the version (tag) for that repo is 'v0.0_init-ref'.
- `astrocats/catalog/`
    - `gitter.py`
        - `git_add_commit_push_all_repos()`
            - [BUG] this command was failing when too many files were being added, because the command line-length was too long.  If many files are being added, break them into separate chunks which are `git add`ed separately.


<a name='v0.3.38'>
### v0.3.38 - 2018/06/23 ###

- Added a new `Analysis` class in [astrocats/catalog/analysis.py](https://github.com/astrocatalogs/astrocats/blob/master/astrocats/catalog/analysis.py).
    - Basic 'count' functionality to report the number of files and tasks in each catalog.
- Added subcommands for git repositories in [astrocats/catalog/catalog.py](https://github.com/astrocatalogs/astrocats/blob/master/astrocats/catalog/catalog.py).
    - The `git-push` subcommand can now be used to add, commit and push all data files in each data repository.  This works in all installed catalogs.
    - `git-clone`: clone all data repositories included in the `repos.json` input file.  *NOTE: the name given in the input file must match the github repository name for this to work.  If the directory name should be different from the github url, then it must be added manually.*
    - 'git-reset' :
        - `git-reset-local`: reset the repository to the local HEAD, i.e. it runs `git reset --hard` in each data repository.
        - `git-reset-origin`: reset the repository to 'origin/master', i.e. it runs `git reset --hard origin/master` in each data repository.
    - 'pull' : There is (currently) no pull command.  Instead, repositories must be pulled manually, or alternatively the `git-reset-[]` commands can be used to hard-reset.
    - `git-status` : print the status of the repo in each directory.
- Loading cached/archived URLs, and the 'refresh' parameters.
    - A new method `Catalog.load_url` has been added to take the place of `Catalog.load_cached_url` and `Task.load_archive` (both of which should now be considered deprecated), as well as some general functionality that has been written into individual tasks.
    - The following command-line arguments and corresponding parameters have been deprecated: {`refresh`, `refresh-list`, `refresh-all`}
- `astrocats/catalog/entry.py`
    - `Entry.add_alias` [new-function]
        - New method to add aliases to an existing entry after first 'cleaning' the alias name - in the same way as the entry names are cleaned by the containing catalog.  In this way, the stored aliases should be guaranteed (in general) to match the corresponding entry names (and naming styles).
    - `Entry._get_save_path`
        - Method now *requires* that an output data repository exists for a file to be saved (instead of just saving to the `output/` directory itself).  New catalogs without output data repos will raise a `RuntimeError` from here when trying to save.
- `astrocats/catalog/source.py`
    - `Source.bibcode_from_url` [new-function]
        - Function extracts the Bibcode from an *ADS-URL* if possible.
- `astrocats/catalog/catalog.py`
    - `Catalog.load_cached_url` [DEPRECATED]
        - Replaced by new method `load_url`
    - `Catalog.load_url` [new-function]
        - This method will load text data from either a URL or a pre-cached file of the data.
        - The detailed behavior of the method depends on whether the code is being run in `archived` or `update` mode, and the settings for the particular task calling the method.
- `astrocats/catalog/task.py`
    - `Task.load_archive` [DEPRECATED]
        - Replaced by functionality in `Catalog.load_url`.


<a name='v0.2.0'>
### v0.2.0 - 2016/07/18 ###

- The code has been completely restructured from a single `import.py` script to load and create the supernovae catalog, into two separate packages (and corresponding repositories): [astrocatalogs/AstroCats](https://github.com/astrocatalogs/astrocats), and the supernova-specific [astrocatalogs/Supernovae](https://github.com/astrocatalogs/supernovae).
- This package, `astrocats`, contains the core machinery for creating any general, astronomical catalog.  In particular, the `astrocats/catalog` directory contains the base set of classes which provide all of the desired functionality.
