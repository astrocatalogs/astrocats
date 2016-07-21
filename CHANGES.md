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

* Create a settings/arguments/parameters class which handles storing parameters (instead of just an `argparse.Namespace`) and has methods for running `argparse`. 
* Create specific 'tasks' for cloning/pulling input/output repositories.
* Add to the `test` task to include coverage of the basic types (e.g. `spectrum`, etc).
* Add fuzzy-string logic to command-line arguments

---

<a name='changelog'>
## Versions ##

### Current ###

- Added a new `Analysis` class in [astrocats/catalog/analysis.py](https://github.com/astrocatalogs/astrocats/blob/master/astrocats/catalog/analysis.py).
    - Basic 'count' functionality to report the number of files and tasks in each catalog.
- Added subcommands for git repositories in [astrocats/catalog/catalog.py](https://github.com/astrocatalogs/astrocats/blob/master/astrocats/catalog/catalog.py).
    - The 'push' subcommand can now be used to add, commit and push all data files in each data repository.  This works in all installed catalogs.
- `astrocats/catalog/entry.py`
    - `Entry.add_alias` [new-function]
        - New method to add aliases to an existing entry after first 'cleaning' the alias name - in the same way as the entry names are cleaned by the containing catalog.  In this way, the stored aliases should be guaranteed (in general) to match the corresponding entry names (and naming styles).
- `astrocats/catalog/source.py`
    - `Source.bibcode_from_url` [new-function]
        - Function extracts the Bibcode from an *ADS-URL* if possible.

<a name='v0.2.0'>
### v0.2.0 - 2016/07/18 ###

- The code has been completely restructured from a single `import.py` script to load and create the supernovae catalog, into two separate packages (and corresponding repositories): [astrocatalogs/AstroCats](https://github.com/astrocatalogs/astrocats), and the supernova-specific [astrocatalogs/Supernovae](https://github.com/astrocatalogs/supernovae).
- This package, `astrocats`, contains the core machinery for creating any general, astronomical catalog.  In particular, the `astrocats/catalog` directory contains the base set of classes which provide all of the desired functionality.
