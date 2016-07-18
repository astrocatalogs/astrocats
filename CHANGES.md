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

<a name='v0.2.0'>
### v0.2.0 - 2016/07/18 ###

- The code has been completely restructured from a single `import.py` script to load and create the supernovae catalog, into two separate packages (and corresponding repositories): [astrocatalogs/AstroCats](https://github.com/astrocatalogs/astrocats), and the supernova-specific [astrocatalogs/Supernovae](https://github.com/astrocatalogs/supernovae).
- This package, `astrocats`, contains the core machinery for creating any general, astronomical catalog.  In particular, the `astrocats/catalog` directory contains the base set of classes which provide all of the desired functionality.
- A basic [AstroCats - wiki](https://github.com/astrocatalogs/astrocats/wiki) has been created with descriptions of the code structure, and how to develop individual catalogs using the **AstroCats** framework.
