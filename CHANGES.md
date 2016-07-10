# AstroCats Release Notes #

[AstroCats on GitHub](https://github.com/astrocatalogs/astrocats)

**AstroCats** maintains two primary active branches: `master` for working, production code; and `dev`, for unstable/changing development code.  Along the `master` branch we try to follow [Semantic Versioning](http://semver.org/) (`MAJOR.MINOR.PATCH`), where changes to the `master` branch should generally correspond to changes to (at least) the `MINOR` version number.

---

### Future / To-Do ###

* Create a settings/arguments/parameters class which handles storing parameters (instead of just an `argparse.Namespace`) and has methods for running `argparse`. 
* Create specific 'tasks' for cloning/pulling input/output repositories.
* Add to the `test` task to include coverage of the basic types (e.g. `spectrum`, etc).
* Add fuzzy-string logic to command-line arguments

---

## Versions ##

### Current ###

- The code has been completely restructured from a single `import.py` script to load and create the supernovae catalog, into two separate packages (and corresponding repositories): [astrocatalogs/AstroCats]( https://github.com/astrocatalogs/astrocats), and the supernova-specific [astrocatalogs/Supernovae]( https://github.com/astrocatalogs/supernovae).
- This package, `astrocats`, contains the core machinery for creating any general, astronomical catalog.
