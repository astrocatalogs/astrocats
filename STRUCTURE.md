# Structure and Code Style #

Unless otherwise noted, everything should be
[PEP8](https://www.python.org/dev/peps/pep-0008) compliant, with 79 characters
maximum per line (we ignore the 72 character docstring maximum).  Please make
variable names descriptive, and add comments describing *everything*.  Avoid
the use of literals whenever possible - using/defining/documenting constants
appropriately.

### Naming Conventions ###
-   Directories: `lowercase`,  
-   Files: `lowercase`  
-   Classes: `TitleCase`  
-   Variables: `lower_case` (i.e. underscores for breaks/spaces)
    -   Constants: `UPPER_CASE`


## File and Directory Structure ##
-   **OpenCatalogs**
-   `LICENSE`
-   `README.md`
-   `requirements.txt`
-   `STRUCTURE.md`
-   `TODO.md`
-   **astrocats**:
    -   `__init__.py`
    -   `__main__.py`
    -   `main.py`:  *the primary entry point for all catalog operations*
    -   **catalog**: *templates and basic machinery for all individual catalogs*
        -   `catalog.py`: *the primary class which handles catalogs*
        -   `entry.py`: *the base class for individual 'entries' in each catalog*
        -   `task.py`: *object to represent the underlying list of operations to construct a catalog (see also `tasks.json`)*
        -   `utils/`: *general purpose utility functions*
    -   **supernovae**
        -   `supernova.py`: *subclass of `entry` specific for the supernova catalog*
        -   `tasks/`: *files associated with operations to build this catalog.*
        -   `input/`: *data files and repositories contributing data/parameters*
            -   `biberrors.json`: *Specific bibcodes with errors*
            -   `non-sne-types.json`: *list of types which are *not* supernovae*
            -   `repos.json`: *list of git repositories in the supernovae catalog*
            -   `tasks.json`: *list of tasks and associated parameters to construct the catalog*
            -   `sne-internal/`: *Hand-made supernovae data files to be loaded and processed*
            -   `sne-external.../`: *Downloaded/externally retrieved data files to be loaded and processed*
            -   ...
        -   `output/`: *data products constructed by the supernova catalog*
            -   `catalog.json`: *output of all supernova data*
            -   `cache/`: *intermediate files used during catalog creation*
                -   `bibauthors.json`
                -   `extinctions.json`
        -   `scripts/`: *individual catalog operations scripts()*
            -   download...
            -   find...
            -   make...
        -   `html/`: *files associated with the supernova catalog webpage*
            -   ...
