# Scripts README #

This directory contains the scripts for loading, building and processing the Open Supernovae Catalog (OSC) database.

The primary script if `import.py` which directs the import (download) of all OCS data.


## Questions ##
-   What is `import_funcs.clear_events` doing?
-   `mtasks.asiago.do_asiago_spectra` `mjd` and `epochstr` arent doing anything.
-   Maybe unify things in 'csv' format, e.g. `do_snls` which seem to start the same way?
-   `mtasks.general_data.do_snls` see note:
    -    "NOTE: Datafiles avail for download suggest diff zeropoints than 30, need to inquire."
