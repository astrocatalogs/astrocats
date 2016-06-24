# Change Log #

## To-Do ##
-   `importer/mtasks/vizier.py:do_vizier`
    -   Break down into separate functions, or classes.  Can use significant templating to cleanup.
-   `importer/mtasks/donations.py:do_donations`
    -   Break down into separate " " " ".
-   `main.py`
    -   Look at '--archived' vs. '--refresh[list]' and '--full-refresh'... what exactly is desired?
-   `name_clean`
    -   FIX: IMPROVE THIS!
-   `EVENTS.add_quantity`
    -   FIX: IMPROVE THIS!
-   Create methods for parsing/cleaning RA & DEC  (e.g. to go into `add_quantity`)
-   `scripts/constants.py`
    -   Combine `OSC_` values into a class
-   Work on fuzzy finding for events, e.g. for 'ASASSN-13ax' perhaps include alias '13ax'?
-   `write_all_events`: save `non_sne_types` data somewhere instead of reloading each time
    -   Move all of the stuff inside the names loop into the EVENT class.
-   Change 'writevents' from a task to an argument parameter?  (e.g. for `journal_events`)
-   Combine `EVENT.check` with `Events.clean_event`
-   `load_cached_url` add warnings for failures
-   Have the different `add_` methods accept lists and add each
    
    
## Questions ##
-   What is `import_funcs.clear_events` doing?
-   `mtasks.asiago.do_asiago_spectra` `mjd` and `epochstr` arent doing anything.
-   Maybe unify things in 'csv' format, e.g. `do_snls` which seem to start the same way?
-   `mtasks.general_data.do_snls` see note:
    -   "NOTE: Datafiles avail for download suggest diff zeropoints than 30, need to inquire."
-   `clean_event`
    -   Why do the sources need to be 'rebuilt', what exactly is being cleaned here?
        -   i.e. why are all of the other attributes (besides 'url', 'bibcode', and 'name') being deleted?
        -   Is 'url' required if no 'bibcode'?
-   `is_erroneous`
    -   What is happening here?
-   What are stubs?  Why?


## Current ##

-   Drastic restructuring of code, primarily the main `import.py` script into modular, functional, and object oriented code.
