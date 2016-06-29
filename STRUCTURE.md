# Structure #

Directories: Titlecase,  
Files: lowercase  
Classes: TitleCase  
Variables: lowercase, underscore (for spaces)

-   **OpenCatalogs**
-   LICENSE
-   README.md
-   requirements.txt
-   STRUCTURE.md
-   TODO.md
-   **astrocats**
    -   \__init__.py
    -   \__main__.py
    -   main.py
    -   **catalog** 
        -   catalog.py
        -   constants.py
        -   entry.py
        -   task.py
        -   utils/
    -   **supernovae**
        -   supernova.py
        -   tasks/
        -   input/
            -   biberrors.json
            -   non-sne-types.json
            -   
            -   repos.json
            -   tasks.json
            -   internal/
            -   external/
            -   ...
        -   output/
            -   catalog.json
            -   cache/  
                -   bibauthors.json
                -   extinctions.json
            -   ...
        -   scripts/
            -   download...
            -   find...
            -   make...
        -   html/
            -   ...
