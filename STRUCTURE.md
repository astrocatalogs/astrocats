# Structure #

Directories: Titlecase,  
Files: lowercase  
Classes: TitleCase  
Variables: lowercase, underscore (for spaces)

-   **Catalogs**
    -   \__init__.py
    -   \__main__.py
    -   main.py
    -   scripts/
    -   **Catalog** 
        -   catalog.py
        -   entry.py
        -   constants.py
        -   task.py
        -   utils/
    -   **SNe**
        -   tasks
        -   supernova.py
        -   input/
            -   tasks.json
            -   internal/
            -   external/
            -   ...
        -   output/
            -   catalog.json
            -   bibauthors.json
            -   ...
        -   scripts/
            -   download...
            -   find...
            -   make...
        -   html/
            -   ...
    -   **BH**
        -   src/
        -   input/
        -   output/
        -   html/
