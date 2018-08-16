# Converting to Astrocats 1.0 (redesign)

## Catalog
- Required class attributes:
    - `RAISE_ERROR_ON_ADDITION_FAILURE`

## Entry
- `Entry.get_filename()` ==> `utils.get_filename()`

## Tasks
- `tasks.json`
    - Make sure the 'module' parameter points to the correct location, i.e. `astrocats.tasks.<TASK>` instead of `catalog.tasks.<TASK>`.
- 
