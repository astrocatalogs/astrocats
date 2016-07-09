# -*- coding: utf-8 -*-
"""Set preferred supernova names.
"""


def sanitize(catalog):
    """This simply calls the catalog function that performs this task.
    """
    catalog.sanitize()
    catalog.save_caches()
    return
