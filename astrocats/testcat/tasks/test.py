"""
"""

from ..testnova import TESTNOVA


def do_test(catalog):
    '''
    TEST_STR = "    TESTING    "
    TEST_STR = "\n" + "="*len(TEST_STR) + "\n" + TEST_STR + "\n" + "="*len(TEST_STR) + "\n"
    print(TEST_STR)

    # -----------------------------
    name = catalog.add_entry("test")
    print("Added `name` = '{}'".format(name))
    src = catalog.entries[name].add_source(arxivid="1605.01054")
    print("Added `src` = '{}'".format(src))
    catalog.entries[name].add_quantity(TESTNOVA.ALIAS, "trial", src)
    catalog.entries[name].add_quantity(
        TESTNOVA.REDSHIFT, "3.14159", src, e_value=False)
    # print("Added `qnt` = '{}'".format(qnt))

    # print(catalog.entries[name])

    catalog.journal_entries()
    # -----------------------------

    print(TEST_STR)
    '''
    return
