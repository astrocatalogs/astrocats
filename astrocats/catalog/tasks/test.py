"""
"""
import os

from astrocats.catalog.catalog import ENTRY
from astrocats.catalog.source import SOURCE

TEST_NAME = 'EN-TEST-AA'

FAKE_NAME_1 = 'Private et al. 2025'
FAKE_BIBCODE_1 = '2025Tst...123..456Z'
FAKE_ALIAS_1 = 'PS-TEST-AB'
FAKE_REDZ_1 = '1.123'

FAKE_NAME_2 = 'Jerk at al. 1925'
FAKE_BIBCODE_2 = '1925Tst...987..654A'
FAKE_ALIAS_2 = 'PTF-TEST-BA'
FAKE_REDZ_2 = '0.987'


def do_test(catalog):
    log = catalog.log
    log.error("do_test()")
    task_str = catalog.get_current_task_str()
    log.error("`task_str`: '{}'".format(task_str))

    if len(catalog.entries) != 0:
        raise RuntimeError("Run test only with empty catalog.")

    # Create a Fake Entry, with some Fake Data
    # ----------------------------------------
    _first_source(catalog)

    #
    log_str = "ADDING SECOND SOURCE"
    log.error("\n\n{}\n{}\n{}\n\n".format("=" * 100, log_str, "=" * 100))

    # Add new Data, from different source, to same fake entry
    # -------------------------------------------------------
    _second_source(catalog)

    # Make sure output file for this test exists
    outdir, filename = catalog.entries[TEST_NAME]._get_save_path()
    save_name = os.path.join(outdir, filename + '.json')
    if not os.path.exists(save_name):
        raise RuntimeError("File not found in '{}'".format(save_name))
    # Delete created test file
    catalog._delete_entry_file(entry_name=TEST_NAME)
    # Make sure it was deleted
    if os.path.exists(save_name):
        raise RuntimeError("File not deleted at '{}'".format(save_name))

    # Delete entry in catalog
    del catalog.entries[TEST_NAME]
    # Make sure entry was deleted
    if len(catalog.entries) != 0:
        raise RuntimeError("Error deleting test entry!")

    # Add entry back catalog to test later tasks
    _first_source(catalog)
    _second_source(catalog)

    return


def _first_source(catalog):
    """Try adding a single source, with some data.
    """
    log = catalog.log
    # Add Entry to Catalog
    log.error("Calling: ``add_entry('{}')``".format(TEST_NAME))
    name = catalog.add_entry(TEST_NAME)
    log.error("\t `name`: '{}'".format(name))
    log.error("\n{}\n".format(repr(catalog.entries[name])))
    # Make sure entry exists
    if TEST_NAME not in catalog.entries:
        raise RuntimeError("`TEST_NAME`: '{}' is not in entries".format(
            TEST_NAME))
    # Make sure entry has the correct name
    stored_name = catalog.entries[TEST_NAME][ENTRY.NAME]
    if stored_name != TEST_NAME:
        raise RuntimeError("`TEST_NAME`[{}]: '{}' does not match".format(
            ENTRY.NAME, stored_name, TEST_NAME))

    # Add source to entry
    log.error("Calling: ``add_source('{}')``".format(FAKE_BIBCODE_1))
    source = catalog.entries[name].add_source(
        name=FAKE_NAME_1, bibcode=FAKE_BIBCODE_1)
    log.error("\t `source`: '{}'".format(source))
    log.error("\n{}\n".format(repr(catalog.entries[name])))
    # Make sure source alias is correct
    if source != '1':
        raise RuntimeError("Returned `source`: '{}' is wrong.".format(source))
    # Make sure source has the right properties
    check_source_1(catalog, name)

    # Add alias
    log.error("Calling: ``add_quantity('alias', '{}', '{}')``".format(
        FAKE_ALIAS_1, source))
    catalog.entries[name].add_quantity(ENTRY.ALIAS, FAKE_ALIAS_1, source)
    log.error("\n{}\n".format(repr(catalog.entries[name])))
    # Make sure source alias is correct
    stored_aliases = catalog.entries[name][ENTRY.ALIAS]
    if ((len(stored_aliases) != 1 or
         stored_aliases[0]['value'] != FAKE_ALIAS_1 or
         stored_aliases[0]['source'] != source)):
        raise RuntimeError("Stored alias: '{}' looks wrong.".format(
            stored_aliases[0]))

    log.error("Calling: ``add_quantity('redshift', '{}', '{}')``".format(
        FAKE_REDZ_1, source))
    catalog.entries[name].add_quantity(
        ENTRY.REDSHIFT, FAKE_REDZ_1, source, kind='spectroscopic')
    log.error("\n{}\n".format(repr(catalog.entries[name])))

    # Add a fake photometric observation
    log.error("Calling: ``add_photometry(...)``")
    catalog.entries[name].add_photometry(
        time='12345', magnitude='20.0', band='g', e_magnitude='0.01',
        telescope='OWELTMT', instrument='UltraCam', observer='I. M. Fake',
        observatory='Mt. Olympus', survey='Zeus Analog Sky Survey',
        source=source)
    log.error("\n{}\n".format(repr(catalog.entries[name])))

    log.error("Calling: ``journal_entries()``")
    catalog.journal_entries()
    log.error("\n{}\n".format(repr(catalog.entries[name])))
    # Make sure the remaining stub looks right
    check_stub(catalog, name)
    return


def _second_source(catalog):
    log = catalog.log

    log.error("Calling: ``add_entry('{}')``".format(FAKE_ALIAS_1))
    name = catalog.add_entry(FAKE_ALIAS_1)
    log.error("\t `name`: '{}'".format(name))
    log.error("\n{}\n".format(repr(catalog.entries[name])))
    # Make sure the proper name is returned (instead of the alias)
    if name != TEST_NAME:
        raise RuntimeError("Returned `name`: '{}' does not match '{}'".format(
            name, TEST_NAME))
    # Make sure previous data was loaded
    check_source_1(catalog, name)

    log.error("Calling: ``add_source('{}')``".format(FAKE_BIBCODE_2))
    source = catalog.entries[name].add_source(
        name=FAKE_NAME_2, bibcode=FAKE_BIBCODE_2)
    log.error("\t `source`: '{}'".format(source))
    log.error("\n{}\n".format(repr(catalog.entries[name])))

    log.error("Calling: ``add_quantity('alias', '{}', '{}')``".format(
        FAKE_ALIAS_2, source))
    catalog.entries[name].add_quantity(ENTRY.ALIAS, FAKE_ALIAS_2, source)
    log.error("\n{}\n".format(repr(catalog.entries[name])))
    check_source_2(catalog, name)

    log.error("Calling: ``add_quantity('redshift', '{}', '{}')``".format(
        FAKE_REDZ_2, source))
    catalog.entries[name].add_quantity(
        ENTRY.REDSHIFT, FAKE_REDZ_2, source, kind='spectroscopic')
    log.error("\n{}\n".format(repr(catalog.entries[name])))

    # Add a fake spectral observation
    log.error("Calling: ``add_spectrum(...)``")
    wavelengths = [str(1.0*x) for x in range(1000, 9000, 100)]
    fluxes = wavelengths
    errors = wavelengths
    catalog.entries[name].add_spectrum(
        u_wavelengths='Angstrom', u_fluxes='erg/s/cm^2/Angstrom',
        u_errors='erg/s/cm^2/Angstrom',
        time='12345',
        wavelengths=wavelengths, fluxes=fluxes, errors=errors,
        telescope='OWELTMT', instrument='MOSICE', observer='I. M. Fake',
        observatory='Mt. Everest', survey='Hillary Transient Factory',
        source=source)
    log.error("\n{}\n".format(repr(catalog.entries[name])))

    log.error("Calling: ``journal_entries()``")
    catalog.journal_entries()
    log.error("\n{}\n".format(repr(catalog.entries[name])))
    check_stub(catalog, name)
    return


def check_source_1(catalog, name):
    stored_sources = catalog.entries[name][ENTRY.SOURCES]
    print(stored_sources)
    if ((len(stored_sources) != 1 or
         stored_sources[0][SOURCE.NAME] != FAKE_NAME_1 or
         stored_sources[0][SOURCE.BIBCODE] != FAKE_BIBCODE_1)):
        raise RuntimeError("Stored source: '{}' looks wrong.".format(
            stored_sources[0]))
    return


def check_source_2(catalog, name):
    stored_sources = catalog.entries[name][ENTRY.SOURCES]
    names = [src[SOURCE.NAME] for src in stored_sources]
    codes = [src[SOURCE.BIBCODE] for src in stored_sources]
    if ((len(stored_sources) != 2 or
         FAKE_NAME_1 not in names or FAKE_NAME_2 not in names or
         FAKE_BIBCODE_1 not in codes or FAKE_BIBCODE_2 not in codes)):
        raise RuntimeError("Stored sources: '{}' look wrong.".format(
            stored_sources))
    return


def check_stub(catalog, name):
    if not catalog.entries[name]._stub:
        raise RuntimeError("Remaining entry is not a stub.")
    if ENTRY.ALIAS not in catalog.entries[name]:
        raise RuntimeError("Remaining entry is missing '{}'.".format(
            ENTRY.ALIAS))
    if ENTRY.SOURCES in catalog.entries[name]:
        raise RuntimeError("Remaining still has '{}'.".format(
            ENTRY.SOURCES))
    return
