"""
"""
import os

from astrocats.catalog.catalog import ENTRY
from astrocats.catalog.source import SOURCE
from astrocats.catalog.quantity import QUANTITY
from astrocats.catalog.photometry import PHOTOMETRY, set_pd_mag_from_counts
from astrocats.catalog.utils import tprint, tq, pbar_strings

FAKE_ALIAS_1 = 'EN-TEST-AA'
FAKE_ALIAS_2 = 'PS-TEST-AB'
FAKE_ALIAS_3 = 'PTF-TEST-BA'
FAKE_ALIAS_4 = 'SN2020abc'
FAKE_ALIAS_5 = 'AT2016omg'

FAKE_NAME_1 = 'Private et al. 2025'
FAKE_BIBCODE_1 = '2025Tst...123..456Z'

FAKE_NAME_2 = 'Jerk et al. 1925'
FAKE_BIBCODE_2 = '1925Tst...987..654A'

FAKE_NAME_3 = 'Ignored et al. 2001'
FAKE_BIBCODE_3 = '2001PASJ..547..364C'

FAKE_REDZ_1 = '1.123'
FAKE_REDZ_2 = '0.987'


def do_test(catalog):
    log = catalog.log
    log.info("do_test()")
    task_str = catalog.get_current_task_str()
    log.info("`task_str`: '{}'".format(task_str))

    if len(catalog.entries) != 0:
        raise RuntimeError("Run test only with empty catalog.")

    # Test URL retrieve functions
    # ---------------------------
    log.info("`args.archived` = '{}', `current_task.archived` = '{}'".format(
        catalog.args.archived, catalog.current_task.archived))

    test_load_url(catalog)

    # Test repo path functions
    # ------------------------
    paths = catalog.PATHS.get_all_repo_folders()
    for path in tq(paths, currenttask='Test tq progress bar.'):
        tprint('Test tprint.')
        log.debug(path)
    paths = catalog.PATHS.get_repo_input_folders()
    for path in pbar_strings(paths, desc='Test pbar_strings progress bar.'):
        log.debug(path)
    boneyard = catalog.PATHS.get_repo_boneyard()
    log.debug(boneyard)

    # Create a Fake Entry, with some Fake Data
    # ----------------------------------------
    _first_event_first_source(catalog)

    log_str = "ADDING SECOND SOURCE"
    log.info("\n\n{}\n{}\n{}\n\n".format("=" * 100, log_str, "=" * 100))

    # Add new Data, from different source, to same fake entry
    # -------------------------------------------------------
    _first_event_second_source(catalog)

    # Make sure output file for this test exists
    outdir, filename = catalog.entries[FAKE_ALIAS_1]._get_save_path()
    save_name = os.path.join(outdir, filename + '.json')
    if not os.path.exists(save_name):
        raise RuntimeError("File not found in '{}'".format(save_name))
    # Delete created test file
    catalog._delete_entry_file(entry_name=FAKE_ALIAS_1)
    # Make sure it was deleted
    if os.path.exists(save_name):
        raise RuntimeError("File not deleted at '{}'".format(save_name))

    # Delete entry in catalog
    del catalog.entries[FAKE_ALIAS_1]
    # Make sure entry was deleted
    if len(catalog.entries) != 0:
        raise RuntimeError("Error deleting test entry!")

    # Add entry back catalog to test later tasks
    _first_event_first_source(catalog)
    _first_event_second_source(catalog)

    # Test some utility functions
    log.debug("Preferred name for 2nd source: " +
              catalog.get_preferred_name(FAKE_ALIAS_2))
    log.debug("Entry exists? " +
              str(catalog.entry_exists(FAKE_ALIAS_2)))
    log.debug("Entry text: " + catalog.entries[FAKE_ALIAS_1].get_entry_text(
        os.path.join(outdir, filename + '.json')))

    # Third source is a duplicate that will be merged
    _first_event_third_source(catalog)

    # Add second event to perform different tests
    _second_event(catalog)

    # Delete name to test name re-addition in sanitize
    for ii, alias in enumerate(
            list(catalog.entries[FAKE_ALIAS_5][ENTRY.ALIAS])):
        if alias[QUANTITY.VALUE] == FAKE_ALIAS_5:
            del catalog.entries[FAKE_ALIAS_1][ENTRY.ALIAS][ii]
            break

    return


def _first_event_first_source(catalog):
    """Try adding a single source, with some data.
    """
    log = catalog.log
    # Add Entry to Catalog
    log.info("Calling: ``add_entry('{}')``".format(FAKE_ALIAS_1))
    name = catalog.add_entry(FAKE_ALIAS_1)
    log.debug("\t `name`: '{}'".format(name))
    log.debug("\n{}\n".format(repr(catalog.entries[name])))
    # Make sure entry exists
    if FAKE_ALIAS_1 not in catalog.entries:
        raise RuntimeError("`FAKE_ALIAS_1`: '{}' is not in entries".format(
            FAKE_ALIAS_1))
    # Make sure entry has the correct name
    stored_name = catalog.entries[FAKE_ALIAS_1][ENTRY.NAME]
    if stored_name != FAKE_ALIAS_1:
        raise RuntimeError("`FAKE_ALIAS_1`[{}]: '{}' does not match".format(
            ENTRY.NAME, stored_name, FAKE_ALIAS_1))

    # Add source to entry
    log.info("Calling: ``add_source('{}')``".format(FAKE_BIBCODE_1))
    source = catalog.entries[name].add_source(
        name=FAKE_NAME_1, bibcode=FAKE_BIBCODE_1)
    log.debug("\t `source`: '{}'".format(source))
    log.debug("\n{}\n".format(repr(catalog.entries[name])))
    # Make sure source alias is correct
    if source != '1':
        raise RuntimeError("Returned `source`: '{}' is wrong.".format(source))
    # Make sure source has the right properties
    check_source_1(catalog, name)

    # Add alias
    log.info("Calling: ``add_quantity('alias', '{}', '{}')``".format(
        FAKE_ALIAS_2, source))
    catalog.entries[name].add_quantity(ENTRY.ALIAS, FAKE_ALIAS_2, source)
    log.debug("\n{}\n".format(repr(catalog.entries[name])))
    # Make sure source alias is correct
    stored_aliases = catalog.entries[name][ENTRY.ALIAS]
    if ((len(stored_aliases) != 1 or
         stored_aliases[0]['value'] != FAKE_ALIAS_2 or
         stored_aliases[0]['source'] != source)):
        raise RuntimeError("Stored alias: '{}' looks wrong.".format(
            stored_aliases[0]))

    log.info("Calling: ``add_quantity('redshift', '{}', '{}')``".format(
        FAKE_REDZ_1, source))
    catalog.entries[name].add_quantity(
        ENTRY.REDSHIFT, FAKE_REDZ_1, source, kind='spectroscopic')
    log.debug("\n{}\n".format(repr(catalog.entries[name])))

    # Add a fake photometric observation
    log.info("Calling: ``add_photometry(...)``")
    catalog.entries[name].add_photometry(
        time='12345', magnitude='20.0', band='g', e_magnitude='0.01',
        telescope='OWELTMT', instrument='UltraCam', observer='I. M. Fake',
        observatory='Mt. Olympus', survey='Zeus Analog Sky Survey',
        source=source)
    log.debug("\n{}\n".format(repr(catalog.entries[name])))

    # Add a fake photometric observation with mangled data
    log.info("Calling: ``add_photometry(...)``")
    catalog.entries[name].add_photometry(
        time='oiasjdqw', magnitude='oihqwr', source=source)
    log.debug("\n{}\n".format(repr(catalog.entries[name])))

    # Add a fake photometric observation without required fields
    log.info("Calling: ``add_photometry(...)``")
    catalog.entries[name].add_photometry(e_magnitude='0.01')
    log.debug("\n{}\n".format(repr(catalog.entries[name])))

    # Add a fake photometric observation by counts
    log.info("Calling: ``add_photometry(...)``")
    photodict = {
        PHOTOMETRY.TIME: '12345',
        PHOTOMETRY.U_TIME: 'MJD',
        PHOTOMETRY.BAND: 'r',
        PHOTOMETRY.SOURCE: source
    }
    set_pd_mag_from_counts(photodict, str(100.0), ec=10.0, zp=30.0)
    catalog.entries[name].add_photometry(**photodict)
    log.debug("\n{}\n".format(repr(catalog.entries[name])))

    log.info("Calling: ``journal_entries()``")
    catalog.journal_entries()
    log.debug("\n{}\n".format(repr(catalog.entries[name])))
    # Make sure the remaining stub looks right
    check_stub(catalog, name)
    return


def _first_event_second_source(catalog):
    log = catalog.log

    log.info("Calling: ``add_entry('{}')``".format(FAKE_ALIAS_2))
    name = catalog.add_entry(FAKE_ALIAS_2)
    log.debug("\t `name`: '{}'".format(name))
    log.debug("\n{}\n".format(repr(catalog.entries[name])))
    # Make sure the proper name is returned (instead of the alias)
    if name != FAKE_ALIAS_1:
        raise RuntimeError("Returned `name`: '{}' does not match '{}'".format(
            name, FAKE_ALIAS_1))
    # Make sure previous data was loaded
    check_source_1(catalog, name)

    log.info("Calling: ``add_source('{}')``".format(FAKE_BIBCODE_2))
    source = catalog.entries[name].add_source(
        name=FAKE_NAME_2, bibcode=FAKE_BIBCODE_2)
    log.debug("\t `source`: '{}'".format(source))
    log.debug("\n{}\n".format(repr(catalog.entries[name])))

    log.info("Calling: ``add_quantity('alias', '{}', '{}')``".format(
        FAKE_ALIAS_3, source))
    catalog.entries[name].add_quantity(ENTRY.ALIAS, FAKE_ALIAS_3, source)
    log.debug("\n{}\n".format(repr(catalog.entries[name])))
    check_source_2(catalog, name)

    log.info("Calling: ``add_quantity('redshift', '{}', '{}')``".format(
        FAKE_REDZ_2, source))
    catalog.entries[name].add_quantity(
        ENTRY.REDSHIFT, FAKE_REDZ_2, source, kind='spectroscopic')
    log.debug("\n{}\n".format(repr(catalog.entries[name])))

    # Add a fake spectral observation
    log.info("Calling: ``add_spectrum(...)``")
    wavelengths = [str(1.0*x) for x in range(1000, 9000, 100)]
    fluxes = wavelengths
    errors = wavelengths
    catalog.entries[name].add_spectrum(
        u_wavelengths='Angstrom', u_fluxes='erg/s/cm^2/Angstrom',
        u_errors='erg/s/cm^2/Angstrom', filename='my_spectrum.txt',
        time='12345', u_time='MJD',
        wavelengths=wavelengths, fluxes=fluxes, errors=errors,
        telescope='OWELTMT', instrument='MOSICE', observer='I. M. Fake',
        observatory='Mt. Everest', survey='Hillary Transient Factory',
        source=source, deredshifted=True)
    log.debug("\n{}\n".format(repr(catalog.entries[name])))

    # Add a duplicate of the above spectrum, shouldn't be added.
    log.info("Calling: ``add_spectrum(...)``")
    wavelengths = [str(1.0*x) for x in range(1000, 9000, 100)]
    fluxes = wavelengths
    errors = wavelengths
    catalog.entries[name].add_spectrum(
        u_wavelengths='Angstrom', u_fluxes='erg/s/cm^2/Angstrom',
        u_errors='erg/s/cm^2/Angstrom', filename='my_spectrum.txt',
        time='12345',
        wavelengths=wavelengths, fluxes=fluxes, source=source)
    log.debug("\n{}\n".format(repr(catalog.entries[name])))

    log.info("Calling: ``journal_entries()``")
    catalog.journal_entries()
    log.debug("\n{}\n".format(repr(catalog.entries[name])))
    check_stub(catalog, name)
    return


def _first_event_third_source(catalog):
    log = catalog.log

    log.info("Calling: ``new_entry('{}')``".format(FAKE_ALIAS_4))
    (name, source) = catalog.new_entry(FAKE_ALIAS_4, srcname=FAKE_NAME_2,
                                       bibcode=FAKE_BIBCODE_2)
    log.debug("\t `name`: '{}', `source`: '{}'".format(name, source))
    log.debug("\n{}\n".format(repr(catalog.entries[name])))

    log.info("Calling: ``add_quantity('alias', '{}', '{}')``".format(
        FAKE_ALIAS_1, source))
    catalog.entries[name].add_quantity(ENTRY.ALIAS, FAKE_ALIAS_1, source)
    log.debug("\n{}\n".format(repr(catalog.entries[name])))

    # Add erroroneous redshift, which should cause add_quantity below to fail.
    log.info("Calling: ``add_error('{}', '{}', '{}')``".format(
        FAKE_BIBCODE_2, SOURCE.BIBCODE, ENTRY.REDSHIFT))
    catalog.entries[name].add_error(
        FAKE_BIBCODE_2, kind=SOURCE.BIBCODE, extra=ENTRY.REDSHIFT)
    log.debug("\n{}\n".format(repr(catalog.entries[name])))

    log.info("Calling: ``add_quantity('redshift', '{}', '{}')``".format(
        FAKE_REDZ_2, source))
    catalog.entries[name].add_quantity(
        ENTRY.REDSHIFT, FAKE_REDZ_2, source, kind='spectroscopic')
    log.debug("\n{}\n".format(repr(catalog.entries[name])))

    return


def _second_event(catalog):
    log = catalog.log

    log.info("Calling: ``new_entry('{}')``".format(FAKE_ALIAS_5))
    (name, source) = catalog.new_entry(FAKE_ALIAS_5, srcname=FAKE_NAME_2,
                                       bibcode=FAKE_BIBCODE_2)
    log.debug("\t `name`: '{}', `source`: '{}'".format(name, source))
    log.debug("\n{}\n".format(repr(catalog.entries[name])))

    # Add an orphan source
    log.info("Calling: ``add_source('{}')``".format(FAKE_BIBCODE_3))
    source = catalog.entries[name].add_source(
        bibcode=FAKE_BIBCODE_3)
    source = catalog.entries[name].add_source(
        name=FAKE_NAME_3)
    log.debug("\t `source`: '{}'".format(source))
    log.debug("\n{}\n".format(repr(catalog.entries[name])))

    return


def check_source_1(catalog, name):
    stored_sources = catalog.entries[name][ENTRY.SOURCES]
    catalog.log.debug(str(stored_sources))
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


def test_load_url(catalog):
    """Test the `Catalog.load_url` methods.
    """
    log = catalog.log
    log.info("Testing `Catalog.load_url`.")

    TEST_URL = "http://www.google.com"       # This should work
    # This shouldnt change much
    TEST_URL_STATIC = "https://www.cfa.harvard.edu/about/contact"
    TEST_BAD_URL = "http://www.BADgooBADgleBAD.com"   # this should fail
    TEST_URL_FNAME_1 = 'google-test-1.txt'   # This should never exist
    TEST_URL_FNAME_2 = 'google-test-2.txt'   # This should be created
    TEST_STATIC_FNAME = 'static-test.txt'

    INPUT_PATH = catalog.get_current_task_repo()

    # Do not Write
    log.info("calling `load_url(...write=False)`")
    test_data = catalog.load_url(TEST_URL, TEST_URL_FNAME_1, write=False)
    if os.path.exists(TEST_URL_FNAME_1):
        err_str = "File '{}' should not have been produced!".format(
            TEST_URL_FNAME_1)
        log_raise(err_str, log)
    if test_data is None:
        err_str = "URL '{}' returned no data!".format(
            TEST_URL)
        log_raise(err_str, log)

    # Do Write
    log.info("calling `load_url(...write=True)`")
    test_path = os.path.join(INPUT_PATH, TEST_URL_FNAME_2)
    # If file exists, delete it to make sure it is recreated
    if os.path.exists(test_path):
        os.remove(test_path)
    test_data = catalog.load_url(TEST_URL, TEST_URL_FNAME_2, write=True)
    if not os.path.exists(test_path):
        err_str = "File '{}' should have been produced!".format(
            test_path)
        log_raise(err_str, log)
    if test_data is None:
        err_str = "URL '{}' returned no data!".format(
            TEST_URL)
        log_raise(err_str, log)

    # Make sure that update mode is working
    log.info("calling `load_url(...update_mode=True)`")
    log.info("WARNING: not testing `update_mode=True`, something is "
             "triggering updates to be detected when none-exist. "
             "Requires further investigation.")
    # test_path = os.path.join(INPUT_PATH, TEST_STATIC_FNAME)
    # # If file exists, delete it to make sure it is recreated
    # if os.path.exists(test_path):
    #     os.remove(test_path)
    # # Create file in update mode
    # test_data = catalog.load_url(TEST_URL_STATIC, TEST_STATIC_FNAME,
    #                              write=True, update_mode=True)
    # if test_data is None:
    #     err_str = "Update mode should Still have created file!"
    #     log_raise(err_str, log)
    # # update mode should return None for unchanged URL
    # test_data = catalog.load_url(TEST_URL_STATIC, TEST_STATIC_FNAME,
    #                              write=False, update_mode=True)
    # if test_data is not None:
    #     err_str = "Update mode should have returned None."
    #     log_raise(err_str, log)

    # Use an invalid URL, but valid file
    log.info("calling `load_url(..., write=False)` with invalid URL")
    test_data = catalog.load_url(TEST_BAD_URL, TEST_URL_FNAME_2, write=False)
    if test_data is None:
        err_str = "Data should still have been retrieved!"
        log_raise(err_str, log)

    # Use invalid URL and invalid file, should return None
    log.info("calling `load_url(..., fail=False)` with invalid URL and File")
    test_data = catalog.load_url(TEST_BAD_URL, TEST_URL_FNAME_1,
                                 write=False, fail=False)
    if test_data is not None:
        err_str = "Return value should have been None!"
        err_str += "\n'{}'".format(test_data)
        log_raise(err_str, log)

    # Use invalid URL and invalid file, with `fail` flag --- should raise
    log.info("calling `load_url(..., fail=True)` with invalid URL and File")
    try:
        test_data = catalog.load_url(TEST_BAD_URL, TEST_URL_FNAME_1,
                                     write=False, fail=True)
    except RuntimeError as err:
        log.debug("Caught error '{}'".format(str(err)))
    # Error should have been raised... so if not, raise error
    else:
        err_str = "An error should have been raised!"
        err_str += "\nInstead, got '{}'".format(test_data)
        log_raise(err_str)

    return

    # Test deprecated function
    test_path = os.path.join(
        catalog.PATHS.get_repo_output_folders()[0], 'test.html')
    catalog.load_cached_url('http://google.com', test_path)

    return


def log_raise(err_str, log):
    log.error(err_str)
    raise RuntimeError(err_str)
