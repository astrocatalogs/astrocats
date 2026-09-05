"""Tests for catalog data structures and the original import test task."""
from astrocats.catalog.analysis import Analysis
from astrocats.catalog.catalog import Catalog
from astrocats.catalog.key import KEY_TYPES, Key, KeyCollection
from astrocats.catalog.tasks.test import FAKE_ALIAS_1, do_test


def test_key_collection_and_key_check():
    class SAMPLE(KeyCollection):
        NAME = Key("name", KEY_TYPES.STRING)
        VALUE = Key("value", KEY_TYPES.NUMERIC)

    assert "NAME" in SAMPLE.keys() or "name" in SAMPLE.vals()
    assert SAMPLE.NAME in SAMPLE.vals()
    assert SAMPLE.NAME.check("hello")
    assert not SAMPLE.NAME.check(1)
    assert SAMPLE.VALUE.check("1.2")
    assert not SAMPLE.VALUE.check("not-a-number")


def test_package_import_and_version():
    import astrocats
    from astrocats.catalog.catalog import Catalog as CatalogCls
    from astrocats.catalog.entry import Entry
    from astrocats.catalog.utils import is_number

    assert astrocats.__version__ == "0.5.0"
    assert CatalogCls is Catalog
    assert Entry is not None
    assert is_number("3.14")


def test_original_import_test_task(catalog_offline):
    """Run the historical `do_test` catalog task without network access."""
    do_test(catalog_offline)
    assert FAKE_ALIAS_1 in catalog_offline.entries


def test_import_data_min_priority_test(catalog_offline):
    """Reproduce the Travis `import --min-task-priority test` flow."""
    catalog_offline.import_data()
    assert FAKE_ALIAS_1 in catalog_offline.entries


def test_analyze_count(catalog_offline):
    lysis = Analysis(catalog_offline, catalog_offline.log)
    result = lysis.count()
    assert result["num_tasks"] == 4
    assert result["num_files"] >= 0
