"""Unit tests for catalog utility helpers."""
from decimal import Decimal

import pytest

from astrocats.catalog.photometry import (
    PHOTOMETRY,
    set_pd_mag_from_counts,
)
from astrocats.catalog.utils.dates import get_source_year, jd_to_mjd, make_date_string
from astrocats.catalog.utils.digits import (
    get_sig_digits,
    is_integer,
    is_number,
    pretty_num,
    round_sig,
    zpad,
)
from astrocats.catalog.utils.imports import compress_gz, uncompress_gz
from astrocats.catalog.utils.lists import listify
from astrocats.catalog.utils.sorting import alias_priority, bib_priority, repo_priority
from astrocats.catalog.utils.strings import (
    dict_to_pretty_string,
    get_entry_filename,
    single_spaces,
    trim_str_arr,
    uniq_cdl,
)


def test_is_number_and_integer():
    assert is_number("1.23")
    assert is_number(["1", "2.5"])
    assert not is_number("1 2")
    assert not is_number(["1 2"])
    assert not is_number("abc")
    assert is_integer("7")
    assert is_integer(["1", "2"])
    assert not is_integer("1.5")
    assert not is_integer(["1", "x"])


def test_rounding_and_sig_digits():
    assert get_sig_digits("1.230") == 3
    assert get_sig_digits("1.230", strip_zeroes=False) == 4
    assert round_sig(1234.56, 3) == 1230.0
    assert round_sig(0.0, 4) == 0.0
    assert pretty_num(1234.56, 3) == "1230"
    assert zpad("3.2") == "03.2"
    assert zpad("12") == "12"


def test_dates():
    assert make_date_string(2020, 3, 9) == "2020/03/09"
    assert make_date_string(2020) == "2020"
    with pytest.raises(ValueError):
        make_date_string("")
    assert get_source_year({"bibcode": "2017ApJ...835...64G"}) == 2017
    assert get_source_year({"bibcode": "XXXX...."}) == -10000
    with pytest.raises(ValueError):
        get_source_year({})


def test_jd_to_mjd_value():
    assert jd_to_mjd(Decimal("2400000.5")) == Decimal("0.0")


def test_strings_and_lists():
    assert listify("a") == ["a"]
    assert listify(["a", "b"]) == ["a", "b"]
    assert single_spaces("a   b  c") == "a b c"
    assert uniq_cdl(["b", "a", "b"]) == "a,b"
    assert get_entry_filename("SN 1999A/B") == "SN 1999A_B"
    assert trim_str_arr(["1.23456789"], length=4)[0].startswith("1.23")
    pretty = dict_to_pretty_string({"a": 1})
    assert '"a"' in pretty


def test_sorting():
    assert alias_priority("SN1", "SN1") == 0
    assert alias_priority("SN1", "other") == 1
    assert bib_priority({"bibcode": "2010ApJ"})[0] == -2010
    assert repo_priority("sne-2015") == 2015
    assert repo_priority("boneyard") == 1000000000


def test_compress_gz_roundtrip(tmp_path):
    path = tmp_path / "sample.txt"
    path.write_bytes(b"hello astrocats")
    gz_path = compress_gz(str(path))
    assert gz_path.endswith(".gz")
    assert not path.exists()
    restored = uncompress_gz(gz_path)
    assert restored == str(path)
    assert path.read_bytes() == b"hello astrocats"


def test_set_pd_mag_from_counts():
    photodict = {}
    set_pd_mag_from_counts(photodict, str(100.0), ec=10.0, zp=30.0)
    assert PHOTOMETRY.MAGNITUDE in photodict
    assert PHOTOMETRY.E_UPPER_MAGNITUDE in photodict
    assert photodict[PHOTOMETRY.ZERO_POINT] == "30.0"
    mag = float(photodict[PHOTOMETRY.MAGNITUDE])
    assert 24.0 < mag < 26.0
