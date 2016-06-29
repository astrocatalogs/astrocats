'''Utility functions related to dates.
'''
from cdecimal import Decimal

__all__ = ['jd_to_mjd', 'make_date_string']


def jd_to_mjd(jd):
    return jd - Decimal(2400000.5)


def make_date_string(year, month='', day=''):
    if not year:
        raise ValueError(
            "At least the year must be specified when constructing date "
            "string")
    datestring = str(year)
    if month:
        datestring = datestring + '/' + str(month).zfill(2)
    if day:
        datestring = datestring + '/' + str(day).zfill(2)

    return datestring
