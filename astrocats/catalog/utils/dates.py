'''Utility functions related to dates.
'''
from .digits import is_number
from decimal import Decimal
import warnings

import astropy as ap
import astropy.time  # noqa
from astropy._erfa.core import ErfaWarning

__all__ = ['astrotime', 'jd_to_mjd', 'make_date_string', 'get_source_year', 'mjd_to_datetime']


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


def get_source_year(source):
    if 'bibcode' in source:
        if is_number(source['bibcode'][:4]):
            return int(source['bibcode'][:4])
        else:
            return -10000
    raise ValueError('No bibcode available for source!')


def mjd_to_datetime(mjd):
    """Convert from MJD to a datetime object using `astropy.time.Time`.

    Suppresses `ErfaWarning` from astropy.

    Arguments
    ---------
    mjd : float

    Returns
    -------
    dt : `datetime.datetime`

    """
    '''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ErfaWarning)
        time = ap.time.Time(mjd, format='mjd').datetime
    '''
    time = astrotime(mjd, input='mjd', output='datetime')
    return time


def astrotime(datetime, input='mjd', output='datetime', to_str=False):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ErfaWarning)
        time = ap.time.Time(datetime, format=input)

    # If `output` is None, return `astropy.time.Time` instance
    if output is None:
        pass
    else:
        try:
            time = getattr(time, output)
        except AttributeError as err:
            add_err = "Unrecognized `output` type = '{}'!".format(output)
            raise AttributeError(add_err + str(err)) from err

    if to_str:
        time = str(time)

    return time
