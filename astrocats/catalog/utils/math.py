"""Math utilities.
"""
import numpy as np

__all__ = ["convert_pm_errors_lin_to_log", "str_sig_figs"]


def convert_pm_errors_lin_to_log(val, err_lo, err_hi=None):
    """Convert from a central value and plus/minus error(s) to log-value and log-errors.

    Values can be strings or floats.  If strings, then their significant figures are preserved,
    and strings are returned as output.

    Example
    -------
    >>> convert_pm_errors_lin_to_log('4.1e6', '1e6', '5.1e6')
    >>> '6.61', '0.1', '0.54'

    """

    # Convert to floats if needed
    if isinstance(val, str):
        is_string = True
        # Get significant figures
        v_sf = str_sig_figs(val)
        l_sf = str_sig_figs(err_lo)
        # Convert to floats
        val = np.float(val)
        err_lo = np.float(err_lo)
        if err_hi is not None:
            h_sf = str_sig_figs(err_hi)
            err_hi = np.float(err_hi)

    # Convert to errors in log-space
    lv = np.log10(val)
    elo = err_lo/val/np.log(10.0)
    if err_hi is not None:
        ehi = err_hi/val/np.log(10.0)

    # Convert back to string as needed
    if is_string:
        lv = round_to_str(lv, v_sf)
        elo = round_to_str(elo, l_sf)
        if err_hi is not None:
            ehi = round_to_str(ehi, h_sf)

    if err_hi is not None:
        return lv, elo, ehi

    return lv, elo


def str_sig_figs(val):
    """Return the number of significant figures of the input digit string
    From: https://codereview.stackexchange.com/a/122320
    """
    integral, _, fractional = val.split('e')[0].partition(".")
    if fractional:
        return len((integral + fractional).lstrip('0'))
    else:
        return len(integral.strip('0'))


def round_to_str(val, ndec):
    rounded = np.around(val, ndec)
    str_val = "{xx:.{sf:}f}".format(xx=rounded, sf=ndec)
    return str_val
