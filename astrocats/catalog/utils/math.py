"""Math utilities.
"""
import numpy as np

__all__ = ["convert_lin_to_log", "str_sig_figs", "round_to_str"]


def convert_lin_to_log(val, errors=None, error_interval=False):
    """Convert from a central value and plus/minus error(s) to log-value and log-errors.

    Values can be strings or floats.  If strings, then their significant figures are preserved,
    and strings are returned as output.

    Example
    -------
    >>> convert_pm_errors_lin_to_log('4.1e6', ['1e6', '5.1e6'])
    >>> '6.61', '0.1', '0.54'

    """

    if errors is not None:
        if np.size(errors) != 2:
            raise ValueError("`errors` must be [lo, hi]")

    # Convert to floats if needed
    if isinstance(val, str):
        is_string = True
        # Get significant figures
        v_sf = str_sig_figs(val)
        # Convert to floats
        val = np.float(val)
        if errors is not None:
            e_sf = [str_sig_figs(ee) if (ee is not None) else None
                    for ee in errors]
            errors = [np.float(ee) if (ee is not None) else None
                      for ee in errors]

    # Convert to errors in log-space
    lv = np.log10(val)
    if errors is not None:
        # If these are confidence intervals, convert to plus/minus
        if error_interval:
            def rev(v1, v2):
                return 1.0 if (v1 > v2) else -1.0

            errors = [np.log10(ee) if (ee is not None) else None
                      for ee in errors]
            errors = [(lv - ee)*rev(lv, ee) if (ee is not None) else None
                      for ee in errors]
        # If these errors are plus/minus errors conver to plus/minus log
        else:
            errors = [ee/val/np.log(10.0) if (ee is not None) else None
                      for ee in errors]

    # Convert back to string as needed
    if is_string:
        lv = round_to_str(lv, v_sf)
        if errors is not None:
            errors = [round_to_str(ee, sf) if (ee is not None) else None
                      for (ee, sf) in zip(errors, e_sf)]

    if errors is None:
        return lv

    return lv, errors


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
