"""Class for representing photometric data.
"""
from collections import OrderedDict
# from .entry import KEYS


class PHOTOMETRY:
    TIME = 'time'
    U_TIME = 'u_time'
    E_TIME = 'e_time'
    TELESCOPE = 'telescope'
    INSTRUMENT = 'instrument'
    BAND = 'band'
    MAGNITUDE = 'magnitude'
    E_MAGNITUDE = 'e_magnitude'
    SOURCE = 'source'
    SYSTEM = 'system'
    SCORRECTED = 'scorrected'
    OBSERVATORY = 'observatory'
    OBSERVER = 'observer'
    SURVEY = 'survey'
    KCORRECTED = 'kcorrected'
    FLUX = 'flux'
    FLUX_DENSITY = 'fluxdensity'
    E_FLUX = 'e_flux'
    E_FLUX_DENSITY = 'e_fluxdensity'
    U_FLUX = 'u_flux'
    U_FLUX_DENSITY = 'u_fluxdensity'
    FREQUENCY = 'frequency'
    U_FREQUENCY = 'u_frequency'
    COUNTS = 'counts'
    E_COUNTS = 'e_counts'
    NHMW = 'nhmw'
    PHOTON_INDEX = 'photonindex'
    UNABSORBED_FLUX = 'unabsorbedflux'
    E_UNABSORBED_FLUX = 'e_unabsorbedflux'
    ENERGY = 'energy'
    U_ENERGY = 'u_energy'
    E_LOWER_MAGNITUDE = 'e_lower_magnitude'
    E_UPPER_MAGNITUDE = 'e_upper_magnitude'
    E_LOWER_TIME = 'e_lower_time'
    E_UPPER_TIME = 'e_upper_time'
    MCORRECTED = 'mcorrected'
    UPPERLIMIT = 'upperlimit'
    HOST = 'host'
    INCLUDESHOST = 'includeshost'

    _keys = sorted([kk for kk in dir() if not kk.startswith('_')])


class Photometry(OrderedDict):
    """
    """

    def __init__(self, **kwargs):
        self._check_kwargs(kwargs)

        photo_keys = PHOTOMETRY._keys
        # Iterate over all passed keys, sanitize and store them
        for key, val in kwargs.items():
            # Make sure the key is one of the allowed keys from `PHOTOMETRY`
            if key not in photo_keys:
                raise KeyError("'{}' is not a valid `PHOTOMETRY` key".format(
                    key))

    def _check_kwargs(self, kwargs):
        REQ_ANY = [
            [PHOTOMETRY.TIME, PHOTOMETRY.HOST],
            [PHOTOMETRY.MAGNITUDE, PHOTOMETRY.FLUX_DENSITY, PHOTOMETRY.FLUX,
             PHOTOMETRY.COUNTS, PHOTOMETRY.UNABSORBED_FLUX]]

        for req_any in REQ_ANY:
            if not any([req_key in kwargs for req_key in req_any]):
                err_str = "Require one of: " + ",".join(
                    "'{}'".format(rk) for rk in req_any)
                raise ValueError(err_str)

        if ((not host and not is_number(time)) or
            (not is_number(magnitude) and not is_number(flux) and not
             is_number(fluxdensity) and not is_number(counts))):
            warnings.warn('Time or brightness not numerical, not adding.')
            tprint('Name : "' + name + '", Time: "' + time + '", Band: "' +
                   band + '", AB magnitude: "' + magnitude + '"')
            return

        if (((e_magnitude and not is_number(e_magnitude)) or
             (e_flux and not is_number(e_flux)) or
             (e_fluxdensity and not is_number(e_fluxdensity)) or
             (e_counts and not is_number(e_counts)))):
            warnings.warn('Brightness error not numerical, not adding.')
            tprint('Name : "' + name + '", Time: "' + time +
                   '", Band: "' + band + '", AB error: "' + e_magnitude + '"')
            return

        if e_time and not is_number(e_time):
            warnings.warn('Time error not numerical, not adding.')
            tprint('Name : "' + name + '", Time: "' +
                   time + '", Time error: "' + e_time + '"')
            return

        if ((flux or fluxdensity) and ((not u_flux and not u_fluxdensity) or
                                       (not frequency and not band and not
                                        energy))):
            warnings.warn(
                "Unit and band/frequency must be set when adding photometry "
                "by flux or flux density, not adding.")
            tprint('Name : "' + name + '", Time: "' + time)
            return

        if not source:
            ValueError('Photometry must have source before being added!')

        if self.is_erroneous('photometry', source):
            return

        return
