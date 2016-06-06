"""Utility functions for OSC import.
"""

import os
import json
import warnings
from cdecimal import Decimal
from collections import OrderedDict
from math import log10, floor, sqrt  # , isnan, ceil
from astropy.time import Time as astrotime

from scripts import _FILENAME_SOURCE_SYNONYMS, _FILENAME_TYPE_SYNONYMS, _FILENAME_ATELS, \
    _FILENAME_CBETS, _FILENAME_IAUCS, _FILENAME_EXTINCT, _FILENAME_NON_SNE_TYPES, \
    _FILENAME_BIBAUTHORS
from .. utils import repo_file_list, is_number, get_repo_folders, get_event_filename, tprint, \
    bandrepf, bandmetaf, get_sig_digits, zpad, pretty_num, round_sig, get_repo_years, tq

from . constants import REPR_BETTER_QUANTITY, OSC_BIBCODE, OSC_NAME, OSC_URL, CLIGHT, PREF_KINDS, \
    KM


def add_event(tasks, args, events, name, load=True, delete=True, source='', loadifempty=True):
    if loadifempty and args.update and not len(events):
        load_stubs(tasks, args, events)

    newname = name_clean(name)
    if newname not in events or 'stub' in events[newname]:
        match = ''
        if newname not in events:
            for event in events:
                aliases = get_aliases(events, event)
                if (len(aliases) > 1 and (newname in aliases) and ('distinctfrom' not in events[event] or newname not in events[event]['distinctfrom'])):
                    match = event
                    break
            # FIX: is this supposed to be here??
            if match:
                newname = match

        if load:
            loadedname = load_event_from_file(events, args, tasks, name=newname, delete=delete)
            if loadedname:
                if 'stub' in events[loadedname]:
                    raise(ValueError('Failed to find event file for stubbed event'))
                return loadedname

        if match:
            return match

        events[newname] = OrderedDict()
        events[newname]['name'] = newname
        if source:
            add_quantity(events, newname, 'alias', newname, source)
        if args.verbose and 'stub' not in events[newname]:
            tprint('Added new event ' + newname)
        return newname
    else:
        return newname


def add_photometry(events, name, time="", u_time="MJD", e_time="", telescope="", instrument="", band="",
                   magnitude="", e_magnitude="", source="", upperlimit=False, system="",
                   observatory="", observer="", host=False, includeshost=False, survey="",
                   flux="", fluxdensity="", e_flux="", e_fluxdensity="", u_flux="", u_fluxdensity="", frequency="",
                   u_frequency="", counts="", e_counts="", nhmw="", photonindex="", unabsorbedflux="",
                   e_unabsorbedflux="", energy="", u_energy="", e_lower_magnitude="", e_upper_magnitude=""):
    if (not time and not host) or (not magnitude and not flux and not fluxdensity and not counts and not unabsorbedflux):
        warnings.warn('Time or brightness not specified when adding photometry, not adding.')
        tprint('Name : "' + name + '", Time: "' + time + '", Band: "' + band + '", AB magnitude: "' + magnitude + '"')
        return

    if (not host and not is_number(time)) or (not is_number(magnitude) and not is_number(flux) and not is_number(fluxdensity) and not is_number(counts)):
        warnings.warn('Time or brightness not numerical, not adding.')
        tprint('Name : "' + name + '", Time: "' + time + '", Band: "' + band + '", AB magnitude: "' + magnitude + '"')
        return

    if (((e_magnitude and not is_number(e_magnitude)) or
         (e_flux and not is_number(e_flux)) or
         (e_fluxdensity and not is_number(e_fluxdensity)) or
         (e_counts and not is_number(e_counts)))):
        warnings.warn('Brightness error not numerical, not adding.')
        tprint('Name : "' + name + '", Time: "' + time + '", Band: "' + band + '", AB error: "' + e_magnitude + '"')
        return

    if e_time and not is_number(e_time):
        warnings.warn('Time error not numerical, not adding.')
        tprint('Name : "' + name + '", Time: "' + time + '", Time error: "' + e_time + '"')
        return

    if (flux or fluxdensity) and ((not u_flux and not u_fluxdensity) or (not frequency and not band and not energy)):
        warnings.warn('Unit and band/frequency must be set when adding photometry by flux or flux density, not adding.')
        tprint('Name : "' + name + '", Time: "' + time)
        return

    if not source:
        ValueError('Photometry must have source before being added!')

    if is_erroneous(events, name, 'photometry', source):
        return

    # Do some basic homogenization
    sband = bandrepf(band)

    sinstrument = instrument
    ssystem = system
    stelescope = telescope

    if not sinstrument:
        sinstrument = bandmetaf(sband, 'instrument')
    if not stelescope:
        stelescope = bandmetaf(sband, 'telescope')
    if not ssystem:
        ssystem = bandmetaf(sband, 'system')

    # Look for duplicate data and don't add if duplicate
    if 'photometry' in events[name]:
        for photo in events[name]['photometry']:
            if ((same_tag_str(photo, sband, 'band') and
                 same_tag_str(photo, u_time, 'u_time') and
                 same_tag_num(photo, time, 'time', canbelist=True) and
                 same_tag_num(photo, magnitude, 'magnitude') and
                 (('host' not in photo and not host) or ('host' in photo and host)) and
                 same_tag_num(photo, flux, 'flux') and
                 same_tag_num(photo, unabsorbedflux, 'unabsorbedflux') and
                 same_tag_num(photo, fluxdensity, 'fluxdensity') and
                 same_tag_num(photo, counts, 'counts') and
                 same_tag_num(photo, energy, 'energy') and
                 same_tag_num(photo, frequency, 'frequency') and
                 same_tag_num(photo, photonindex, 'photonindex') and
                 same_tag_num(photo, e_magnitude, 'e_magnitude') and
                 same_tag_num(photo, e_lower_magnitude, 'e_lower_magnitude') and
                 same_tag_num(photo, e_upper_magnitude, 'e_upper_magnitude') and
                 same_tag_num(photo, e_flux, 'e_flux') and
                 same_tag_num(photo, e_unabsorbedflux, 'e_unabsorbedflux') and
                 same_tag_num(photo, e_fluxdensity, 'e_fluxdensity') and
                 same_tag_num(photo, e_counts, 'e_counts') and
                 same_tag_num(photo, u_flux, 'u_flux') and
                 same_tag_num(photo, u_fluxdensity, 'u_fluxdensity') and
                 same_tag_num(photo, u_frequency, 'u_frequency') and
                 same_tag_num(photo, u_energy, 'u_energy') and
                 same_tag_str(photo, ssystem, 'system'))):
                    return

    photoentry = OrderedDict()
    if time:
        photoentry['time'] = time if isinstance(time, list) or isinstance(time, str) else str(time)
    if e_time:
        photoentry['e_time'] = str(e_time)
    if u_time:
        photoentry['u_time'] = u_time
    if sband:
        photoentry['band'] = sband
    if ssystem:
        photoentry['system'] = ssystem
    if magnitude:
        photoentry['magnitude'] = str(magnitude)
    if e_magnitude:
        photoentry['e_magnitude'] = str(e_magnitude)
    if e_lower_magnitude:
        photoentry['e_lower_magnitude'] = str(e_lower_magnitude)
    if e_upper_magnitude:
        photoentry['e_upper_magnitude'] = str(e_upper_magnitude)
    if frequency:
        photoentry['frequency'] = frequency if isinstance(frequency, list) or isinstance(frequency, str) else str(frequency)
    if u_frequency:
        photoentry['u_frequency'] = u_frequency
    if energy:
        photoentry['energy'] = energy if isinstance(energy, list) or isinstance(energy, str) else str(energy)
    if u_energy:
        photoentry['u_energy'] = u_energy
    if flux:
        photoentry['flux'] = str(flux)
    if e_flux:
        photoentry['e_flux'] = str(e_flux)
    if unabsorbedflux:
        photoentry['unabsorbedflux'] = str(unabsorbedflux)
    if e_unabsorbedflux:
        photoentry['e_unabsorbedflux'] = str(e_unabsorbedflux)
    if u_flux:
        photoentry['u_flux'] = str(u_flux)
    if photonindex:
        photoentry['photonindex'] = str(photonindex)
    if fluxdensity:
        photoentry['fluxdensity'] = str(fluxdensity)
    if e_fluxdensity:
        photoentry['e_fluxdensity'] = str(e_fluxdensity)
    if u_fluxdensity:
        photoentry['u_fluxdensity'] = str(u_fluxdensity)
    if counts:
        photoentry['counts'] = str(counts)
    if e_counts:
        photoentry['e_counts'] = str(e_counts)
    if upperlimit:
        photoentry['upperlimit'] = upperlimit
    if host:
        photoentry['host'] = host
    if includeshost:
        photoentry['includeshost'] = includeshost
    if observer:
        photoentry['observer'] = observer
    if survey:
        photoentry['survey'] = survey
    if observatory:
        photoentry['observatory'] = observatory
    if stelescope:
        photoentry['telescope'] = stelescope
    if sinstrument:
        photoentry['instrument'] = sinstrument
    if nhmw:
        photoentry['nhmw'] = nhmw
    if source:
        photoentry['source'] = source
    events[name].setdefault('photometry', []).append(photoentry)


def add_quantity(events, name, quantity, value, sources, forcereplacebetter=False,
                 lowerlimit='', upperlimit='', error='', unit='', kind='', extra=''):
    if not quantity:
        raise(ValueError('Quantity must be specified for add_quantity.'))
    if not sources:
        raise(ValueError('Source must be specified for quantity before it is added.'))
    if not isinstance(value, str) and (not isinstance(value, list) or not isinstance(value[0], str)):
        raise(ValueError('Quantity must be a string or an array of strings.'))

    with open(_FILENAME_TYPE_SYNONYMS, 'r') as f:
        typereps = json.loads(f.read(), object_pairs_hook=OrderedDict)

    if is_erroneous(events, name, quantity, sources):
        return

    svalue = value.strip()
    serror = error.strip()
    skind = kind.strip()
    sunit = ''

    if not svalue or svalue == '--' or svalue == '-':
        return
    if serror and (not is_number(serror) or float(serror) < 0.):
        raise(ValueError('Quanta error value must be a number and positive.'))

    # Set default units
    if not unit and quantity == 'velocity':
        unit = 'KM/s'
    if not unit and quantity == 'ra':
        unit = 'hours'
    if not unit and quantity == 'dec':
        unit = 'degrees'
    if not unit and quantity in ['lumdist', 'comovingdist']:
        unit = 'Mpc'

    # Handle certain quantity
    if quantity == 'alias':
        svalue = name_clean(svalue)
        if 'distinctfrom' in events[name]:
            if svalue in [x['value'] for x in events[name]['distinctfrom']]:
                return
    if quantity in ['velocity', 'redshift', 'ebv', 'lumdist', 'comovingdist']:
        if not is_number(svalue):
            return
    if quantity == 'host':
        if is_number(svalue):
            return
        if svalue.lower() in ['anonymous', 'anon.', 'anon', 'intergalactic']:
            return
        if svalue.startswith('M ') and is_number(svalue[2:]):
            svalue.replace('M ', 'M', 1)
        svalue = svalue.strip("()").replace('  ', ' ', 1)
        svalue = svalue.replace("Abell", "Abell ", 1)
        svalue = svalue.replace("APMUKS(BJ)", "APMUKS(BJ) ", 1)
        svalue = svalue.replace("ARP", "ARP ", 1)
        svalue = svalue.replace("CGCG", "CGCG ", 1)
        svalue = svalue.replace("HOLM", "HOLM ", 1)
        svalue = svalue.replace("IC", "IC ", 1)
        svalue = svalue.replace("Intergal.", "Intergalactic", 1)
        svalue = svalue.replace("MCG+", "MCG +", 1)
        svalue = svalue.replace("MCG-", "MCG -", 1)
        svalue = svalue.replace("M+", "MCG +", 1)
        svalue = svalue.replace("M-", "MCG -", 1)
        svalue = svalue.replace("MGC ", "MCG ", 1)
        svalue = svalue.replace("Mrk", "MRK", 1)
        svalue = svalue.replace("MRK", "MRK ", 1)
        svalue = svalue.replace("NGC", "NGC ", 1)
        svalue = svalue.replace("PGC", "PGC ", 1)
        svalue = svalue.replace("SDSS", "SDSS ", 1)
        svalue = svalue.replace("UGC", "UGC ", 1)
        if len(svalue) > 4 and svalue.startswith("PGC "):
            svalue = svalue[:4] + svalue[4:].lstrip(" 0")
        if len(svalue) > 4 and svalue.startswith("UGC "):
            svalue = svalue[:4] + svalue[4:].lstrip(" 0")
        if len(svalue) > 5 and svalue.startswith(("MCG +", "MCG -")):
            svalue = svalue[:5] + '-'.join([x.zfill(2) for x in svalue[5:].strip().split("-")])
        if len(svalue) > 5 and svalue.startswith("CGCG "):
            svalue = svalue[:5] + '-'.join([x.zfill(3) for x in svalue[5:].strip().split("-")])
        if (len(svalue) > 1 and svalue.startswith("E")) or (len(svalue) > 3 and svalue.startswith('ESO')):
            if svalue[0] == "E":
                esplit = svalue[1:].split("-")
            else:
                esplit = svalue[3:].split("-")
            if len(esplit) == 2 and is_number(esplit[0].strip()):
                if esplit[1].strip()[0] == 'G':
                    parttwo = esplit[1][1:].strip()
                else:
                    parttwo = esplit[1].strip()
                if is_number(parttwo.strip()):
                    svalue = 'ESO ' + esplit[0].lstrip('0') + '-G' + parttwo.lstrip('0')
        svalue = ' '.join(svalue.split())

        is_abell = svalue.lower().startswith('abell') and is_number(svalue[5:].strip())
        if (not skind and (is_abell or 'cluster' in svalue.lower())):
            skind = 'cluster'

    elif quantity == 'claimedtype':
        isq = False
        svalue = svalue.replace('young', '')
        if '?' in svalue:
            isq = True
            svalue = svalue.strip(' ?')
        for rep in typereps:
            if svalue in typereps[rep]:
                svalue = rep
                break
        if isq:
            svalue = svalue + '?'

    elif quantity in ['ra', 'dec', 'hostra', 'hostdec']:
        if unit == 'floatdegrees':
            deg = float('%g' % Decimal(svalue))
            sig = get_sig_digits(svalue)
            if 'ra' in quantity:
                flhours = deg / 360.0 * 24.0
                hours = floor(flhours)
                minutes = floor((flhours - hours) * 60.0)
                seconds = (flhours * 60.0 - (hours * 60.0 + minutes)) * 60.0
                if seconds > 60.0:
                    raise(ValueError('Invalid seconds value for ' + quantity))
                svalue = str(hours).zfill(2) + ':' + str(minutes).zfill(2) + ':' + zpad(pretty_num(seconds, sig=sig-1))
            elif 'dec' in quantity:
                fldeg = abs(deg)
                degree = floor(fldeg)
                minutes = floor((fldeg - degree) * 60.0)
                seconds = (fldeg * 60.0 - (degree * 60.0 + minutes)) * 60.0
                if seconds > 60.0:
                    raise(ValueError('Invalid seconds value for ' + quantity))
                svalue = (('+' if deg >= 0.0 else '-') + str(degree).strip('+-').zfill(2) + ':' +
                          str(minutes).zfill(2) + ':' + zpad(pretty_num(seconds, sig=sig-1)))

        elif unit == 'nospace' and 'ra' in quantity:
            svalue = svalue[:2] + ':' + svalue[2:4] + ((':' + zpad(svalue[4:])) if len(svalue) > 4 else '')

        elif unit == 'nospace' and 'dec' in quantity:
            if svalue.startswith(('+', '-')):
                svalue = svalue[:3] + ':' + svalue[3:5] + ((':' + zpad(svalue[5:])) if len(svalue) > 5 else '')
            else:
                svalue = '+' + svalue[:2] + ':' + svalue[2:4] + ((':' + zpad(svalue[4:])) if len(svalue) > 4 else '')

        else:
            svalue = svalue.replace(' ', ':')
            if 'dec' in quantity:
                valuesplit = svalue.split(':')
                svalue = (('-' if valuesplit[0].startswith('-') else '+') + valuesplit[0].strip('+-').zfill(2) +
                          (':' + valuesplit[1].zfill(2) if len(valuesplit) > 1 else '') +
                          (':' + zpad(valuesplit[2]) if len(valuesplit) > 2 else ''))

        if 'ra' in quantity:
            sunit = 'hours'
        elif 'dec' in quantity:
            sunit = 'degrees'

        # Correct case of arcseconds = 60.0.
        valuesplit = svalue.split(':')
        if len(valuesplit) == 3 and valuesplit[-1] in ["60.0", "60.", "60"]:
            svalue = valuesplit[0] + ':' + str(Decimal(valuesplit[1]) + Decimal(1.0)) + ':' + "00.0"

        # Strip trailing dots.
        svalue = svalue.rstrip('.')
    elif quantity == 'maxdate' or quantity == 'discoverdate':
        # Make sure month and day have leading zeroes
        sparts = svalue.split('/')
        if len(sparts) >= 2:
            svalue = sparts[0] + '/' + sparts[1].zfill(2)
        if len(sparts) == 3:
            svalue = svalue + '/' + sparts[2].zfill(2)

        if quantity in events[name]:
            for i, ct in enumerate(events[name][quantity]):
                # Only add dates if they have more information
                if len(ct['value'].split('/')) > len(svalue.split('/')):
                    return

    if is_number(svalue):
        svalue = '%g' % Decimal(svalue)
    if serror:
        serror = '%g' % Decimal(serror)

    if quantity in events[name]:
        for i, ct in enumerate(events[name][quantity]):
            if ct['value'] == svalue and sources:
                if 'kind' in ct and skind and ct['kind'] != skind:
                    return
                for source in sources.split(','):
                    if source not in events[name][quantity][i]['source'].split(','):
                        events[name][quantity][i]['source'] += ',' + source
                        if serror and 'error' not in events[name][quantity][i]:
                            events[name][quantity][i]['error'] = serror
                return

    if not sunit:
        sunit = unit

    quantaentry = OrderedDict()
    quantaentry['value'] = svalue
    if serror:
        quantaentry['error'] = serror
    if sources:
        quantaentry['source'] = sources
    if skind:
        quantaentry['kind'] = skind
    if sunit:
        quantaentry['unit'] = sunit
    if lowerlimit:
        quantaentry['lowerlimit'] = lowerlimit
    if upperlimit:
        quantaentry['upperlimit'] = upperlimit
    if extra:
        quantaentry['extra'] = extra
    if (forcereplacebetter or quantity in REPR_BETTER_QUANTITY) and quantity in events[name]:
        newquantities = []
        isworse = True
        if quantity in ['discoverdate', 'maxdate']:
            for ct in events[name][quantity]:
                ctsplit = ct['value'].split('/')
                svsplit = svalue.split('/')
                if len(ctsplit) < len(svsplit):
                    isworse = False
                    continue
                elif len(ctsplit) < len(svsplit) and len(svsplit) == 3:
                    if max(2, get_sig_digits(ctsplit[-1].lstrip('0'))) < max(2, get_sig_digits(svsplit[-1].lstrip('0'))):
                        isworse = False
                        continue
                newquantities.append(ct)
        else:
            newsig = get_sig_digits(svalue)
            for ct in events[name][quantity]:
                if 'error' in ct:
                    if serror:
                        if float(serror) < float(ct['error']):
                            isworse = False
                            continue
                    newquantities.append(ct)
                else:
                    if serror:
                        isworse = False
                        continue
                    oldsig = get_sig_digits(ct['value'])
                    if oldsig >= newsig:
                        newquantities.append(ct)
                    if newsig >= oldsig:
                        isworse = False
        if not isworse:
            newquantities.append(quantaentry)
        events[name][quantity] = newquantities
    else:
        events[name].setdefault(quantity, []).append(quantaentry)


def add_source(events, name, refname='', reference='', url='', bibcode='', secondary='', acknowledgment=''):
    with open(_FILENAME_SOURCE_SYNONYMS, 'r') as f:
        sourcereps = json.loads(f.read(), object_pairs_hook=OrderedDict)

    nsources = len(events[name]['sources']) if 'sources' in events[name] else 0
    if not refname:
        if not bibcode:
            raise(ValueError('Bibcode must be specified if name is not.'))

        if bibcode and len(bibcode) != 19:
            raise(ValueError('Bibcode "' + bibcode + '" must be exactly 19 characters long'))

        refname = bibcode

    if refname.upper().startswith('ATEL') and not bibcode:
        atels_dict = get_atels_dict()
        refname = refname.replace('ATEL', 'ATel').replace('Atel', 'ATel').replace('ATel #', 'ATel ').replace('ATel#', 'ATel').replace('ATel', 'ATel ')
        refname = ' '.join(refname.split())
        atelnum = refname.split()[-1]
        if is_number(atelnum) and atelnum in atels_dict:
            bibcode = atels_dict[atelnum]

    if refname.upper().startswith('CBET') and not bibcode:
        cbets_dict = get_cbets_dict()
        refname = refname.replace('CBET', 'CBET ')
        refname = ' '.join(refname.split())
        cbetnum = refname.split()[-1]
        if is_number(cbetnum) and cbetnum in cbets_dict:
            bibcode = cbets_dict[cbetnum]

    if refname.upper().startswith('IAUC') and not bibcode:
        iaucs_dict = get_iaucs_dict()
        refname = refname.replace('IAUC', 'IAUC ')
        refname = ' '.join(refname.split())
        iaucnum = refname.split()[-1]
        if is_number(iaucnum) and iaucnum in iaucs_dict:
            bibcode = iaucs_dict[iaucnum]

    for rep in sourcereps:
        if refname in sourcereps[rep]:
            refname = rep
            break

    no_source = 'sources' not in events[name]
    no_ref = refname not in [x['name'] for x in events[name]['sources']]
    no_bib = (not bibcode or
              bibcode not in [x['bibcode'] if 'bibcode' in x else '' for x in events[name]['sources']])
    if no_source or (no_ref and no_bib):
        source = str(nsources + 1)
        newsource = OrderedDict()
        newsource['name'] = refname
        if url:
            newsource['url'] = url
        if reference:
            newsource['reference'] = reference
        if bibcode:
            newsource['bibcode'] = bibcode
        if acknowledgment:
            newsource['acknowledgment'] = acknowledgment
        newsource['alias'] = source
        if secondary:
            newsource['secondary'] = True
        events[name].setdefault('sources', []).append(newsource)
    else:
        if refname in [x['name'] for x in events[name]['sources']]:
            source = [x['alias'] for x in events[name]['sources']][
                [x['name'] for x in events[name]['sources']].index(refname)]
        elif bibcode and bibcode in [x['bibcode'] if 'bibcode' in x else '' for x in events[name]['sources']]:
            source = [x['alias'] for x in events[name]['sources']][
                [x['bibcode'] if 'bibcode' in x else '' for x in events[name]['sources']].index(bibcode)]
        else:
            raise(ValueError("Couldn't find source that should exist!"))
    return source


def add_spectrum(events, name, waveunit, fluxunit, wavelengths="", fluxes="", u_time="", time="",
                 instrument="", deredshifted="", dereddened="", errorunit="", errors="", source="",
                 snr="", telescope="", observer="", reducer="", filename="", observatory="",
                 data=""):

    if is_erroneous(events, name, 'spectra', source):
        return

    spectrumentry = OrderedDict()

    if 'spectra' in events[name]:
        for si, spectrum in enumerate(events[name]['spectra']):
            if 'filename' in spectrum and spectrum['filename'] == filename:
                # Copy exclude info
                if 'exclude' in spectrum:
                    spectrumentry['exclude'] = spectrum['exclude']
                # Don't add duplicate spectra
                if 'data' in spectrum:
                    return
                del(events[name]['spectra'][si])
                break

    if not waveunit:
        warnings.warn('No error unit specified, not adding spectrum.')
        return
    if not fluxunit:
        warnings.warn('No flux unit specified, not adding spectrum.')
        return

    if not data or (not wavelengths or not fluxes):
        ValueError('Spectrum must have wavelengths and fluxes set, or data set.')

    if not source:
        ValueError('Spectrum must have source before being added!')

    if deredshifted != '':
        spectrumentry['deredshifted'] = deredshifted
    if dereddened != '':
        spectrumentry['dereddened'] = dereddened
    if instrument:
        spectrumentry['instrument'] = instrument
    if telescope:
        spectrumentry['telescope'] = telescope
    if observatory:
        spectrumentry['observatory'] = observatory
    if u_time:
        spectrumentry['u_time'] = u_time
    if time:
        spectrumentry['time'] = time
    if snr:
        spectrumentry['snr'] = snr
    if observer:
        spectrumentry['observer'] = observer
    if reducer:
        spectrumentry['reducer'] = reducer
    if filename:
        spectrumentry['filename'] = filename

    spectrumentry['waveunit'] = waveunit
    spectrumentry['fluxunit'] = fluxunit
    if data:
        spectrumentry['data'] = data
    else:
        if errors and max([float(x) for x in errors]) > 0.:
            if not errorunit:
                warnings.warn('No error unit specified, not adding spectrum.')
                return
            spectrumentry['errorunit'] = errorunit
            data = [trim_str_arr(wavelengths), trim_str_arr(fluxes), trim_str_arr(errors)]
        else:
            data = [trim_str_arr(wavelengths), trim_str_arr(fluxes)]
        spectrumentry['data'] = [list(i) for i in zip(*data)]
    if source:
        spectrumentry['source'] = source
    events[name].setdefault('spectra', []).append(spectrumentry)


def alias_priority(name, attr):
    if name == attr:
        return 0
    return 1


def archived_task(tasks, args, atask):
    if 'archived' in tasks[atask] and args.archived:
        return True
    if ('archived' in tasks[atask] and tasks[atask]['archived'] and atask not in args.refreshlist.split(',') and not args.fullrefresh):
        return True
    return False


def clean_event(events, dirtyevent):
    bibcodes = []
    name = next(reversed(dirtyevent))

    # This is very hacky and is only necessary because we don't have a proper 'Event' object yet.
    events['temp'] = dirtyevent[name]

    if 'name' not in events['temp']:
        events['temp']['name'] = name
    if 'sources' in events['temp']:
        # Rebuild the sources
        # newsources = []
        oldsources = events['temp']['sources']
        del(events['temp']['sources'])
        for s, source in enumerate(oldsources):
            if 'bibcode' in source:
                bibcodes.append(source['bibcode'])
                add_source(events, 'temp', bibcode=source['bibcode'])
            else:
                add_source(events, 'temp', refname=source['name'], url=source['url'])

    # Clean some legacy fields
    if 'aliases' in events['temp'] and isinstance(events['temp']['aliases'], list):
        source = add_source(events, 'temp', bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
        for alias in events['temp']['aliases']:
            add_quantity(events, 'temp', 'alias', alias, source)
        del(events['temp']['aliases'])

    if (('distinctfrom' in events['temp'] and isinstance(events['temp']['distinctfrom'], list) and
         isinstance(events['temp']['distinctfrom'][0], str))):
            distinctfroms = [x for x in events['temp']['distinctfrom']]
            del(events['temp']['distinctfrom'])
            source = add_source(events, 'temp', bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
            for df in distinctfroms:
                add_quantity(events, 'temp', 'distinctfrom', df, source)

    if (('errors' in events['temp'] and isinstance(events['temp']['errors'], list) and
         'sourcekind' in events['temp']['errors'][0])):
            source = add_source(events, 'temp', bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
            for err in events['temp']['errors']:
                add_quantity(events, 'temp', 'error', err['quantity'], source, kind=err['sourcekind'], extra=err['id'])
            del(events['temp']['errors'])

    if not bibcodes:
        add_source(events, 'temp', bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
        bibcodes = [OSC_BIBCODE]

    for key in list(events['temp'].keys()):
        if key in ['name', 'sources']:
            pass
        elif key == 'photometry':
            for p, photo in enumerate(events['temp']['photometry']):
                if photo['u_time'] == 'JD':
                    events['temp']['photometry'][p]['u_time'] = 'MJD'
                    events['temp']['photometry'][p]['time'] = str(jd_to_mjd(Decimal(photo['time'])))
                if bibcodes and 'source' not in photo:
                    source = add_source(events, 'temp', bibcode=bibcodes[0])
                    events['temp']['photometry'][p]['source'] = source
        else:
            for qi, quantity in enumerate(events['temp'][key]):
                if bibcodes and 'source' not in quantity:
                    source = add_source(events, 'temp', bibcode=bibcodes[0])
                    events['temp'][key][qi]['source'] = source

    cleanevent = events['temp']
    del (events['temp'])
    return OrderedDict([[name, cleanevent]])


def clear_events(events):
    events = OrderedDict((k, OrderedDict([['name', events[k]['name']]] + ([['alias', events[k]['alias']]] if 'alias' in events[k] else []) + [['stub', True]])) for k in events)


def convert_aq_output(row):
    return OrderedDict([(x, str(row[x]) if is_number(row[x]) else row[x]) for x in row.colnames])


def copy_to_event(events, fromname, destname):
    tprint('Copying ' + fromname + ' to event ' + destname)
    newsourcealiases = {}
    keys = list(sorted(events[fromname].keys(), key=lambda key: event_attr_priority(key)))

    if 'sources' in events[fromname]:
        for source in events[fromname]['sources']:
            newsourcealiases[source['alias']] = add_source(
                events, destname, bibcode=source['bibcode'] if 'bibcode' in source else '',
                refname=source['name'] if 'name' in source else '',
                reference=source['reference'] if 'reference' in source else '',
                url=source['url'] if 'url' in source else '')

    for key in keys:
        if key not in ['name', 'sources']:
            for item in events[fromname][key]:
                # isd = False
                sources = []
                if 'source' not in item:
                    ValueError("Item has no source!")
                for sid in item['source'].split(','):
                    if sid == 'D':
                        sources.append('D')
                    elif sid in newsourcealiases:
                        sources.append(newsourcealiases[sid])
                    else:
                        ValueError("Couldn't find source alias!")
                sources = uniq_cdl(sources)

                if key == 'photometry':
                    add_photometry(
                        events, destname, u_time=null_field(item, "u_time"), time=null_field(item, "time"),
                        e_time=null_field(item, "e_time"), telescope=null_field(item, "telescope"),
                        instrument=null_field(item, "instrument"), band=null_field(item, "band"),
                        magnitude=null_field(item, "magnitude"), e_magnitude=null_field(item, "e_magnitude"),
                        source=sources, upperlimit=null_field(item, "upperlimit"), system=null_field(item, "system"),
                        observatory=null_field(item, "observatory"), observer=null_field(item, "observer"),
                        host=null_field(item, "host"), survey=null_field(item, "survey"))
                elif key == 'spectra':
                    add_spectrum(
                        events, destname, null_field(item, "waveunit"), null_field(item, "fluxunit"), data=null_field(item, "data"),
                        u_time=null_field(item, "u_time"), time=null_field(item, "time"),
                        instrument=null_field(item, "instrument"), deredshifted=null_field(item, "deredshifted"),
                        dereddened=null_field(item, "dereddened"), errorunit=null_field(item, "errorunit"),
                        source=sources, snr=null_field(item, "snr"),
                        telescope=null_field(item, "telescope"), observer=null_field(item, "observer"),
                        reducer=null_field(item, "reducer"), filename=null_field(item, "filename"),
                        observatory=null_field(item, "observatory"))
                elif key == 'errors':
                    add_quantity(
                        events, destname, key, item['value'], sources,
                        kind=null_field(item, "kind"), extra=null_field(item, "extra"))
                else:
                    add_quantity(
                        events, destname, key, item['value'], sources, error=null_field(item, "error"),
                        unit=null_field(item, "unit"), kind=null_field(item, "kind"))


def ct_priority(name, attr):
    aliases = attr['source'].split(',')
    max_source_year = -10000
    vaguetypes = ['CC', 'I']
    if attr['value'] in vaguetypes:
        return -max_source_year
    for alias in aliases:
        if alias == 'D':
            continue
        source = get_source_by_alias(name, alias)
        if 'bibcode' in source:
            source_year = get_source_year(source)
            if source_year > max_source_year:
                max_source_year = source_year
    return -max_source_year


def derive_and_sanitize(tasks, args, events):
    biberrordict = {
        "2012Sci..337..942D": "2012Sci...337..942D",
        "2012MNRAS.420.1135": "2012MNRAS.420.1135S",
        "2014MNRAS.438,368": "2014MNRAS.438..368T",
        "2006ApJ...636...400Q": "2006ApJ...636..400Q",
        "0609268": "2007AJ....133...58K",
        "2004MNRAS.tmp..131P": "2004MNRAS.352..457P",
        "2013MNRAS.tmp.1499F": "2013MNRAS.433.1312F",
        "1991MNRAS.247P.410B": "1991A&A...247..410B",
        "2011Sci.333..856S": "2011Sci...333..856S"
    }

    # Calculate some columns based on imported data, sanitize some fields
    for name in events:
        aliases = get_aliases(events, name, includename=False)
        if name not in aliases:
            if 'sources' in events[name]:
                add_quantity(events, name, 'alias', name, '1')
            else:
                source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
                add_quantity(events, name, 'alias', name, source)

        if ((name.startswith('SN') and is_number(name[2:6]) and 'discoverdate' in events[name] and
             int(events[name]['discoverdate'][0]['value'].split('/')[0]) >= 2016 and
             not any(['AT' in x for x in aliases]))):
            source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
            add_quantity(events, name, 'alias', 'AT' + name[2:], source)

        events[name]['alias'] = list(sorted(events[name]['alias'], key=lambda key: alias_priority(name, key)))
        aliases = get_aliases(events, name)

        set_first_max_light(events, name)

        if 'claimedtype' in events[name]:
            events[name]['claimedtype'] = list(sorted(events[name]['claimedtype'], key=lambda key: ct_priority(name, key)))
        if 'discoverdate' not in events[name]:
            prefixes = ['MLS', 'SSS', 'CSS']
            for alias in aliases:
                for prefix in prefixes:
                    if alias.startswith(prefix) and is_number(alias.replace(prefix, '')[:2]):
                        discoverdate = '/'.join(['20' + alias.replace(prefix, '')[:2],
                                                alias.replace(prefix, '')[2:4],
                                                alias.replace(prefix, '')[4:6]])
                        if args.verbose:
                            tprint('Added discoverdate from name: ' + discoverdate)
                        source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
                        add_quantity(events, name, 'discoverdate', discoverdate, source)
                        break
                if 'discoverdate' in events[name]:
                    break
        if 'discoverdate' not in events[name]:
            prefixes = ['ASASSN-', 'PS1-', 'PS1', 'PS', 'iPTF', 'PTF', 'SCP-']
            for alias in aliases:
                for prefix in prefixes:
                    if alias.startswith(prefix) and is_number(alias.replace(prefix, '')[:2]):
                        discoverdate = '20' + alias.replace(prefix, '')[:2]
                        if args.verbose:
                            tprint('Added discoverdate from name: ' + discoverdate)
                        source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
                        add_quantity(events, name, 'discoverdate', discoverdate, source)
                        break
                if 'discoverdate' in events[name]:
                    break
        if 'discoverdate' not in events[name]:
            prefixes = ['SNF']
            for alias in aliases:
                for prefix in prefixes:
                    if alias.startswith(prefix) and is_number(alias.replace(prefix, '')[:4]):
                        discoverdate = '/'.join([alias.replace(prefix, '')[:4],
                                                 alias.replace(prefix, '')[4:6],
                                                 alias.replace(prefix, '')[6:8]])
                        if args.verbose:
                            tprint('Added discoverdate from name: ' + discoverdate)
                        source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
                        add_quantity(events, name, 'discoverdate', discoverdate, source)
                        break
                if 'discoverdate' in events[name]:
                    break
        if 'discoverdate' not in events[name]:
            prefixes = ['AT', 'SN']
            for alias in aliases:
                for prefix in prefixes:
                    if alias.startswith(prefix) and is_number(alias.replace(prefix, '')[:4]):
                        discoverdate = alias.replace(prefix, '')[:4]
                        if args.verbose:
                            tprint('Added discoverdate from name: ' + discoverdate)
                        source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
                        add_quantity(events, name, 'discoverdate', discoverdate, source)
                        break
                if 'discoverdate' in events[name]:
                    break
        if 'ra' not in events[name] or 'dec' not in events[name]:
            prefixes = ['PSN J', 'MASJ', 'CSS', 'SSS', 'MASTER OT J']
            for alias in aliases:
                for prefix in prefixes:
                    if alias.startswith(prefix) and is_number(alias.replace(prefix, '')[:6]):
                        noprefix = alias.split(':')[-1].replace(prefix, '').replace('.', '')
                        decsign = '+' if '+' in noprefix else '-'
                        noprefix = noprefix.replace('+', '|').replace('-', '|')
                        nops = noprefix.split('|')
                        if len(nops) < 2:
                            continue
                        rastr = nops[0]
                        decstr = nops[1]
                        ra = ':'.join([rastr[:2], rastr[2:4], rastr[4:6]]) + ('.' + rastr[6:] if len(rastr) > 6 else '')
                        dec = decsign + ':'.join([decstr[:2], decstr[2:4], decstr[4:6]]) + ('.' + decstr[6:] if len(decstr) > 6 else '')
                        if args.verbose:
                            tprint('Added ra/dec from name: ' + ra + ' ' + dec)
                        source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
                        add_quantity(events, name, 'ra', ra, source)
                        add_quantity(events, name, 'dec', dec, source)
                        break
                if 'ra' in events[name]:
                    break

        no_host = ('host' not in events[name] or
                   not any([x['value'] == 'Milky Way' for x in events[name]['host']]))
        if ('ra' in events[name] and 'dec' in events[name] and no_host):
            from astroquery.irsa_dust import IrsaDust
            extinctions_dict = get_extinctions_dict()
            if name not in extinctions_dict:
                try:
                    ra_dec = events[name]['ra'][0]['value'] + " " + events[name]['dec'][0]['value']
                    result = IrsaDust.get_query_table(ra_dec, section='ebv')
                except:
                    warnings.warn("Coordinate lookup for " + name + " failed in IRSA.")
                else:
                    ebv = result['ext SandF mean'][0]
                    ebverr = result['ext SandF std'][0]
                    extinctions_dict[name] = [ebv, ebverr]
            if name in extinctions_dict:
                source = add_source(events, name, bibcode='2011ApJ...737..103S')
                add_quantity(events, name, 'ebv', str(extinctions_dict[name][0]), source, error=str(extinctions_dict[name][1]))
        if 'claimedtype' in events[name]:
            events[name]['claimedtype'][:] = [ct for ct in events[name]['claimedtype'] if (ct['value'] != '?' and ct['value'] != '-')]
        if 'claimedtype' not in events[name] and name.startswith('AT'):
            source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
            add_quantity(events, name, 'claimedtype', 'Candidate', source)
        if 'redshift' not in events[name] and 'velocity' in events[name]:
            # Find the "best" velocity to use for this
            bestsig = 0
            for hv in events[name]['velocity']:
                sig = get_sig_digits(hv['value'])
                if sig > bestsig:
                    besthv = hv['value']
                    bestsig = sig
            if bestsig > 0 and is_number(besthv):
                voc = float(besthv)*1.e5/CLIGHT
                source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
                add_quantity(events, name, 'redshift', pretty_num(sqrt((1. + voc)/(1. - voc)) - 1., sig=bestsig), source, kind='heliocentric')
        if 'redshift' not in events[name] and has_task(tasks, args, 'nedd') and 'host' in events[name]:
            from astropy.cosmology import Planck15 as cosmo, z_at_value
            import statistics
            reference = "NED-D"
            refurl = "http://ned.ipac.caltech.edu/Library/Distances/"
            for host in events[name]['host']:
                if host['value'] in nedddict:
                    secondarysource = add_source(events, name, refname=reference, url=refurl, secondary=True)
                    meddist = statistics.median(nedddict[host['value']])
                    redshift = pretty_num(z_at_value(cosmo.comoving_distance, float(meddist) * un.Mpc), sig=get_sig_digits(str(meddist)))
                    add_quantity(events, name, 'redshift', redshift, secondarysource, kind='host')
        if 'maxabsmag' not in events[name] and 'maxappmag' in events[name] and 'lumdist' in events[name]:
            # Find the "best" distance to use for this
            bestsig = 0
            for ld in events[name]['lumdist']:
                sig = get_sig_digits(ld['value'])
                if sig > bestsig:
                    bestld = ld['value']
                    bestsig = sig
            if bestsig > 0 and is_number(bestld) and float(bestld) > 0.:
                source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
                pnum = float(events[name]['maxappmag'][0]['value']) - 5.0*(log10(float(bestld)*1.0e6) - 1.0)
                pnum = pretty_num(pnum, sig=bestsig)
                add_quantity(events, name, 'maxabsmag', pnum, source)
        if 'redshift' in events[name]:
            # Find the "best" redshift to use for this
            (bestz, bestkind, bestsig) = get_best_redshift(events, name)
            if bestsig > 0:
                bestz = float(bestz)
                if 'velocity' not in events[name]:
                    source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
                    pnum = CLIGHT/KM*((bestz + 1.)**2. - 1.)/((bestz + 1.)**2. + 1.)
                    pnum = pretty_num(pnum, sig=bestsig)
                    add_quantity(events, name, 'velocity', pnum, source, kind=PREF_KINDS[bestkind])
                if bestz > 0.:
                    from astropy.cosmology import Planck15 as cosmo
                    if 'lumdist' not in events[name]:
                        dl = cosmo.luminosity_distance(bestz)
                        source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
                        add_quantity(events, name, 'lumdist', pretty_num(dl.value, sig=bestsig), source, kind=PREF_KINDS[bestkind])
                        if 'maxabsmag' not in events[name] and 'maxappmag' in events[name]:
                            source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
                            pnum = pretty_num(float(events[name]['maxappmag'][0]['value']) - 5.0*(log10(dl.to('pc').value) - 1.0), sig=bestsig)
                            add_quantity(events, name, 'maxabsmag', pnum, source)
                    if 'comovingdist' not in events[name]:
                        dl = cosmo.comoving_distance(bestz)
                        source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
                        add_quantity(events, name, 'comovingdist', pretty_num(dl.value, sig=bestsig), source)
        if 'photometry' in events[name]:
            events[name]['photometry'].sort(key=lambda x: ((float(x['time']) if isinstance(x['time'], str) else
                min([float(y) for y in x['time']])) if 'time' in x else 0.0,
                x['band'] if 'band' in x else '', float(x['magnitude']) if 'magnitude' in x else ''))
        if 'spectra' in events[name] and list(filter(None, ['time' in x for x in events[name]['spectra']])):
            events[name]['spectra'].sort(key=lambda x: (float(x['time']) if 'time' in x else 0.0))
        if 'sources' in events[name]:
            bibauthor_dict = get_bibauthor_dict()
            for source in events[name]['sources']:
                if 'bibcode' in source:
                    import urllib
                    from html import unescape
                    # First sanitize the bibcode
                    if len(source['bibcode']) != 19:
                        source['bibcode'] = urllib.parse.unquote(unescape(source['bibcode'])).replace('A.A.', 'A&A')
                    if source['bibcode'] in biberrordict:
                        source['bibcode'] = biberrordict[source['bibcode']]

                    if source['bibcode'] not in bibauthor_dict:
                        bibcode = source['bibcode']
                        adsquery = ('http://adsabs.harvard.edu/cgi-bin/nph-abs_connect?db_key=ALL&version=1&bibcode=' +
                                    urllib.parse.quote(bibcode) + '&data_type=Custom&format=%253m%20%25(y)')
                        response = urllib.request.urlopen(adsquery)
                        html = response.read().decode('utf-8')
                        hsplit = html.split("\n")
                        if len(hsplit) > 5:
                            bibcodeauthor = hsplit[5]
                        else:
                            bibcodeauthor = ''

                        if not bibcodeauthor:
                            warnings.warn("Bibcode didn't return authors, not converting this bibcode.")

                        bibauthor_dict[bibcode] = unescape(bibcodeauthor).strip()

            for source in events[name]['sources']:
                if 'bibcode' in source and source['bibcode'] in bibauthor_dict and bibauthor_dict[source['bibcode']]:
                    source['reference'] = bibauthor_dict[source['bibcode']]
                    if 'name' not in source and source['bibcode']:
                        source['name'] = source['bibcode']
        if 'redshift' in events[name]:
            events[name]['redshift'] = list(sorted(events[name]['redshift'], key=lambda key: frame_priority(key)))
        if 'velocity' in events[name]:
            events[name]['velocity'] = list(sorted(events[name]['velocity'], key=lambda key: frame_priority(key)))
        if 'claimedtype' in events[name]:
            events[name]['claimedtype'] = list(sorted(events[name]['claimedtype'], key=lambda key: ct_priority(name, key)))

        events[name] = OrderedDict(sorted(events[name].items(), key=lambda key: event_attr_priority(key[0])))


def delete_old_event_files():
    # Delete all old event JSON files
    files = repo_file_list()
    for f in files:
        os.remove(f)


def do_task(tasks, args, checktask, task, quiet=False):
    global currenttask
    dotask = has_task(tasks, args, task) and checktask == task
    if dotask and not quiet:
        currenttask = (tasks[task]['nicename'] if tasks[task]['nicename'] else task).replace('%pre', 'Updating' if args.update else 'Loading')
    return dotask


def event_attr_priority(attr):
    if attr == 'photometry':
        return 'zzzzzzzy'
    if attr == 'spectra':
        return 'zzzzzzzz'
    if attr == 'name':
        return 'aaaaaaaa'
    if attr == 'sources':
        return 'aaaaaaab'
    if attr == 'alias':
        return 'aaaaaaac'
    return attr


def event_exists(events, name):
    if name in events:
        return True
    for ev in events:
        if name in get_aliases(events, ev):
            return True
    return False


def frame_priority(attr):
    if 'kind' in attr:
        if attr['kind'] in PREF_KINDS:
            return PREF_KINDS.index(attr['kind'])
        else:
            return len(PREF_KINDS)
    return len(PREF_KINDS)


def get_aliases(events, name, includename=True):
    if 'alias' in events[name]:
        aliases = [x['value'] for x in events[name]['alias']]
        if includename and name not in aliases:
            return [name] + aliases
        return aliases
    if includename:
        return [name]
    return []


def get_atels_dict():
    # path = '../atels.json'
    if os.path.isfile(_FILENAME_ATELS):
        with open(_FILENAME_ATELS, 'r') as f:
            atels_dict = json.loads(f.read(), object_pairs_hook=OrderedDict)
    else:
        atels_dict = OrderedDict()
    return atels_dict


def get_bibauthor_dict():
    # path = '../bibauthors.json'
    if os.path.isfile(_FILENAME_BIBAUTHORS):
        with open(_FILENAME_BIBAUTHORS, 'r') as f:
            bibauthor_dict = json.loads(f.read(), object_pairs_hook=OrderedDict)
    else:
        bibauthor_dict = OrderedDict()
    return bibauthor_dict


def get_cbets_dict():
    # path = '../cbets.json'
    if os.path.isfile(_FILENAME_CBETS):
        with open(_FILENAME_CBETS, 'r') as f:
            cbets_dict = json.loads(f.read(), object_pairs_hook=OrderedDict)
    else:
        cbets_dict = OrderedDict()
    return cbets_dict


def get_extinctions_dict():
    # path = '../extinctions.json'
    if os.path.isfile(_FILENAME_EXTINCT):
        with open(_FILENAME_EXTINCT, 'r') as f:
            extinctions_dict = json.loads(f.read(), object_pairs_hook=OrderedDict)
    else:
        extinctions_dict = OrderedDict()
    return extinctions_dict


def get_iaucs_dict():
    # path = '../iaucs.json'
    if os.path.isfile(_FILENAME_IAUCS):
        with open(_FILENAME_IAUCS, 'r') as f:
            iaucs_dict = json.loads(f.read(), object_pairs_hook=OrderedDict)
    else:
        iaucs_dict = OrderedDict()
    return iaucs_dict


def get_best_redshift(events, name):
    bestsig = -1
    bestkind = 10
    for z in events[name]['redshift']:
        kind = PREF_KINDS.index(z['kind'] if 'kind' in z else '')
        sig = get_sig_digits(z['value'])
        if sig > bestsig and kind <= bestkind:
            bestz = z['value']
            bestkind = kind
            bestsig = sig

    return (bestz, bestkind, bestsig)


def get_first_light(events, name):
    if 'photometry' not in events[name]:
        return (None, None)

    eventphoto = [(Decimal(x['time']) if isinstance(x['time'], str) else Decimal(min(float(y) for y in x['time'])),
                  x['source']) for x in events[name]['photometry'] if 'upperlimit' not in x and
                  'time' in x and 'u_time' in x and x['u_time'] == 'MJD']
    if not eventphoto:
        return (None, None)
    flmjd = min([x[0] for x in eventphoto])
    flindex = [x[0] for x in eventphoto].index(flmjd)
    flmjd = float(flmjd)
    flsource = eventphoto[flindex][1]
    return (astrotime(flmjd, format='mjd').datetime, flsource)


def get_max_light(name):
    if 'photometry' not in events[name]:
        return (None, None, None, None)

    eventphoto = [(x['u_time'], x['time'], Decimal(x['magnitude']), x['band'] if 'band' in x else '', x['source']) for x in events[name]['photometry'] if
                  ('magnitude' in x and 'time' in x and 'u_time' in x and 'upperlimit' not in x)]
    if not eventphoto:
        return (None, None, None, None)

    mlmag = None
    for mb in MAX_BANDS:
        leventphoto = [x for x in eventphoto if x[3] in mb]
        if leventphoto:
            mlmag = min([x[2] for x in leventphoto])
            eventphoto = leventphoto
            break

    if not mlmag:
        mlmag = min([x[2] for x in eventphoto])

    mlindex = [x[2] for x in eventphoto].index(mlmag)
    mlband = eventphoto[mlindex][3]
    mlsource = eventphoto[mlindex][4]

    if eventphoto[mlindex][0] == 'MJD':
        mlmjd = float(eventphoto[mlindex][1])
        return (astrotime(mlmjd, format='mjd').datetime, mlmag, mlband, mlsource)
    else:
        return (None, mlmag, mlband, mlsource)


def get_preferred_name(name):
    if name not in events:
        matches = []
        for event in events:
            aliases = get_aliases(events, event)
            if len(aliases) > 1 and name in aliases:
                return event
        return name
    else:
        return name


def get_source_by_alias(name, alias):
    for source in events[name]['sources']:
        if source['alias'] == alias:
            return source
    raise(ValueError('Source alias not found!'))


def get_source_year(source):
    if 'bibcode' in source:
        if is_number(source['bibcode'][:4]):
            return int(source['bibcode'][:4])
        else:
            return -10000
    raise(ValueError('No bibcode available for source!'))


def has_task(tasks, args, task):
    return task in tasks and (not args.update or tasks[task]['update'])


def is_erroneous(events, name, field, sources):
    if 'errors' in events[name]:
        for alias in sources.split(','):
            source = get_source_by_alias(name, alias)
            if (('bibcode' in source and source['bibcode'] in
                [x['value'] for x in events[name]['errors'] if x['kind'] == 'bibcode' and x['extra'] == field])):
                    return True
            if (('name' in source and source['name'] in
                [x['value'] for x in events[name]['errors'] if x['kind'] == 'name' and x['extra'] == field])):
                    return True
    return False


def jd_to_mjd(jd):
    return jd - Decimal(2400000.5)


def journal_events(tasks, args, events, clear=True):
    if 'writeevents' in tasks:
        write_all_events(events, args)
    if clear:
        clear_events(events)


def load_event_from_file(events, args, tasks, name='', location='', clean=False, delete=True, append=False):
    if not name and not location:
        raise ValueError('Either event name or location must be specified to load event')

    path = ''
    namepath = ''
    repo_folders = get_repo_folders()
    if location:
        path = location
    if name:
        indir = '../'
        for rep in repo_folders:
            filename = get_event_filename(name)
            newpath = indir + rep + '/' + filename + '.json'
            if os.path.isfile(newpath):
                namepath = newpath

    if not path and not namepath:
        return False
    else:
        newevent = ''
        newevent2 = ''
        if path or namepath:
            if name in events:
                del events[name]

        if path and namepath:
            with open(path, 'r') as f, open(namepath, 'r') as nf:
                newevent = json.loads(f.read(), object_pairs_hook=OrderedDict)
                newevent2 = json.loads(nf.read(), object_pairs_hook=OrderedDict)
        elif path:
            with open(path, 'r') as f:
                newevent = json.loads(f.read(), object_pairs_hook=OrderedDict)
        elif namepath:
            with open(namepath, 'r') as f:
                newevent = json.loads(f.read(), object_pairs_hook=OrderedDict)

        if newevent:
            if clean:
                newevent = clean_event(events, newevent)
            name = next(reversed(newevent))
            if append:
                indir = '../'
                for rep in repo_folders:
                    filename = get_event_filename(name)
                    newpath = indir + rep + '/' + filename + '.json'
                    if os.path.isfile(newpath):
                        namepath = newpath
                if namepath:
                    with open(namepath, 'r') as f:
                        newevent2 = json.loads(f.read(), object_pairs_hook=OrderedDict)
                        namename = next(reversed(newevent2))

            if newevent2:
                # Needs to be fixed
                newevent = OrderedDict([['temp', newevent[name]]])
                copy_to_event(events, 'temp', namename)
            else:
                events.update(newevent)

            if args.verbose and not args.travis:
                tprint('Loaded ' + name)

        if 'writeevents' in tasks and delete and namepath:
            os.remove(namepath)
        return name


def load_stubs(tasks, args, events):
    global currenttask
    currenttask = 'Loading event stubs'
    files = repo_file_list()

    # try:
    #    namepath = '../names.min.json'
    #    with open(namepath, 'r') as f:
    #        names = json.loads(f.read(), object_pairs_hook=OrderedDict)
    #    for fi in tq(files):
    #        name = os.path.basename(os.path.splitext(fi)[0])
    #        if name not in names:
    #            name = name.replace("_", "/")
    #        events[name] = OrderedDict(([['name', name], ['alias', [OrderedDict(([['value', x]])) for x in names[name]]], ['stub', True]]))
    # except:
    #    events = OrderedDict()
    for fi in tq(files, currenttask):
        fname = fi
        if '.gz' in fi:
            import shutil
            import gzip
            fname = fi.replace('.gz', '')
            with gzip.open(fi, 'rb') as f_in, open(fname, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(fi)
        name = os.path.basename(os.path.splitext(fname)[0]).replace('.json', '')
        name = add_event(tasks, args, events, name, delete=False, loadifempty=False)
        events[name] = OrderedDict(([['name', events[name]['name']]] + ([['alias', events[name]['alias']]] if 'alias' in events[name] else []) + [['stub', True]]))


def load_cached_url(args, url, filepath, timeout=120, write=True):
    import codecs
    from hashlib import md5
    filemd5 = ''
    filetxt = ''
    if not args.refresh and os.path.isfile(filepath):
        with codecs.open(filepath, 'r', encoding='utf8') as f:
            filetxt = f.read()
            if args.update:
                filemd5 = md5(filetxt.encode('utf-8')).hexdigest()

    try:
        import requests
        session = requests.Session()
        response = session.get(url, timeout=timeout)
        if any([x.status_code == 307 for x in response.history]):
            raise
        txt = response.text
        newmd5 = md5(txt.encode('utf-8')).hexdigest()
        # tprint(filemd5 + ": " + newmd5)
        if args.update and newmd5 == filemd5:
            tprint('Skipping file in "' + currenttask + '," local and remote copies identical [' + newmd5 + '].')
            return False
    except:
        return filetxt
    else:
        if write:
            with codecs.open(filepath, 'w', encoding='utf8') as f:
                f.write(txt)
    return txt


def make_date_string(year, month='', day=''):
    if not year:
        raise ValueError('At least the year must be specified when constructing date string')
    datestring = str(year)
    if month:
        datestring = datestring + '/' + str(month).zfill(2)
    if day:
        datestring = datestring + '/' + str(day).zfill(2)

    return datestring


# Merge and remove duplicate events
def merge_duplicates(tasks, args, events):
    if not len(events):
        load_stubs(tasks, args, events)
    currenttask = 'Merging duplicate events'
    keys = list(sorted(list(events.keys())))
    for n1, name1 in enumerate(tq(keys[:], currenttask)):
        if name1 not in events:
            continue
        allnames1 = get_aliases(events, name1) + (['AT' + name1[2:]] if (name1.startswith('SN') and is_number(name1[2:6])) else [])
        for name2 in keys[n1+1:]:
            if name2 not in events or name1 == name2:
                continue
            allnames2 = get_aliases(events, name2) + (['AT' + name2[2:]] if (name2.startswith('SN') and is_number(name2[2:6])) else [])
            if any(i in allnames1 for i in allnames2):
                tprint('Found single event with multiple entries (' + name1 + ' and ' + name2 + '), merging.')
                load1 = load_event_from_file(events, args, tasks, name1, delete=True)
                load2 = load_event_from_file(events, args, tasks, name2, delete=True)
                if load1 and load2:
                    priority1 = 0
                    priority2 = 0
                    for an in allnames1:
                        if len(an) >= 2 and an.startswith(('SN', 'AT')):
                            priority1 = priority1 + 1
                    for an in allnames2:
                        if len(an) >= 2 and an.startswith(('SN', 'AT')):
                            priority2 = priority2 + 1

                    if priority1 > priority2:
                        copy_to_event(events, name2, name1)
                        keys.append(name1)
                        del(events[name2])
                    else:
                        copy_to_event(events, name1, name2)
                        keys.append(name2)
                        del(events[name1])
                else:
                    print('Duplicate already deleted')
                journal_events(tasks, args, events)


def name_clean(name):
    newname = name.strip(' ;,*')
    if newname.startswith('MASJ'):
        newname = newname.replace('MASJ', 'MASTER OT J', 1)
    if newname.startswith('MASTER') and is_number(newname[7]):
        newname = newname.replace('MASTER', 'MASTER OT J', 1)
    if newname.startswith('MASTER OT J '):
        newname = newname.replace('MASTER OT J ', 'MASTER OT J', 1)
    if newname.startswith('Psn'):
        newname = newname.replace('Psn', 'PSN', 1)
    if newname.startswith('PSNJ'):
        newname = newname.replace('PSNJ', 'PSN J', 1)
    if newname.startswith('TCPJ'):
        newname = newname.replace('TCPJ', 'TCP J', 1)
    if newname.startswith('SMTJ'):
        newname = newname.replace('SMTJ', 'SMT J', 1)
    if newname.startswith('PSN20J'):
        newname = newname.replace('PSN20J', 'PSN J', 1)
    if newname.startswith('ASASSN') and newname[6] != '-':
        newname = newname.replace('ASASSN', 'ASASSN-', 1)
    if newname.startswith('ROTSE3J'):
        newname = newname.replace('ROTSE3J', 'ROTSE3 J', 1)
    if newname.startswith('SNHunt'):
        newname = newname.replace('SNHunt', 'SNhunt', 1)
    if newname.startswith('ptf'):
        newname = newname.replace('ptf', 'PTF', 1)
    if newname.startswith('PTF '):
        newname = newname.replace('PTF ', 'PTF', 1)
    if newname.startswith('iPTF '):
        newname = newname.replace('iPTF ', 'iPTF', 1)
    if newname.startswith('SNHunt'):
        newname = newname.replace('SNHunt', 'SNhunt', 1)
    if newname.startswith('PESSTOESO'):
        newname = newname.replace('PESSTOESO', 'PESSTO ESO ', 1)
    if newname.startswith('snf'):
        newname = newname.replace('snf', 'SNF', 1)
    if newname.startswith('SNF') and is_number(newname[3:]) and len(newname) >= 12:
        newname = 'SNF' + newname[3:11] + '-' + newname[11:]
    if newname.startswith(('MASTER OT J', 'ROTSE3 J')):
        prefix = newname.split('J')[0]
        coords = newname.split('J')[-1].strip()
        decsign = '+' if '+' in coords else '-'
        coordsplit = coords.replace('+', '-').split('-')
        if '.' not in coordsplit[0] and len(coordsplit[0]) > 6 and '.' not in coordsplit[1] and len(coordsplit[1]) > 6:
            newname = (prefix + 'J' + coordsplit[0][:6] + '.' + coordsplit[0][6:] + decsign + coordsplit[1][:6] + '.' + coordsplit[1][6:])
    if newname.startswith('Gaia ') and is_number(newname[3:4]) and len(newname) > 5:
        newname = newname.replace('Gaia ', 'Gaia', 1)
    if len(newname) <= 4 and is_number(newname):
        newname = 'SN' + newname + 'A'
    if len(newname) > 4 and is_number(newname[:4]) and not is_number(newname[4:]):
        newname = 'SN' + newname
    if newname.startswith('sn') and is_number(newname[2:6]) and len(newname) > 6:
        newname = newname.replace('sn', 'SN', 1)
    if newname.startswith('SN ') and is_number(newname[3:7]) and len(newname) > 7:
        newname = newname.replace('SN ', 'SN', 1)
    if newname.startswith('SN') and is_number(newname[2:6]) and len(newname) == 7 and newname[6].islower():
        newname = 'SN' + newname[2:6] + newname[6].upper()
    elif (newname.startswith('SN') and is_number(newname[2:6]) and (len(newname) == 8 or len(newname) == 9) and newname[6:].isupper()):
        newname = 'SN' + newname[2:6] + newname[6:].lower()

    newname = (' '.join(newname.split())).strip()
    return newname


def null_field(obj, field):
    return obj[field] if field in obj else ''


def same_tag_num(photo, val, tag, canbelist=False):
    issame = (
        (tag not in photo and not val) or
        (tag in photo and not val) or
        (tag in photo and
         ((not canbelist and Decimal(photo[tag]) == Decimal(val)) or
          (canbelist and
           ((isinstance(photo[tag], str) and isinstance(val, str) and Decimal(photo[tag]) == Decimal(val)) or
            (isinstance(photo[tag], list) and isinstance(val, list) and photo[tag] == val))))))
    return issame


def same_tag_str(photo, val, tag):
    issame = ((tag not in photo and not val) or (tag in photo and not val) or (tag in photo and photo[tag] == val))
    return issame


def set_first_max_light(events, name):
    if 'maxappmag' not in events[name]:
        (mldt, mlmag, mlband, mlsource) = get_max_light(name)
        if mldt:
            source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
            add_quantity(events, name, 'maxdate', make_date_string(mldt.year, mldt.month, mldt.day), uniq_cdl([source, mlsource]))
        if mlmag:
            source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
            add_quantity(events, name, 'maxappmag', pretty_num(mlmag), uniq_cdl([source, mlsource]))
        if mlband:
            source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
            add_quantity(events, name, 'maxband', mlband, uniq_cdl([source, mlsource]))

    if 'discoverdate' not in events[name] or max([len(x['value'].split('/')) for x in events[name]['discoverdate']]) < 3:
        (fldt, flsource) = get_first_light(events, name)
        if fldt:
            source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
            add_quantity(events, name, 'discoverdate', make_date_string(fldt.year, fldt.month, fldt.day), uniq_cdl([source, flsource]))

    if 'discoverdate' not in events[name] and 'spectra' in events[name]:
        minspecmjd = float("+inf")
        for spectrum in events[name]['spectra']:
            if 'time' in spectrum and 'u_time' in spectrum:
                if spectrum['u_time'] == 'MJD':
                    mjd = float(spectrum['time'])
                elif spectrum['u_time'] == 'JD':
                    mjd = float(jd_to_mjd(Decimal(spectrum['time'])))
                else:
                    continue

                if mjd < minspecmjd:
                    minspecmjd = mjd
                    minspecsource = spectrum['source']

        if minspecmjd < float("+inf"):
            fldt = astrotime(minspecmjd, format='mjd').datetime
            source = add_source(events, name, bibcode=OSC_BIBCODE, refname=OSC_NAME, url=OSC_URL, secondary=True)
            add_quantity(events, name, 'discoverdate', make_date_string(fldt.year, fldt.month, fldt.day), 'D,' + minspecsource)


def set_preferred_names(tasks, args, events):
    if not len(events):
        load_stubs(tasks, args, events)
    for name in list(sorted(list(events.keys()))):
        if name not in events:
            continue
        newname = ''
        aliases = get_aliases(events, name)
        if len(aliases) <= 1:
            continue
        if (name.startswith('SN') and ((is_number(name[2:6]) and not is_number(name[6:])) or
                                       (is_number(name[2:5]) and not is_number(name[5:])))):
            continue
        for alias in aliases:
            if (alias[:2] == 'SN' and ((is_number(alias[2:6]) and not is_number(alias[6:])) or
                                       (is_number(alias[2:5]) and not is_number(alias[5:])))):
                newname = alias
                break
        if not newname and 'discoverer' in events[name]:
            discoverer = ','.join([x['value'].upper() for x in events[name]['discoverer']])
            if 'ASAS' in discoverer:
                for alias in aliases:
                    if 'ASASSN' in alias.upper():
                        newname = alias
                        break
            if not newname and 'OGLE' in discoverer:
                for alias in aliases:
                    if 'OGLE' in alias.upper():
                        newname = alias
                        break
            if not newname and 'CRTS' in discoverer:
                for alias in aliases:
                    if True in [x in alias.upper() for x in ['CSS', 'MLS', 'SSS', 'SNHUNT']]:
                        newname = alias
                        break
            if not newname and 'PS1' in discoverer:
                for alias in aliases:
                    if 'PS1' in alias.upper():
                        newname = alias
                        break
            if not newname and 'PTF' in discoverer:
                for alias in aliases:
                    if 'PTF' in alias.upper():
                        newname = alias
                        break
            if not newname and 'GAIA' in discoverer:
                for alias in aliases:
                    if 'GAIA' in alias.upper():
                        newname = alias
                        break
        if not newname:
            for alias in aliases:
                # Always prefer another alias over PSN
                if name.startswith('PSN'):
                    newname = alias
                    break
        if newname and name != newname:
            # Make sure new name doesn't already exist
            if load_event_from_file(events, args, tasks, newname):
                continue
            if load_event_from_file(events, args, tasks, name, delete=True):
                tprint('Changing event name (' + name + ') to preferred name (' + newname + ').')
                events[newname] = events[name]
                events[newname]['name'] = newname
                del(events[name])
                journal_events(tasks, args, events)


def snname(string):
    newstring = string.replace(' ', '').upper()
    if (newstring[:2] == "SN"):
        head = newstring[:6]
        tail = newstring[6:]
        if len(tail) >= 2 and tail[1] != '?':
            tail = tail.lower()
        newstring = head + tail

    return newstring


def trim_str_arr(arr, length=10):
    return [str(round_sig(float(x), length)) if
            (len(x) > length and len(str(round_sig(float(x), length))) < len(x))
            else x for x in arr]


def uniq_cdl(values):
    return ','.join(list(OrderedDict.fromkeys(values).keys()))


def utf8(x):
    return str(x, 'utf-8')


def write_all_events(events, args, empty=False, gz=False, bury=False):
    import codecs
    repo_folders = get_repo_folders()
    non_sne_types = None
    if bury:
        with open(_FILENAME_NON_SNE_TYPES, 'r') as f:
            non_sne_types = json.loads(f.read(), object_pairs_hook=OrderedDict)
            non_sne_types = [x.upper() for x in non_sne_types]

    # Write it all out!
    for name in events:
        if 'stub' in events[name]:
            if not empty:
                continue
            else:
                del(events[name]['stub'])
        if args.verbose and not args.travis:
            tprint('Writing ' + name)
        filename = get_event_filename(name)

        outdir = '../'
        if 'discoverdate' in events[name]:
            repo_years = get_repo_years()
            for r, year in enumerate(repo_years):
                if int(events[name]['discoverdate'][0]['value'].split('/')[0]) <= year:
                    outdir += repo_folders[r]
                    break
        else:
            outdir += str(repo_folders[0])

        # Delete non-SN events here without IAU designations (those with only banned types)
        if bury:
            buryevent = False
            nonsneprefixes = ('PNVJ', 'PNV J', 'OGLE-2013-NOVA')
            if name.startswith(nonsneprefixes):
                tprint('Burying ' + name + ', non-SNe prefix.')
                continue
            if 'claimedtype' in events[name] and not (name.startswith('SN') and is_number(name[2:6])):
                for ct in events[name]['claimedtype']:
                    if ct['value'].upper() not in non_sne_types and ct['value'].upper() != 'CANDIDATE':
                        buryevent = False
                        break
                    if ct['value'].upper() in non_sne_types:
                        buryevent = True
                if buryevent:
                    tprint('Burying ' + name + ' (' + ct['value'] + ').')
                    outdir = '../sne-boneyard'

        jsonstring = json.dumps({name: events[name]}, indent='\t', separators=(',', ':'), ensure_ascii=False)

        path = outdir + '/' + filename + '.json'
        with codecs.open(path, 'w', encoding='utf8') as f:
            f.write(jsonstring)

        if gz:
            if os.path.getsize(path) > 90000000:
                import shutil
                import gzip
                if not args.travis:
                    tprint('Compressing ' + name)
                with open(path, 'rb') as f_in, gzip.open(path + '.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                os.remove(path)
                os.system('cd ' + outdir + '; git rm ' + filename + '.json; git add -f ' + filename + '.json.gz; cd ' + '../scripts')
