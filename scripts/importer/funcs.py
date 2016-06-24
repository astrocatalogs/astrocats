"""Utility functions for OSC import.
"""

import os
import json
import warnings
# import sys
from cdecimal import Decimal
from collections import OrderedDict
from math import log10, sqrt, floor
from astropy.time import Time as astrotime
from astropy import units

from scripts import FILENAME
from . constants import OSC_BIBCODE, OSC_NAME, OSC_URL, CLIGHT, PREF_KINDS, \
    KM, MAX_BANDS
from .. utils import bandrepf, bandmetaf, is_number, \
    get_sig_digits, pretty_num, round_sig, tprint, zpad


def add_photometry(events, name, time = "", u_time = "MJD", e_time = "", telescope = "", instrument = "", band = "",
                   magnitude = "", e_magnitude = "", source = "", upperlimit = False, system = "", scorrected = "",
                   observatory = "", observer = "", host = False, includeshost = False, survey = "", kcorrected = "",
                   flux = "", fluxdensity = "", e_flux = "", e_fluxdensity = "", u_flux = "", u_fluxdensity = "", frequency = "",
                   u_frequency = "", counts = "", e_counts = "", nhmw = "", photonindex = "", unabsorbedflux = "",
                   e_unabsorbedflux = "", energy = "", u_energy = "", e_lower_magnitude = "", e_upper_magnitude = "",
                   e_lower_time = "", e_upper_time = "", mcorrected = ""):
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

    if events[name].is_erroneous('photometry', source):
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
                 same_tag_num(photo, energy, 'energy', canbelist = True) and
                 same_tag_num(photo, frequency, 'frequency') and
                 same_tag_num(photo, photonindex, 'photonindex') and
                 same_tag_num(photo, e_magnitude, 'e_magnitude') and
                 same_tag_num(photo, e_lower_time, 'e_lower_time') and
                 same_tag_num(photo, e_upper_time, 'e_upper_time') and
                 same_tag_num(photo, e_lower_magnitude, 'e_lower_magnitude') and
                 same_tag_num(photo, e_upper_magnitude, 'e_upper_magnitude') and
                 same_tag_num(photo, e_flux, 'e_flux') and
                 same_tag_num(photo, e_unabsorbedflux, 'e_unabsorbedflux') and
                 same_tag_num(photo, e_fluxdensity, 'e_fluxdensity') and
                 same_tag_num(photo, e_counts, 'e_counts') and
                 same_tag_str(photo, u_flux, 'u_flux') and
                 same_tag_str(photo, u_fluxdensity, 'u_fluxdensity') and
                 same_tag_str(photo, u_frequency, 'u_frequency') and
                 same_tag_str(photo, u_energy, 'u_energy') and
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
    if e_lower_time:
        photoentry['e_lower_time'] = str(e_lower_time)
    if e_upper_time:
        photoentry['e_upper_time'] = str(e_upper_time)
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
    if kcorrected:
        photoentry['kcorrected'] = kcorrected
    if scorrected:
        photoentry['scorrected'] = scorrected
    if mcorrected:
        photoentry['mcorrected'] = mcorrected
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


def add_spectrum(events, name, waveunit, fluxunit, wavelengths="", fluxes="", u_time="", time="",
                 instrument="", deredshifted="", dereddened="", errorunit="", errors="", source="",
                 snr="", telescope="", observer="", survey="", reducer="", filename="", observatory="",
                 data=""):

    if events[name].is_erroneous('spectra', source):
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
                del events[name]['spectra'][si]
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
    if survey:
        spectrumentry['survey'] = survey
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


'''
def archived_task(tasks, args, atask):
    if 'archived' in tasks[atask] and args.archived:
        return True
    if (('archived' in tasks[atask] and tasks[atask]['archived'] and
         atask not in args.refresh_list.split(',') and not args.full_refresh)):
        return True
    return False
'''


def convert_aq_output(row):
    return OrderedDict([(x, str(row[x]) if is_number(row[x]) else row[x]) for x in row.colnames])


def copy_to_event(events, fromname, destname):
    tprint('Copying ' + fromname + ' to event ' + destname)
    newsourcealiases = {}
    keys = list(sorted(events[fromname].keys(), key=lambda xx: event_attr_priority(xx)))

    if 'sources' in events[fromname]:
        for source in events[fromname]['sources']:
            newsourcealiases[source['alias']] = events[destname].add_source(
                bibcode=source['bibcode'] if 'bibcode' in source else '',
                refname=source['name'] if 'name' in source else '',
                reference=source['reference'] if 'reference' in source else '',
                url=source['url'] if 'url' in source else '')

    if 'errors' in events[fromname]:
        for err in events[fromname]['errors']:
            events[destname].setdefault('errors',[]).append(err)

    for key in keys:
        if key not in ['schema', 'name', 'sources', 'errors']:
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
                    events[destname].add_quantity(
                        key, item['value'], sources,
                        kind=null_field(item, "kind"), extra=null_field(item, "extra"))
                else:
                    events[destname].add_quantity(
                        key, item['value'], sources, error=null_field(item, "error"),
                        unit = null_field(item, "unit"), probability=null_field(item, "probability"), kind=null_field(item, "kind"))


def ct_priority(events, name, attr):
    aliases = attr['source'].split(',')
    max_source_year = -10000
    vaguetypes = ['CC', 'I']
    if attr['value'] in vaguetypes:
        return -max_source_year
    for alias in aliases:
        if alias == 'D':
            continue
        source = events[name].get_source_by_alias(alias)
        if 'bibcode' in source:
            source_year = get_source_year(source)
            if source_year > max_source_year:
                max_source_year = source_year
    return -max_source_year


def derive_and_sanitize(tasks, args, events, extinctions_dict, bibauthor_dict, nedd_dict):
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
        aliases = events[name].get_aliases(includename=False)
        if name not in aliases:
            if 'sources' in events[name]:
                events[name].add_quantity('alias', name, '1')
            else:
                source = events[name].add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
                events[name].add_quantity('alias', name, source)

        if ((name.startswith('SN') and is_number(name[2:6]) and 'discoverdate' in events[name] and
             int(events[name]['discoverdate'][0]['value'].split('/')[0]) >= 2016 and
             not any(['AT' in x for x in aliases]))):
            source = events[name].add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
            events[name].add_quantity('alias', 'AT' + name[2:], source)

        events[name]['alias'] = list(sorted(events[name]['alias'], key=lambda key: alias_priority(name, key)))
        aliases = events[name].get_aliases()

        set_first_max_light(events, name)

        if 'claimedtype' in events[name]:
            events[name]['claimedtype'] = list(sorted(events[name]['claimedtype'], key=lambda key: ct_priority(events, name, key)))
        if 'discoverdate' not in events[name]:
            prefixes = ['MLS', 'SSS', 'CSS', 'GRB ']
            for alias in aliases:
                for prefix in prefixes:
                    if alias.startswith(prefix) and is_number(alias.replace(prefix, '')[:2]):
                        discoverdate = '/'.join(['20' + alias.replace(prefix, '')[:2],
                                                alias.replace(prefix, '')[2:4],
                                                alias.replace(prefix, '')[4:6]])
                        if args.verbose:
                            tprint ('Added discoverdate from name [' + alias + ']: ' + discoverdate)
                        source = events[name].add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
                        events[name].add_quantity('discoverdate', discoverdate, source, derived = True)
                        break
                if 'discoverdate' in events[name]:
                    break
        if 'discoverdate' not in events[name]:
            prefixes = ['ASASSN-', 'PS1-', 'PS1', 'PS', 'iPTF', 'PTF', 'SCP-', 'SNLS-', 'SPIRITS', 'LSQ', 'DES', 'SNHiTS',
                'GND', 'GNW', 'GSD', 'GSW', 'EGS', 'COS']
            for alias in aliases:
                for prefix in prefixes:
                    if alias.startswith(prefix) and is_number(alias.replace(prefix, '')[:2]):
                        discoverdate = '20' + alias.replace(prefix, '')[:2]
                        if args.verbose:
                            tprint ('Added discoverdate from name [' + alias + ']: ' + discoverdate)
                        source = events[name].add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
                        events[name].add_quantity('discoverdate', discoverdate, source, derived = True)
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
                            tprint ('Added discoverdate from name [' + alias + ']: ' + discoverdate)
                        source = events[name].add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
                        events[name].add_quantity('discoverdate', discoverdate, source, derived = True)
                        break
                if 'discoverdate' in events[name]:
                    break
        if 'discoverdate' not in events[name]:
            prefixes = ['PTFS', 'SNSDF']
            for alias in aliases:
                for prefix in prefixes:
                    if alias.startswith(prefix) and is_number(alias.replace(prefix, '')[:2]):
                        discoverdate = '/'.join(['20' + alias.replace(prefix, '')[:2],
                            alias.replace(prefix, '')[2:4]])
                        if args.verbose:
                            tprint ('Added discoverdate from name [' + alias + ']: ' + discoverdate)
                        source = events[name].add_source(bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
                        events[name].add_quantity('discoverdate', discoverdate, source, derived = True)
                        break
                if 'discoverdate' in events[name]:
                    break
        if 'discoverdate' not in events[name]:
            prefixes = ['AT', 'SN', 'OGLE-', 'SM ', 'KSN-']
            for alias in aliases:
                for prefix in prefixes:
                    if (alias.startswith(prefix) and is_number(alias.replace(prefix, '')[:4])
			and '.' not in alias.replace(prefix, '')[:4]):
                        discoverdate = alias.replace(prefix, '')[:4]
                        if args.verbose:
                            tprint ('Added discoverdate from name [' + alias + ']: ' + discoverdate)
                        source = events[name].add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
                        events[name].add_quantity('discoverdate', discoverdate, source, derived = True)
                        break
                if 'discoverdate' in events[name]:
                    break
        if 'ra' not in events[name] or 'dec' not in events[name]:
            prefixes = ['PSN J', 'MASJ', 'CSS', 'SSS', 'MASTER OT J', 'HST J', 'TCP J', 'MACS J', '2MASS J', 'EQ J', 'CRTS J', 'SMT J']
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
                        source = events[name].add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
                        events[name].add_quantity('ra', ra, source, derived = True)
                        events[name].add_quantity('dec', dec, source, derived = True)
                        break
                if 'ra' in events[name]:
                    break

        no_host = ('host' not in events[name] or
                   not any([x['value'] == 'Milky Way' for x in events[name]['host']]))
        if ('ra' in events[name] and 'dec' in events[name] and no_host):
            from astroquery.irsa_dust import IrsaDust
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
                sources = uniq_cdl([events[name].add_source(bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True),
                    events[name].add_source(bibcode = '2011ApJ...737..103S')])
                events[name].add_quantity('ebv', str(extinctionsdict[name][0]), sources, error = str(extinctionsdict[name][1]), derived = True)
        if 'host' in events[name] and ('hostra' not in events[name] or 'hostdec' not in events[name]):
            for host in events[name]['host']:
                alias = host['value']
                if ' J' in alias and is_number(alias.split(' J')[-1][:6]):
                    noprefix = alias.split(' J')[-1].split(':')[-1].replace('.', '')
                    decsign = '+' if '+' in noprefix else '-'
                    noprefix = noprefix.replace('+','|').replace('-','|')
                    nops = noprefix.split('|')
                    if len(nops) < 2:
                        continue
                    rastr = nops[0]
                    decstr = nops[1]
                    hostra = ':'.join([rastr[:2], rastr[2:4], rastr[4:6]]) + ('.' + rastr[6:] if len(rastr) > 6 else '')
                    hostdec = decsign + ':'.join([decstr[:2], decstr[2:4], decstr[4:6]]) + ('.' + decstr[6:] if len(decstr) > 6 else '')
                    if args.verbose:
                        tprint ('Added hostra/hostdec from name: ' + hostra + ' ' + hostdec)
                    source = events[name].add_source(bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
                    events[name].add_quantity('hostra', hostra, source, derived = True)
                    events[name].add_quantity('hostdec', hostdec, source, derived = True)
                    break
                if 'hostra' in events[name]:
                    break
        if 'claimedtype' in events[name]:
            events[name]['claimedtype'][:] = [ct for ct in events[name]['claimedtype'] if (ct['value'] != '?' and ct['value'] != '-')]
            if not len(events[name]['claimedtype']):
                del(events[name]['claimedtype'])
        if 'claimedtype' not in events[name] and name.startswith('AT'):
            source = events[name].add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
            events[name].add_quantity('claimedtype', 'Candidate', source)
        if 'redshift' not in events[name] and 'velocity' in events[name]:
            # Find the "best" velocity to use for this
            bestsig = 0
            for hv in events[name]['velocity']:
                sig = get_sig_digits(hv['value'])
                if sig > bestsig:
                    besthv = hv['value']
                    bestsrc = hv['source']
                    bestsig = sig
            if bestsig > 0 and is_number(besthv):
                voc = float(besthv)*1.e5/CLIGHT
                source = events[name].add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
                sources = uniq_cdl([source] + bestsrc.split(','))
                events[name].add_quantity('redshift', pretty_num(sqrt((1. + voc)/(1. - voc)) - 1., sig = bestsig), sources,
                    kind = 'heliocentric', derived = True)
        if 'redshift' not in events[name] and has_task(tasks, args, 'nedd') and 'host' in events[name]:
            from astropy.cosmology import Planck15 as cosmo, z_at_value
            import statistics
            reference = "NED-D"
            refurl = "http://ned.ipac.caltech.edu/Library/Distances/"
            for host in events[name]['host']:
                if host['value'] in nedd_dict:
                    source = events[name].add_source(bibcode = '2015arXiv150201589P')
                    secondarysource = events[name].add_source(srcname=reference, url=refurl, secondary=True)
                    meddist = statistics.median(nedd_dict[host['value']])
                    redshift = pretty_num(z_at_value(cosmo.comoving_distance, float(meddist) * units.Mpc), sig=get_sig_digits(str(meddist)))
                    events[name].add_quantity(name, 'redshift', redshift, uniq_cdl([source,secondarysource]), kind = 'host', derived = True)
        if 'maxabsmag' not in events[name] and 'maxappmag' in events[name] and 'lumdist' in events[name]:
            # Find the "best" distance to use for this
            bestsig = 0
            for ld in events[name]['lumdist']:
                sig = get_sig_digits(ld['value'])
                if sig > bestsig:
                    bestld = ld['value']
                    bestsrc = ld['source']
                    bestsig = sig
            if bestsig > 0 and is_number(bestld) and float(bestld) > 0.:
                source = events[name].add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
                sources = uniq_cdl([source] + bestsrc.split(','))
                pnum = float(events[name]['maxappmag'][0]['value']) - 5.0*(log10(float(bestld)*1.0e6) - 1.0)
                pnum = pretty_num(pnum, sig=bestsig)
                events[name].add_quantity('maxabsmag', pnum, sources, derived = True)
        if 'redshift' in events[name]:
            # Find the "best" redshift to use for this
            (bestz, bestkind, bestsig) = get_best_redshift(events, name)
            if bestsig > 0:
                bestz = float(bestz)
                if 'velocity' not in events[name]:
                    source = events[name].add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
                    pnum = CLIGHT/KM*((bestz + 1.)**2. - 1.)/((bestz + 1.)**2. + 1.)
                    pnum = pretty_num(pnum, sig=bestsig)
                    events[name].add_quantity('velocity', pnum, source, kind=PREF_KINDS[bestkind])
                if bestz > 0.:
                    from astropy.cosmology import Planck15 as cosmo
                    if 'lumdist' not in events[name]:
                        dl = cosmo.luminosity_distance(bestz)
                        sources = [events[name].add_source(bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True),
                            events[name].add_source(bibcode = '2015arXiv150201589P')]
                        sources = uniq_cdl(sources + bestsrc.split(','))
                        events[name].add_quantity('lumdist', pretty_num(dl.value, sig = bestsig), sources,
                            kind = PREF_KINDS[bestkind], derived = True)
                        if 'maxabsmag' not in events[name] and 'maxappmag' in events[name]:
                            source = events[name].add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
                            pnum = pretty_num(float(events[name]['maxappmag'][0]['value']) - 5.0*(log10(dl.to('pc').value) - 1.0), sig=bestsig)
                            events[name].add_quantity('maxabsmag', pnum, sources, derived = True)
                    if 'comovingdist' not in events[name]:
                        cd = cosmo.comoving_distance(bestz)
                        sources = [events[name].add_source(bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True),
                            events[name].add_source(bibcode = '2015arXiv150201589P')]
                        sources = uniq_cdl(sources + bestsrc.split(','))
                        events[name].add_quantity('comovingdist', pretty_num(cd.value, sig = bestsig), sources, derived = True)
        if all([x in events[name] for x in ['ra', 'dec', 'hostra', 'hostdec']]):
            # For now just using first coordinates that appear in entry
            try:
                c1 = coord(ra=events[name]['ra'][0]['value'], dec=events[name]['dec'][0]['value'], unit=(un.hourangle, un.deg))
                c2 = coord(ra=events[name]['hostra'][0]['value'], dec=events[name]['hostdec'][0]['value'], unit=(un.hourangle, un.deg))
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                pass
            else:
                sources = uniq_cdl([events[name].add_source(bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)] +
                    events[name]['ra'][0]['source'].split(',') + events[name]['dec'][0]['source'].split(',') +
                    events[name]['hostra'][0]['source'].split(',') + events[name]['hostdec'][0]['source'].split(','))
                if 'hostoffsetang' not in events[name]:
                    events[name].add_quantity('hostoffsetang', pretty_num(Decimal(hypot(c1.ra.degree - c2.ra.degree,
                        c1.dec.degree - c2.dec.degree))*Decimal(3600.)), sources, derived = True, unit = 'arcseconds')
                if 'comovingdist' in events[name] and 'redshift' in events[name] and 'hostoffsetdist' not in events[name]:
                    offsetsig = get_sig_digits(events[name]['hostoffsetang'][0]['value'])
                    sources = uniq_cdl(sources.split(',') +
                        events[name]['comovingdist'][0]['source'].split(',') + events[name]['redshift'][0]['source'].split(','))
                    events[name].add_quantity('hostoffsetdist',
                        pretty_num(float(events[name]['hostoffsetang'][0]['value']) / 3600. * (pi / 180.) *
                        float(events[name]['comovingdist'][0]['value']) * 1000. / (1.0 + float(events[name]['redshift'][0]['value'])),
                        sig = offsetsig), sources)
        if 'photometry' in events[name]:
            events[name]['photometry'].sort(
                key=lambda x: ((float(x['time']) if isinstance(x['time'], str) else
                                min([float(y) for y in x['time']])) if 'time' in x else 0.0,
                               x['band'] if 'band' in x else '', float(x['magnitude']) if 'magnitude' in x else ''))
        if 'spectra' in events[name] and list(filter(None, ['time' in x for x in events[name]['spectra']])):
            events[name]['spectra'].sort(key=lambda x: (float(x['time']) if 'time' in x else 0.0))
        if 'sources' in events[name]:
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
            events[name]['claimedtype'] = list(sorted(events[name]['claimedtype'], key=lambda key: ct_priority(events, name, key)))

        events[name] = OrderedDict(sorted(events[name].items(), key=lambda key: event_attr_priority(key[0])))

    return events, extinctions_dict, bibauthor_dict

'''
def do_task(tasks, args, checktask, task, quiet=False):
    """
    """
    global currenttask
    dotask = has_task(tasks, args, task) and checktask == task
    if dotask and not quiet:
        currenttask = (tasks[task]['nicename'] if tasks[task]['nicename'] else task).replace('%pre', 'Updating' if args.update else 'Loading')
    return dotask
'''


def event_attr_priority(attr):
    if attr == 'photometry':
        return 'zzy'
    if attr == 'spectra':
        return 'zzz'
    if attr == 'schema':
        return 'aaa'
    if attr == 'name':
        return 'aab'
    if attr == 'sources':
        return 'aac'
    if attr == 'alias':
        return 'aad'
    return attr


def event_exists(events, name):
    if name in events:
        return True
    for ev in events:
        if name in events[ev].get_aliases():
            return True
    return False


def frame_priority(attr):
    if 'kind' in attr:
        if attr['kind'] in PREF_KINDS:
            return PREF_KINDS.index(attr['kind'])
        else:
            return len(PREF_KINDS)
    return len(PREF_KINDS)


def get_atels_dict():
    # path = '../atels.json'
    if os.path.isfile(FILENAME.ATELS):
        with open(FILENAME.ATELS, 'r') as f:
            atels_dict = json.loads(f.read(), object_pairs_hook=OrderedDict)
    else:
        atels_dict = OrderedDict()
    return atels_dict


def get_bibauthor_dict():
    # path = '../bibauthors.json'
    if os.path.isfile(FILENAME.BIBAUTHORS):
        with open(FILENAME.BIBAUTHORS, 'r') as f:
            bibauthor_dict = json.loads(f.read(), object_pairs_hook=OrderedDict)
    else:
        bibauthor_dict = OrderedDict()
    return bibauthor_dict


def get_cbets_dict():
    # path = '../cbets.json'
    if os.path.isfile(FILENAME.CBETS):
        with open(FILENAME.CBETS, 'r') as f:
            cbets_dict = json.loads(f.read(), object_pairs_hook=OrderedDict)
    else:
        cbets_dict = OrderedDict()
    return cbets_dict


def get_extinctions_dict():
    # path = '../extinctions.json'
    if os.path.isfile(FILENAME.EXTINCT):
        with open(FILENAME.EXTINCT, 'r') as f:
            extinctions_dict = json.loads(f.read(), object_pairs_hook=OrderedDict)
    else:
        extinctions_dict = OrderedDict()
    return extinctions_dict


def get_iaucs_dict():
    # path = '../iaucs.json'
    if os.path.isfile(FILENAME.IAUCS):
        with open(FILENAME.IAUCS, 'r') as f:
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
            bestsrc = z['source']

    return (bestz, bestkind, bestsig, bestsrc)


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


def get_max_light(events, name):
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


def get_preferred_name(events, name):
    if name not in events:
        # matches = []
        for event in events:
            aliases = events[event].get_aliases()
            if len(aliases) > 1 and name in aliases:
                return event
        return name
    else:
        return name


def get_source_year(source):
    if 'bibcode' in source:
        if is_number(source['bibcode'][:4]):
            return int(source['bibcode'][:4])
        else:
            return -10000
    raise ValueError('No bibcode available for source!')


def has_task(tasks, args, task):
    return task in tasks and (not args.update or tasks[task]['update'])


def jd_to_mjd(jd):
    return jd - Decimal(2400000.5)


def load_cached_url(args, current_task, url, filepath, timeout=120, write=True, failhard = False):
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
        response.raise_for_status()
        for x in response.history:
            x.raise_for_status()
            if x.status_code == 500 or x.status_code == 307 or x.status_code == 404:
                raise
        txt = response.text
        newmd5 = md5(txt.encode('utf-8')).hexdigest()
        # tprint(filemd5 + ": " + newmd5)
        if args.update and newmd5 == filemd5:
            tprint('Skipping file in "' + current_task + '," local and remote copies identical [' + newmd5 + '].')
            return False
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        if failhard:
            return ''
        return filetxt
    else:
        if write:
            with codecs.open(filepath, 'w', encoding='utf8') as f:
                f.write(txt if txt else filetxt)
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


def name_clean(name):
    newname = name.strip(' ;,*')
    if newname.startswith('NAME '):
        newname = newname.replace('NAME ', '', 1)
    if newname.endswith(' SN'):
        newname = newname.replace(' SN', '')
    if newname.endswith(':SN'):
        newname = newname.replace(':SN', '')
    if newname.startswith('MASJ'):
        newname = newname.replace('MASJ', 'MASTER OT J', 1)
    if newname.startswith('MASTER') and is_number(newname[7]):
        newname = newname.replace('MASTER', 'MASTER OT J', 1)
    if newname.startswith('MASTER OT J '):
        newname = newname.replace('MASTER OT J ', 'MASTER OT J', 1)
    if newname.startswith('OGLE '):
        newname = newname.replace('OGLE ', 'OGLE-', 1)
    if newname.startswith('OGLE-') and len(newname) != 16:
        namesp = newname.split('-')
        if len(namesp[1]) == 4 and is_number(namesp[1]) and is_number(namesp[3]):
            newname = 'OGLE-' + namesp[1] + '-SN-' + namesp[3].zfill(3)
    if newname.startswith('SN SDSS'):
        newname = newname.replace('SN SDSS ', 'SDSS', 1)
    if newname.startswith('SDSS '):
        newname = newname.replace('SDSS ', 'SDSS', 1)
    if newname.startswith('SDSS'):
        namesp = newname.split('-')
        if len(namesp) == 3 and is_number(namesp[0][4:]) and is_number(namesp[1]) and is_number(namesp[2]):
            newname = namesp[0] + '-' + namesp[1] + '-' + namesp[2].zfill(3)
    if newname.startswith('SDSS-II SN'):
        namesp = newname.split()
        if len(namesp) == 3 and is_number(namesp[2]):
            newname = 'SDSS-II SN ' + namesp[2].lstrip('0')
    if newname.startswith('SN CL'):
        newname = newname.replace('SN CL', 'CL', 1)
    if newname.startswith('SN HiTS '):
        newname = newname.replace('SN HiTS ', 'SNHiTS', 1)
    if newname.startswith('GAIA'):
        newname = newname.replace('GAIA', 'Gaia', 1)
    if newname.startswith('Gaia '):
        newname = newname.replace('Gaia ', 'Gaia', 1)
    if newname.startswith('Gaia'):
        newname = 'Gaia' + newname[4:].lower()
    if newname.startswith('GRB'):
        newname = newname.replace('GRB', 'GRB ', 1)
    if newname.startswith('GRB ') and is_number(newname[4:].strip()):
        newname = 'GRB ' + newname[4:].strip() + 'A'
    if newname.startswith('LSQ '):
        newname = newname.replace('LSQ ', 'LSQ', 1)
    if newname.startswith('KSN '):
        newname = newname.replace('KSN ', 'KSN-', 1)
    if newname.startswith('SNSDF '):
        newname = newname.replace(' ', '')
    if newname.startswith('SNSDF'):
        namesp = newname.split('.')
        if len(namesp[0]) == 9:
            newname = namesp[0] + '-' + namesp[1].zfill(2)
    if newname.startswith('HFF '):
        newname = newname.replace(' ', '')
    if newname.startswith('SN HST'):
        newname = newname.replace('SN HST', 'HST', 1)
    if newname.startswith('HST ') and newname[4] != 'J':
        newname = newname.replace('HST ', 'HST J', 1)
    if newname.startswith('SNLS') and newname[4] != '-':
        newname = newname.replace('SNLS', 'SNLS-', 1)
    if newname.startswith('SNLS- '):
        newname = newname.replace('SNLS- ', 'SNLS-', 1)
    if newname.startswith('CRTS CSS'):
        newname = newname.replace('CRTS CSS', 'CSS', 1)
    if newname.startswith('CRTS MLS'):
        newname = newname.replace('CRTS MLS', 'MLS', 1)
    if newname.startswith('CRTS SSS'):
        newname = newname.replace('CRTS SSS', 'SSS', 1)
    if newname.startswith(('CSS', 'MLS', 'SSS')):
        newname = newname.replace(' ', ':').replace('J', '')
    if newname.startswith('SN HFF'):
        newname = newname.replace('SN HFF', 'HFF', 1)
    if newname.startswith('SN GND'):
        newname = newname.replace('SN GND', 'GND', 1)
    if newname.startswith('SN SCP'):
        newname = newname.replace('SN SCP', 'SCP', 1)
    if newname.startswith('SN UDS'):
        newname = newname.replace('SN UDS', 'UDS', 1)
    if newname.startswith('SCP') and newname[3] != '-':
        newname = newname.replace('SCP', 'SCP-', 1)
    if newname.startswith('SCP- '):
        newname = newname.replace('SCP- ', 'SCP-', 1)
    if newname.startswith('PS 1'):
        newname = newname.replace('PS 1', 'PS1', 1)
    if newname.startswith('PS1 SN PS'):
        newname = newname.replace('PS1 SN PS', 'PS', 1)
    if newname.startswith('PS1 SN'):
        newname = newname.replace('PS1 SN', 'PS1', 1)
    if newname.startswith('PSN K'):
        newname = newname.replace('PSN K', 'K', 1)
    if newname.startswith('K') and is_number(newname[1:5]):
        namesp = newname.split('-')
        if len(namesp[0]) == 5:
            newname = namesp[0] + '-' + namesp[1].zfill(3)
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
    if newname.startswith('SN ASASSN'):
        newname = newname.replace('SN ASASSN', 'ASASSN', 1)
    if newname.startswith('ASASSN '):
        newname = newname.replace('ASASSN ', 'ASASSN-', 1).replace('--', '-')
    if newname.startswith('ASASSN') and newname[6] != '-':
        newname = newname.replace('ASASSN', 'ASASSN-', 1)
    if newname.startswith('ROTSE3J'):
        newname = newname.replace('ROTSE3J', 'ROTSE3 J', 1)
    if newname.startswith('MACSJ'):
        newname = newname.replace('MACSJ', 'MACS J', 1)
    if newname.startswith('MWSNR'):
        newname = newname.replace('MWSNR', 'MWSNR ', 1)
    if newname.startswith('SN HUNT'):
        newname = newname.replace('SN HUNT', 'SNhunt', 1)
    if newname.startswith('SN Hunt'):
        newname = newname.replace(' ', '')
    if newname.startswith('SNHunt'):
        newname = newname.replace('SNHunt', 'SNhunt', 1)
    if newname.startswith('SNhunt '):
        newname = newname.replace('SNhunt ', 'SNhunt', 1)
    if newname.startswith('ptf'):
        newname = newname.replace('ptf', 'PTF', 1)
    if newname.startswith('SN PTF'):
        newname = newname.replace('SN PTF', 'PTF', 1)
    if newname.startswith('PTF '):
        newname = newname.replace('PTF ', 'PTF', 1)
    if newname.startswith('iPTF '):
        newname = newname.replace('iPTF ', 'iPTF', 1)
    if newname.startswith('PESSTOESO'):
        newname = newname.replace('PESSTOESO', 'PESSTO ESO ', 1)
    if newname.startswith('snf'):
        newname = newname.replace('snf', 'SNF', 1)
    if newname.startswith('SNF '):
        newname = newname.replace('SNF ', 'SNF', 1)
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
    if newname.startswith('Sn ') and is_number(newname[3:7]) and len(newname) > 7:
        newname = newname.replace('Sn ', 'SN', 1)
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


def radec_clean(svalue, quantity, unit = ''):
    if unit == 'floatdegrees':
        if not is_number(svalue):
            return (svalue, unit)
        deg = float('%g' % Decimal(svalue))
        sig = get_sig_digits(svalue)
        if 'ra' in quantity:
            flhours = deg / 360.0 * 24.0
            hours = floor(flhours)
            minutes = floor((flhours - hours) * 60.0)
            seconds = (flhours * 60.0 - (hours * 60.0 + minutes)) * 60.0
            hours = 0 if hours < 1.e-6 else hours
            minutes = 0 if minutes < 1.e-6 else minutes
            seconds = 0.0 if seconds < 1.e-6 else seconds
            if seconds > 60.0:
                raise(ValueError('Invalid seconds value for ' + quantity))
            svalue = str(hours).zfill(2) + ':' + str(minutes).zfill(2) + ':' + zpad(pretty_num(seconds, sig = sig - 1))
        elif 'dec' in quantity:
            fldeg = abs(deg)
            degree = floor(fldeg)
            minutes = floor((fldeg - degree) * 60.0)
            seconds = (fldeg * 60.0 - (degree * 60.0 + minutes)) * 60.0
            if seconds > 60.0:
                raise(ValueError('Invalid seconds value for ' + quantity))
            svalue = (('+' if deg >= 0.0 else '-') + str(degree).strip('+-').zfill(2) + ':' +
                str(minutes).zfill(2) + ':' + zpad(pretty_num(seconds, sig = sig - 1)))
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

    return (svalue, sunit)


def host_clean(name):
    newname = name.strip(' ;,*')

    # Handle some special cases
    hostcases = {'M051a':'M51A', 'M051b':'M51B'}
    for k in hostcases:
        if newname == k:
            newname = hostcases[k]

    # Some general cases
    newname = newname.strip("()").replace('  ', ' ', 1)
    newname = newname.replace("ABELL", "Abell", 1)
    newname = newname.replace("Abell", "Abell ", 1)
    newname = newname.replace("APMUKS(BJ)", "APMUKS(BJ) ", 1)
    newname = newname.replace("ARP", "ARP ", 1)
    newname = newname.replace("CGCG", "CGCG ", 1)
    newname = newname.replace("HOLM", "HOLM ", 1)
    newname = newname.replace("IC", "IC ", 1)
    newname = newname.replace("Intergal.", "Intergalactic", 1)
    newname = newname.replace("MCG+", "MCG +", 1)
    newname = newname.replace("MCG-", "MCG -", 1)
    newname = newname.replace("M+", "MCG +", 1)
    newname = newname.replace("M-", "MCG -", 1)
    newname = newname.replace("MGC ", "MCG ", 1)
    newname = newname.replace("Mrk", "MRK", 1)
    newname = newname.replace("MRK", "MRK ", 1)
    newname = newname.replace("NGC", "NGC ", 1)
    newname = newname.replace("PGC", "PGC ", 1)
    newname = newname.replace("SDSS", "SDSS ", 1)
    newname = newname.replace("UGC", "UGC ", 1)
    if newname.startswith('MESSIER '):
        newname = newname.replace('MESSIER ', 'M', 1)
    if newname.startswith('M ') and is_number(newname[2:]):
        newname = newname.replace('M ', 'M', 1)
    if newname.startswith('M') and is_number(newname[1:]):
        newname = 'M' + newname[1:].lstrip(" 0")
    if len(newname) > 4 and newname.startswith("PGC "):
        newname = newname[:4] + newname[4:].lstrip(" 0")
    if len(newname) > 4 and newname.startswith("UGC "):
        newname = newname[:4] + newname[4:].lstrip(" 0")
    if len(newname) > 5 and newname.startswith(("MCG +", "MCG -")):
        newname = newname[:5] + '-'.join([x.zfill(2) for x in newname[5:].strip().split("-")])
    if len(newname) > 5 and newname.startswith("CGCG "):
        newname = newname[:5] + '-'.join([x.zfill(3) for x in newname[5:].strip().split("-")])
    if (len(newname) > 1 and newname.startswith("E")) or (len(newname) > 3 and newname.startswith('ESO')):
        if newname[0] == "E":
            esplit = newname[1:].split("-")
        else:
            esplit = newname[3:].split("-")
        if len(esplit) == 2 and is_number(esplit[0].strip()):
            if esplit[1].strip()[0] == 'G':
                parttwo = esplit[1][1:].strip()
            else:
                parttwo = esplit[1].strip()
            if is_number(parttwo.strip()):
                newname = 'ESO ' + esplit[0].lstrip('0') + '-G' + parttwo.lstrip('0')
    newname = ' '.join(newname.split())
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
        (mldt, mlmag, mlband, mlsource) = get_max_light(events, name)
        if mldt:
            source = events[name].add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
            events[name].add_quantity('maxdate', make_date_string(mldt.year, mldt.month, mldt.day),
                uniq_cdl([source]+mlsource.split(',')), derived = True)
        if mlmag:
            source = events[name].add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
            events[name].add_quantity('maxappmag', pretty_num(mlmag),
                uniq_cdl([source]+mlsource.split(',')), derived = True)
        if mlband:
            source = events[name].add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
            events[name].add_quantity('maxband', mlband,
                uniq_cdl([source]+mlsource.split(',')), derived = True)

    if 'discoverdate' not in events[name] or max([len(x['value'].split('/')) for x in events[name]['discoverdate']]) < 3:
        (fldt, flsource) = get_first_light(events, name)
        if fldt:
            source = events[name].add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
            events[name].add_quantity('discoverdate', make_date_string(fldt.year, fldt.month, fldt.day),
                uniq_cdl([source]+flsource.split(',')), derived = True)

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
            source = events[name].add_source(bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL, secondary=True)
            events[name].add_quantity('discoverdate', make_date_string(fldt.year, fldt.month, fldt.day),
                uniq_cdl([source]+minspecsource.split(',')), derived = True)


def clean_snname(string):
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
    return ','.join(sorted(list(set(values))))


def utf8(x):
    return str(x, 'utf-8')
