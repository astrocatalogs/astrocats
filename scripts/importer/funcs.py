"""Utility functions for OSC import.
"""

import json
import os
import statistics
import warnings
from collections import OrderedDict
from math import floor, hypot, log10, pi, sqrt

from astropy import units as un
from astropy.coordinates import SkyCoord as coord
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value

from cdecimal import Decimal

from ..utils import (get_sig_digits, is_number, pretty_num, round_sig, tprint,
                     zpad)
from .constants import (ADS_BIB_URL, CLIGHT, KM, OSC_BIBCODE, OSC_NAME,
                        OSC_URL, PREF_KINDS)


def alias_priority(name, attr):
    if name == attr:
        return 0
    return 1


def convert_aq_output(row):
    return OrderedDict([(x, str(row[x]) if is_number(row[x]) else row[x])
                        for x in row.colnames])


def derive_and_sanitize(catalog):
    # Calculate some columns based on imported data, sanitize some fields
    for name, event in catalog.events.items():
        aliases = event.get_aliases(includename=False)
        if name not in aliases:
            if 'sources' in event:
                event.add_quantity('alias', name, '1')
            else:
                source = event.add_source(
                    bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL,
                    secondary=True)
                event.add_quantity('alias', name, source)

        if ((name.startswith('SN') and is_number(name[2:6]) and
             'discoverdate' in event and
             int(event['discoverdate'][0]['value'].
                 split('/')[0]) >= 2016 and
             not any(['AT' in x for x in aliases]))):
            source = event.add_source(
                bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL,
                secondary=True)
            event.add_quantity('alias', 'AT' + name[2:], source)

        event['alias'] = list(
            sorted(event['alias'],
                   key=lambda key: alias_priority(name, key)))
        aliases = event.get_aliases()

        event.set_first_max_light()

        if 'claimedtype' in event:
            # FIX: this is something that should be done completely internally
            #      i.e. add it to `clean` or something??
            event['claimedtype'] = event.ct_list_prioritized()
        if 'discoverdate' not in event:
            prefixes = ['MLS', 'SSS', 'CSS', 'GRB ']
            for alias in aliases:
                for prefix in prefixes:
                    if (alias.startswith(prefix) and
                            is_number(alias.replace(prefix, '')[:2])):
                        discoverdate = ('/'.
                                        join(['20' +
                                              alias.replace(prefix, '')[:2],
                                              alias.replace(prefix, '')[2:4],
                                              alias.replace(prefix, '')[4:6]]))
                        if catalog.args.verbose:
                            tprint(
                                'Added discoverdate from name [' +
                                alias + ']: ' + discoverdate)
                        source = event.add_source(
                            bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL,
                            secondary=True)
                        event.add_quantity(
                            'discoverdate', discoverdate, source, derived=True)
                        break
                if 'discoverdate' in event:
                    break
        if 'discoverdate' not in event:
            prefixes = ['ASASSN-', 'PS1-', 'PS1', 'PS', 'iPTF', 'PTF', 'SCP-',
                        'SNLS-', 'SPIRITS', 'LSQ', 'DES', 'SNHiTS',
                        'GND', 'GNW', 'GSD', 'GSW', 'EGS', 'COS']
            for alias in aliases:
                for prefix in prefixes:
                    if (alias.startswith(prefix) and
                            is_number(alias.replace(prefix, '')[:2])):
                        discoverdate = '20' + alias.replace(prefix, '')[:2]
                        if catalog.args.verbose:
                            tprint(
                                'Added discoverdate from name [' +
                                alias + ']: ' + discoverdate)
                        source = event.add_source(
                            bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL,
                            secondary=True)
                        event.add_quantity(
                            'discoverdate', discoverdate, source, derived=True)
                        break
                if 'discoverdate' in event:
                    break
        if 'discoverdate' not in event:
            prefixes = ['SNF']
            for alias in aliases:
                for prefix in prefixes:
                    if (alias.startswith(prefix) and
                            is_number(alias.replace(prefix, '')[:4])):
                        discoverdate = ('/'
                                        .join(
                                            [alias.replace(prefix, '')[:4],
                                             alias.replace(prefix, '')[4:6],
                                             alias.replace(prefix, '')[6:8]]))
                        if catalog.args.verbose:
                            tprint(
                                'Added discoverdate from name [' +
                                alias + ']: ' + discoverdate)
                        source = event.add_source(
                            bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL,
                            secondary=True)
                        event.add_quantity(
                            'discoverdate', discoverdate, source, derived=True)
                        break
                if 'discoverdate' in event:
                    break
        if 'discoverdate' not in event:
            prefixes = ['PTFS', 'SNSDF']
            for alias in aliases:
                for prefix in prefixes:
                    if (alias.startswith(prefix) and
                            is_number(alias.replace(prefix, '')[:2])):
                        discoverdate = ('/'
                                        .join(
                                            ['20' +
                                             alias.replace(prefix, '')[:2],
                                             alias.replace(prefix, '')[2:4]]))
                        if catalog.args.verbose:
                            tprint(
                                'Added discoverdate from name [' +
                                alias + ']: ' + discoverdate)
                        source = event.add_source(
                            bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL,
                            secondary=True)
                        event.add_quantity(
                            'discoverdate', discoverdate, source, derived=True)
                        break
                if 'discoverdate' in event:
                    break
        if 'discoverdate' not in event:
            prefixes = ['AT', 'SN', 'OGLE-', 'SM ', 'KSN-']
            for alias in aliases:
                for prefix in prefixes:
                    if (alias.startswith(prefix) and
                            is_number(alias.replace(prefix, '')[:4]) and
                            '.' not in alias.replace(prefix, '')[:4]):
                        discoverdate = alias.replace(prefix, '')[:4]
                        if catalog.args.verbose:
                            tprint(
                                'Added discoverdate from name [' +
                                alias + ']: ' + discoverdate)
                        source = event.add_source(
                            bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL,
                            secondary=True)
                        event.add_quantity(
                            'discoverdate', discoverdate, source, derived=True)
                        break
                if 'discoverdate' in event:
                    break
        if 'ra' not in event or 'dec' not in event:
            prefixes = ['PSN J', 'MASJ', 'CSS', 'SSS', 'MASTER OT J', 'HST J',
                        'TCP J', 'MACS J', '2MASS J', 'EQ J', 'CRTS J',
                        'SMT J']
            for alias in aliases:
                for prefix in prefixes:
                    if (alias.startswith(prefix) and
                            is_number(alias.replace(prefix, '')[:6])):
                        noprefix = alias.split(
                            ':')[-1].replace(prefix, '').replace('.', '')
                        decsign = '+' if '+' in noprefix else '-'
                        noprefix = noprefix.replace('+', '|').replace('-', '|')
                        nops = noprefix.split('|')
                        if len(nops) < 2:
                            continue
                        rastr = nops[0]
                        decstr = nops[1]
                        ra = ':'.join([rastr[:2], rastr[2:4], rastr[4:6]]) + \
                            ('.' + rastr[6:] if len(rastr) > 6 else '')
                        dec = (decsign + ':'
                               .join([decstr[:2], decstr[2:4], decstr[4:6]]) +
                               ('.' + decstr[6:] if len(decstr) > 6 else ''))
                        if catalog.args.verbose:
                            tprint('Added ra/dec from name: ' + ra + ' ' + dec)
                        source = event.add_source(
                            bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL,
                            secondary=True)
                        event.add_quantity(
                            'ra', ra, source, derived=True)
                        event.add_quantity(
                            'dec', dec, source, derived=True)
                        break
                if 'ra' in event:
                    break

        no_host = ('host' not in event or
                   not any([x['value'] == 'Milky Way' for x in
                            event['host']]))
        if ('ra' in event and 'dec' in event and no_host):
            from astroquery.irsa_dust import IrsaDust
            if name not in catalog.extinctions_dict:
                try:
                    ra_dec = event['ra'][0]['value'] + \
                        " " + event['dec'][0]['value']
                    result = IrsaDust.get_query_table(ra_dec, section='ebv')
                except (KeyboardInterrupt, SystemExit):
                    raise
                except:
                    warnings.warn("Coordinate lookup for " +
                                  name + " failed in IRSA.")
                else:
                    ebv = result['ext SandF mean'][0]
                    ebverr = result['ext SandF std'][0]
                    catalog.extinctions_dict[name] = [ebv, ebverr]
            if name in catalog.extinctions_dict:
                sources = uniq_cdl(
                    [event.add_source(
                        bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL,
                        secondary=True),
                     event.add_source(bibcode='2011ApJ...737..103S')])
                (event
                 .add_quantity('ebv',
                               str(catalog
                                   .extinctions_dict[name][0]),
                               sources,
                               error=str(catalog
                                         .extinctions_dict[name][1]),
                               derived=True))
        if (('host' in event and
             ('hostra' not in event or 'hostdec' not in event))):
            for host in event['host']:
                alias = host['value']
                if ' J' in alias and is_number(alias.split(' J')[-1][:6]):
                    noprefix = alias.split(
                        ' J')[-1].split(':')[-1].replace('.', '')
                    decsign = '+' if '+' in noprefix else '-'
                    noprefix = noprefix.replace('+', '|').replace('-', '|')
                    nops = noprefix.split('|')
                    if len(nops) < 2:
                        continue
                    rastr = nops[0]
                    decstr = nops[1]
                    hostra = (':'.join([rastr[:2], rastr[2:4], rastr[4:6]]) +
                              ('.' + rastr[6:] if len(rastr) > 6 else ''))
                    hostdec = decsign + ':'.join([decstr[:2], decstr[2:4],
                                                  decstr[4:6]]) + (
                        '.' + decstr[6:] if len(decstr) > 6 else '')
                    if catalog.args.verbose:
                        tprint('Added hostra/hostdec from name: ' +
                               hostra + ' ' + hostdec)
                    source = event.add_source(
                        bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL,
                        secondary=True)
                    event.add_quantity(
                        'hostra', hostra, source, derived=True)
                    event.add_quantity(
                        'hostdec', hostdec, source, derived=True)
                    break
                if 'hostra' in event:
                    break
        if 'claimedtype' in event:
            event['claimedtype'][:] = [ct for ct in event[
                'claimedtype'] if (ct['value'] != '?' and ct['value'] != '-')]
            if not len(event['claimedtype']):
                del(event['claimedtype'])
        if 'claimedtype' not in event and name.startswith('AT'):
            source = event.add_source(
                bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL,
                secondary=True)
            event.add_quantity('claimedtype', 'Candidate', source)
        if 'redshift' not in event and 'velocity' in event:
            # Find the "best" velocity to use for this
            bestsig = 0
            for hv in event['velocity']:
                sig = get_sig_digits(hv['value'])
                if sig > bestsig:
                    besthv = hv['value']
                    bestsrc = hv['source']
                    bestsig = sig
            if bestsig > 0 and is_number(besthv):
                voc = float(besthv) * 1.e5 / CLIGHT
                source = event.add_source(
                    bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL,
                    secondary=True)
                sources = uniq_cdl([source] + bestsrc.split(','))
                (event
                 .add_quantity('redshift',
                               pretty_num(sqrt((1. + voc) / (1. - voc)) - 1.,
                                          sig=bestsig),
                               sources, kind='heliocentric',
                               derived=True))
        if ('redshift' not in event and len(catalog.nedd_dict) > 0 and
                'host' in event):
            reference = "NED-D"
            refurl = "http://ned.ipac.caltech.edu/Library/Distances/"
            for host in event['host']:
                if host['value'] in catalog.nedd_dict:
                    source = event.add_source(
                        bibcode='2015arXiv150201589P')
                    secondarysource = event.add_source(
                        srcname=reference, url=refurl, secondary=True)
                    meddist = statistics.median(
                        catalog.nedd_dict[host['value']])
                    redz = z_at_value(
                        cosmo.comoving_distance, float(meddist) * un.Mpc)
                    redshift = pretty_num(
                        redz, sig=get_sig_digits(str(meddist)))
                    event.add_quantity(
                        name, 'redshift', redshift,
                        uniq_cdl([source, secondarysource]),
                        kind='host', derived=True)
        if ('maxabsmag' not in event and 'maxappmag' in event and
                'lumdist' in event):
            # Find the "best" distance to use for this
            bestsig = 0
            for ld in event['lumdist']:
                sig = get_sig_digits(ld['value'])
                if sig > bestsig:
                    bestld = ld['value']
                    bestsrc = ld['source']
                    bestsig = sig
            if bestsig > 0 and is_number(bestld) and float(bestld) > 0.:
                source = event.add_source(
                    bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL,
                    secondary=True)
                sources = uniq_cdl([source] + bestsrc.split(','))
                # FIX: what's happening here?!
                pnum = (float(event['maxappmag'][0]['value']) -
                        5.0 * (log10(float(bestld) * 1.0e6) - 1.0))
                pnum = pretty_num(pnum, sig=bestsig)
                event.add_quantity(
                    'maxabsmag', pnum, sources, derived=True)
        if 'redshift' in event:
            # Find the "best" redshift to use for this
            bestz, bestkind, bestsig, bestsrc = event.get_best_redshift()
            if bestsig > 0:
                bestz = float(bestz)
                if 'velocity' not in event:
                    source = event.add_source(
                        bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL,
                        secondary=True)
                    # FIX: what's happening here?!
                    pnum = CLIGHT / KM * \
                        ((bestz + 1.)**2. - 1.) / ((bestz + 1.)**2. + 1.)
                    pnum = pretty_num(pnum, sig=bestsig)
                    event.add_quantity(
                        'velocity', pnum, source, kind=PREF_KINDS[bestkind])
                if bestz > 0.:
                    from astropy.cosmology import Planck15 as cosmo
                    if 'lumdist' not in event:
                        dl = cosmo.luminosity_distance(bestz)
                        sources = [
                            event.add_source(
                                bibcode=OSC_BIBCODE, srcname=OSC_NAME,
                                url=OSC_URL, secondary=True),
                            event.add_source(bibcode='2015arXiv150201589P')]
                        sources = uniq_cdl(sources + bestsrc.split(','))
                        event.add_quantity(
                            'lumdist', pretty_num(dl.value, sig=bestsig),
                            sources, kind=PREF_KINDS[bestkind],
                            derived=True)
                        if ('maxabsmag' not in event and
                                'maxappmag' in event):
                            source = event.add_source(
                                bibcode=OSC_BIBCODE, srcname=OSC_NAME,
                                url=OSC_URL, secondary=True)
                            pnum = pretty_num(
                                float(event['maxappmag'][0]['value']) -
                                5.0 * (log10(dl.to('pc').value) - 1.0),
                                sig=bestsig)
                            event.add_quantity(
                                'maxabsmag', pnum, sources, derived=True)
                    if 'comovingdist' not in event:
                        cd = cosmo.comoving_distance(bestz)
                        sources = [
                            event.add_source(
                                bibcode=OSC_BIBCODE, srcname=OSC_NAME,
                                url=OSC_URL,
                                secondary=True),
                            event.add_source(bibcode='2015arXiv150201589P')]
                        sources = uniq_cdl(sources + bestsrc.split(','))
                        event.add_quantity(
                            'comovingdist', pretty_num(cd.value, sig=bestsig),
                            sources, derived=True)
        if all([x in event for x in ['ra', 'dec', 'hostra', 'hostdec']]):
            # For now just using first coordinates that appear in entry
            try:
                c1 = coord(
                    ra=event['ra'][0]['value'], dec=event['dec'][0]['value'],
                    unit=(un.hourangle, un.deg))
                c2 = coord(
                    ra=event['hostra'][0]['value'],
                    dec=event['hostdec'][0]['value'],
                    unit=(un.hourangle, un.deg))
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                pass
            else:
                sources = uniq_cdl(
                    [event.add_source(
                        bibcode=OSC_BIBCODE, srcname=OSC_NAME, url=OSC_URL,
                        secondary=True)] +
                    event['ra'][0]['source'].split(',') +
                    event['dec'][0]['source'].split(',') +
                    event['hostra'][0]['source'].split(',') +
                    event['hostdec'][0]['source'].split(','))
                if 'hostoffsetang' not in event:
                    hosa = Decimal(hypot(c1.ra.degree - c2.ra.degree,
                                         c1.dec.degree - c2.dec.degree))
                    hosa = pretty_num(hosa * Decimal(3600.))
                    event.add_quantity(
                        'hostoffsetang', hosa, sources,
                        derived=True, unit='arcseconds')
                if ('comovingdist' in event and
                        'redshift' in event and
                        'hostoffsetdist' not in event):
                    offsetsig = get_sig_digits(
                        event['hostoffsetang'][0]['value'])
                    sources = uniq_cdl(sources.split(',') +
                                       (event['comovingdist']
                                        [0]['source']).split(',') +
                                       (event['redshift']
                                        [0]['source']).split(','))
                    (event
                     .add_quantity('hostoffsetdist',
                                   pretty_num(
                                       float(event['hostoffsetang']
                                             [0]['value']) /
                                       3600. * (pi / 180.) *
                                       float(event['comovingdist']
                                             [0]['value']) *
                                       1000. / (1.0 +
                                                float(event['redshift']
                                                      [0]['value'])),
                                       sig=offsetsig), sources))

        if 'photometry' in event:
            event['photometry'].sort(
                key=lambda x: ((float(x['time']) if isinstance(x['time'], str)
                                else min([float(y) for y in x['time']])) if
                               'time' in x else 0.0,
                               x['band'] if 'band' in x else '',
                               float(x['magnitude']) if
                               'magnitude' in x else ''))
        if ('spectra' in event and
                list(filter(None, ['time' in x
                                   for x in event['spectra']]))):
            event['spectra'].sort(key=lambda x: (
                float(x['time']) if 'time' in x else 0.0))
        if 'sources' in event:
            for source in event['sources']:
                if 'bibcode' in source:
                    import urllib
                    from html import unescape
                    # First sanitize the bibcode
                    if len(source['bibcode']) != 19:
                        source['bibcode'] = urllib.parse.unquote(
                            unescape(source['bibcode'])).replace('A.A.', 'A&A')
                    if source['bibcode'] in catalog.biberror_dict:
                        source['bibcode'] = \
                            catalog.biberror_dict[source['bibcode']]

                    if source['bibcode'] not in catalog.bibauthor_dict:
                        bibcode = source['bibcode']
                        adsquery = (ADS_BIB_URL +
                                    urllib.parse.quote(bibcode) +
                                    '&data_type=Custom&format=%253m%20%25(y)')
                        response = urllib.request.urlopen(adsquery)
                        html = response.read().decode('utf-8')
                        hsplit = html.split("\n")
                        if len(hsplit) > 5:
                            bibcodeauthor = hsplit[5]
                        else:
                            bibcodeauthor = ''

                        if not bibcodeauthor:
                            warnings.warn(
                                "Bibcode didn't return authors, not converting"
                                "this bibcode.")

                        catalog.bibauthor_dict[bibcode] = unescape(
                            bibcodeauthor).strip()

            for source in event['sources']:
                if ('bibcode' in source and
                        source['bibcode'] in catalog.bibauthor_dict and
                        catalog.bibauthor_dict[source['bibcode']]):
                    source['reference'] = catalog.bibauthor_dict[
                        source['bibcode']]
                    if 'name' not in source and source['bibcode']:
                        source['name'] = source['bibcode']
        if 'redshift' in event:
            event['redshift'] = list(
                sorted(event['redshift'], key=lambda key:
                       frame_priority(key)))
        if 'velocity' in event:
            event['velocity'] = list(
                sorted(event['velocity'], key=lambda key:
                       frame_priority(key)))
        if 'claimedtype' in event:
            event['claimedtype'] = event.ct_list_prioritized()

        # event = OrderedDict(
        #    sorted(event.items(), key=lambda key:
        #           event_attr_priority(key[0])))

    return catalog.extinctions_dict, catalog.bibauthor_dict


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


def frame_priority(attr):
    if 'kind' in attr:
        if attr['kind'] in PREF_KINDS:
            return PREF_KINDS.index(attr['kind'])
        else:
            return len(PREF_KINDS)
    return len(PREF_KINDS)


def read_json_dict(filename):
    # path = '../atels.json'
    if os.path.isfile(filename):
        with open(filename, 'r') as f:
            mydict = json.loads(f.read(), object_pairs_hook=OrderedDict)
    else:
        mydict = OrderedDict()
    return mydict


def read_json_arr(filename):
    if os.path.isfile(filename):
        with open(filename, 'r') as f:
            myarr = json.loads(f.read())
    else:
        myarr = []
    return myarr


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


def load_cached_url(args, current_task, url, filepath, timeout=120, write=True,
                    failhard=False):
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
            if (x.status_code == 500 or x.status_code == 307 or
                    x.status_code == 404):
                raise
        txt = response.text
        newmd5 = md5(txt.encode('utf-8')).hexdigest()
        # tprint(filemd5 + ": " + newmd5)
        if args.update and newmd5 == filemd5:
            tprint('Skipping file in "' + current_task +
                   '," local and remote copies identical [' + newmd5 + '].')
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
        raise ValueError(
            "At least the year must be specified when constructing date "
            "string")
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
        if (len(namesp[1]) == 4 and is_number(namesp[1]) and
                is_number(namesp[3])):
            newname = 'OGLE-' + namesp[1] + '-SN-' + namesp[3].zfill(3)
    if newname.startswith('SN SDSS'):
        newname = newname.replace('SN SDSS ', 'SDSS', 1)
    if newname.startswith('SDSS '):
        newname = newname.replace('SDSS ', 'SDSS', 1)
    if newname.startswith('SDSS'):
        namesp = newname.split('-')
        if (len(namesp) == 3 and is_number(namesp[0][4:]) and
                is_number(namesp[1]) and is_number(namesp[2])):
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
    if (newname.startswith('SNF') and
            is_number(newname[3:]) and len(newname) >= 12):
        newname = 'SNF' + newname[3:11] + '-' + newname[11:]
    if newname.startswith(('MASTER OT J', 'ROTSE3 J')):
        prefix = newname.split('J')[0]
        coords = newname.split('J')[-1].strip()
        decsign = '+' if '+' in coords else '-'
        coordsplit = coords.replace('+', '-').split('-')
        if ('.' not in coordsplit[0] and
                len(coordsplit[0]) > 6 and '.' not in coordsplit[1] and
                len(coordsplit[1]) > 6):
            newname = (prefix + 'J' + coordsplit[0][:6] + '.' +
                       coordsplit[0][6:] + decsign + coordsplit[1][:6] +
                       '.' + coordsplit[1][6:])
    if (newname.startswith('Gaia ') and
            is_number(newname[3:4]) and len(newname) > 5):
        newname = newname.replace('Gaia ', 'Gaia', 1)
    if len(newname) <= 4 and is_number(newname):
        newname = 'SN' + newname + 'A'
    if (len(newname) > 4 and is_number(newname[:4]) and not
            is_number(newname[4:])):
        newname = 'SN' + newname
    if (newname.startswith('Sn ') and
            is_number(newname[3:7]) and len(newname) > 7):
        newname = newname.replace('Sn ', 'SN', 1)
    if (newname.startswith('sn') and
            is_number(newname[2:6]) and len(newname) > 6):
        newname = newname.replace('sn', 'SN', 1)
    if (newname.startswith('SN ') and
            is_number(newname[3:7]) and len(newname) > 7):
        newname = newname.replace('SN ', 'SN', 1)
    if (newname.startswith('SN') and
            is_number(newname[2:6]) and len(newname) == 7 and
            newname[6].islower()):
        newname = 'SN' + newname[2:6] + newname[6].upper()
    elif (newname.startswith('SN') and is_number(newname[2:6]) and
          (len(newname) == 8 or len(newname) == 9) and newname[6:].isupper()):
        newname = 'SN' + newname[2:6] + newname[6:].lower()

    newname = (' '.join(newname.split())).strip()
    return newname


def radec_clean(svalue, quantity, unit=''):
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
            svalue = str(hours).zfill(2) + ':' + str(minutes).zfill(2) + \
                ':' + zpad(pretty_num(seconds, sig=sig - 1))
        elif 'dec' in quantity:
            fldeg = abs(deg)
            degree = floor(fldeg)
            minutes = floor((fldeg - degree) * 60.0)
            seconds = (fldeg * 60.0 - (degree * 60.0 + minutes)) * 60.0
            if seconds > 60.0:
                raise(ValueError('Invalid seconds value for ' + quantity))
            svalue = (('+' if deg >= 0.0 else '-') +
                      str(degree).strip('+-').zfill(2) + ':' +
                      str(minutes).zfill(2) + ':' +
                      zpad(pretty_num(seconds, sig=sig - 1)))
    elif unit == 'nospace' and 'ra' in quantity:
        svalue = svalue[:2] + ':' + svalue[2:4] + \
            ((':' + zpad(svalue[4:])) if len(svalue) > 4 else '')
    elif unit == 'nospace' and 'dec' in quantity:
        if svalue.startswith(('+', '-')):
            svalue = svalue[:3] + ':' + svalue[3:5] + \
                ((':' + zpad(svalue[5:])) if len(svalue) > 5 else '')
        else:
            svalue = '+' + svalue[:2] + ':' + svalue[2:4] + \
                ((':' + zpad(svalue[4:])) if len(svalue) > 4 else '')
    else:
        svalue = svalue.replace(' ', ':')
        if 'dec' in quantity:
            valuesplit = svalue.split(':')
            svalue = (('-' if valuesplit[0].startswith('-') else '+') +
                      valuesplit[0].strip('+-').zfill(2) +
                      (':' + valuesplit[1].zfill(2) if
                       len(valuesplit) > 1 else '') +
                      (':' + zpad(valuesplit[2]) if
                       len(valuesplit) > 2 else ''))

    if 'ra' in quantity:
        sunit = 'hours'
    elif 'dec' in quantity:
        sunit = 'degrees'

    # Correct case of arcseconds = 60.0.
    valuesplit = svalue.split(':')
    if len(valuesplit) == 3 and valuesplit[-1] in ["60.0", "60.", "60"]:
        svalue = valuesplit[0] + ':' + str(Decimal(valuesplit[1]) +
                                           Decimal(1.0)) + ':' + "00.0"

    # Strip trailing dots.
    svalue = svalue.rstrip('.')

    return (svalue, sunit)


def host_clean(name):
    newname = name.strip(' ;,*')

    # Handle some special cases
    hostcases = {'M051a': 'M51A', 'M051b': 'M51B'}
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
        newname = newname[:5] + '-'.join([x.zfill(2)
                                          for x in
                                          newname[5:].strip().split("-")])
    if len(newname) > 5 and newname.startswith("CGCG "):
        newname = newname[:5] + '-'.join([x.zfill(3)
                                          for x in
                                          newname[5:].strip().split("-")])
    if ((len(newname) > 1 and newname.startswith("E")) or
            (len(newname) > 3 and newname.startswith('ESO'))):
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
                newname = 'ESO ' + \
                    esplit[0].lstrip('0') + '-G' + parttwo.lstrip('0')
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
           ((isinstance(photo[tag], str) and isinstance(val, str) and
             Decimal(photo[tag]) == Decimal(val)) or
            (isinstance(photo[tag], list) and isinstance(val, list) and
             photo[tag] == val))))))
    return issame


def same_tag_str(photo, val, tag):
    issame = ((tag not in photo and not val) or (
        tag in photo and not val) or (tag in photo and photo[tag] == val))
    return issame


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
            (len(x) > length and
             len(str(round_sig(float(x), length))) < len(x))
            else x for x in arr]


def uniq_cdl(values):
    return ','.join(sorted(list(set(values))))


def utf8(x):
    return str(x, 'utf-8')
