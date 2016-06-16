"""Import data from OGLE.
"""
# from astropy.time import Time as astrotime
from bs4 import BeautifulSoup
import calendar
# from math import floor
import os
import re
import urllib

from scripts import PATH
from .. funcs import add_event, add_source, add_quantity, \
    clean_snname, journal_events, load_cached_url, uniq_cdl, utf8
from ... utils import pbar, is_number


def do_asiago_photo(events, args, tasks):
    current_task = 'asiago: photometry'
    # response = urllib.request.urlopen('http://graspa.oapd.inaf.it/cgi-bin/sncat.php')
    path = os.path.abspath(os.path.join(PATH.REPO_EXTERNAL, 'asiago-cat.php'))
    response = urllib.request.urlopen('file://' + path)
    html = response.read().decode('utf-8')
    html = html.replace('\r', "")

    soup = BeautifulSoup(html, 'html5lib')
    table = soup.find('table')

    records = []
    for r, row in enumerate(table.findAll('tr')):
        if r == 0:
            continue

        col = row.findAll('td')
        records.append([utf8(x.renderContents()) for x in col])

    for record in pbar(records, current_task):
        if len(record) > 1 and record[1] != '':
            name = clean_snname('SN' + record[1]).strip('?')
            name = add_event(tasks, args, events, name)

            reference = 'Asiago Supernova Catalogue'
            refurl = 'http://graspa.oapd.inaf.it/cgi-bin/sncat.php'
            refbib = '1989A&AS...81..421B'
            source = add_source(
                events, name, refname=reference, url=refurl, bibcode=refbib, secondary=True)
            add_quantity(events, name, 'alias', name, source)

            year = re.findall(r'\d+', name)[0]
            add_quantity(events, name, 'discoverdate', year, source)

            hostname = record[2]
            hostra = record[3]
            hostdec = record[4]
            ra = record[5].strip(':')
            dec = record[6].strip(':')
            redvel = record[11].strip(':')
            discoverer = record[19]

            datestring = year

            monthday = record[18]
            if "*" in monthday:
                datekey = 'discover'
            else:
                datekey = 'max'

            if monthday.strip() != '':
                monthstr = ''.join(re.findall('[a-zA-Z]+', monthday))
                monthstr = str(list(calendar.month_abbr).index(monthstr))
                datestring = datestring + '/' + monthstr

                dayarr = re.findall(r'\d+', monthday)
                if dayarr:
                    daystr = dayarr[0]
                    datestring = datestring + '/' + daystr

            add_quantity(events, name, datekey + 'date', datestring, source)

            velocity = ''
            redshift = ''
            if redvel != '':
                if round(float(redvel)) == float(redvel):
                    velocity = int(redvel)
                else:
                    redshift = float(redvel)
                redshift = str(redshift)
                velocity = str(velocity)

            claimedtype = record[17].replace(':', '').replace('*', '').strip()

            if (hostname != ''):
                add_quantity(events, name, 'host', hostname, source)
            if (claimedtype != ''):
                add_quantity(events, name, 'claimedtype', claimedtype, source)
            if (redshift != ''):
                add_quantity(events, name, 'redshift', redshift, source, kind='host')
            if (velocity != ''):
                add_quantity(events, name, 'velocity', velocity, source, kind='host')
            if (hostra != ''):
                add_quantity(events, name, 'hostra', hostra, source, unit='nospace')
            if (hostdec != ''):
                add_quantity(events, name, 'hostdec', hostdec, source, unit='nospace')
            if (ra != ''):
                add_quantity(events, name, 'ra', ra, source, unit='nospace')
            if (dec != ''):
                add_quantity(events, name, 'dec', dec, source, unit='nospace')
            if (discoverer != ''):
                add_quantity(events, name, 'discoverer', discoverer, source)

    events = journal_events(tasks, args, events)
    return events


def do_asiago_spectra(events, args, tasks):
    current_task = 'asiago: spectra'
    html = load_cached_url(args, 'http://sngroup.oapd.inaf.it./cgi-bin/output_class.cgi?sn=1990',
                           os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'Asiago/spectra.html'))
    if not html:
        return events

    bs = BeautifulSoup(html, 'html5lib')
    trs = bs.findAll('tr')
    for tr in pbar(trs, current_task):
        tds = tr.findAll('td')
        name = ''
        host = ''
        # fitsurl = ''
        source = ''
        reference = ''
        for tdi, td in enumerate(tds):
            if tdi == 0:
                butt = td.find('button')
                if not butt:
                    break
                alias = butt.text.strip()
                alias = alias.replace('PSNJ', 'PSN J').replace('GAIA', 'Gaia')
            elif tdi == 1:
                name = td.text.strip().replace('PSNJ', 'PSN J').replace('GAIA', 'Gaia')
                if name.startswith('SN '):
                    name = 'SN' + name[3:]
                if not name:
                    name = alias
                if is_number(name[:4]):
                    name = 'SN' + name
                name = add_event(tasks, args, events, name)
                reference = 'Asiago Supernova Catalogue'
                refurl = 'http://graspa.oapd.inaf.it/cgi-bin/sncat.php'
                secondarysource = add_source(
                    events, name, refname=reference, url=refurl, secondary=True)
                add_quantity(events, name, 'alias', name, secondarysource)
                if alias != name:
                    add_quantity(events, name, 'alias', alias, secondarysource)
            elif tdi == 2:
                host = td.text.strip()
                if host == 'anonymous':
                    host = ''
            elif tdi == 3:
                discoverer = td.text.strip()
            elif tdi == 5:
                ra = td.text.strip()
            elif tdi == 6:
                dec = td.text.strip()
            elif tdi == 7:
                claimedtype = td.text.strip()
            elif tdi == 8:
                redshift = td.text.strip()
            # elif tdi == 9:
            #     epochstr = td.text.strip()
            #     if epochstr:
            #         mjd = (astrotime(epochstr[:4] + '-' + epochstr[4:6] + '-' +
            #                str(floor(float(epochstr[6:]))).zfill(2)).mjd +
            #                float(epochstr[6:]) - floor(float(epochstr[6:])))
            #     else:
            #         mjd = ''
            elif tdi == 10:
                refs = td.findAll('a')
                source = ''
                reference = ''
                refurl = ''
                for ref in refs:
                    if ref.text != 'REF':
                        reference = ref.text
                        refurl = ref['href']
                if reference:
                    source = add_source(events, name, refname=reference, url=refurl)
                add_quantity(events, name, 'alias', name, secondarysource)
                sources = uniq_cdl(list(filter(None, [source, secondarysource])))
            elif tdi == 12:
                pass
                # fitslink = td.find('a')
                # if fitslink:
                #     fitsurl = fitslink['href']
        if name:
            add_quantity(events, name, 'claimedtype', claimedtype, sources)
            add_quantity(events, name, 'ra', ra, sources)
            add_quantity(events, name, 'dec', dec, sources)
            add_quantity(events, name, 'redshift', redshift, sources)
            add_quantity(events, name, 'discoverer', discoverer, sources)
            add_quantity(events, name, 'host', host, sources)

            # if fitsurl:
            #    response = urllib.request.urlopen('http://sngroup.oapd.inaf.it./' + fitsurl)
            #    compressed = io.BytesIO(response.read())
            #    decompressed = gzip.GzipFile(fileobj=compressed)
            #    hdulist = fits.open(decompressed)
            #    scidata = hdulist[0].data
            #    print(hdulist[0].header)
            #
            #    print(scidata[3])
            #    sys.exit()

    events = journal_events(tasks, args, events)
    return events
