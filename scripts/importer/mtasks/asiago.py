"""Import data from OGLE.
"""
import calendar
# from math import floor
import os
import re
import urllib

# from astropy.time import Time as astrotime
from bs4 import BeautifulSoup

from scripts import PATH

from ...utils import is_number, pbar
from ..funcs import clean_snname, load_cached_url, uniq_cdl, utf8


def do_asiago_photo(catalog):
    current_task = catalog.current_task
    # response = (urllib.request
    # .urlopen('http://graspa.oapd.inaf.it/cgi-bin/sncat.php'))
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
            oldname = clean_snname("SN" + record[1]).strip('?')

            reference = 'Asiago Supernova Catalogue'
            refurl = 'http://graspa.oapd.inaf.it/cgi-bin/sncat.php'
            refbib = '1989A&AS...81..421B'

            (name,
             source) = catalog.new_event(oldname,
                                         srcname=reference, url=refurl,
                                         bibcode=refbib, secondary=True)

            year = re.findall(r'\d+', oldname)[0]
            catalog.events[name].add_quantity('discoverdate', year, source)

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

            catalog.events[name].add_quantity(datekey + 'date', datestring,
                                              source)

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
                catalog.events[name].add_quantity('host', hostname, source)
            if (claimedtype != ''):
                catalog.events[name].add_quantity('claimedtype', claimedtype,
                                                  source)
            if (redshift != ''):
                catalog.events[name].add_quantity(
                    'redshift', redshift, source, kind='host')
            if (velocity != ''):
                catalog.events[name].add_quantity(
                    'velocity', velocity, source, kind='host')
            if (hostra != ''):
                catalog.events[name].add_quantity(
                    'hostra', hostra, source, unit='nospace')
            if (hostdec != ''):
                catalog.events[name].add_quantity(
                    'hostdec', hostdec, source, unit='nospace')
            if (ra != ''):
                catalog.events[name].add_quantity('ra', ra, source,
                                                  unit='nospace')
            if (dec != ''):
                catalog.events[name].add_quantity('dec', dec, source,
                                                  unit='nospace')
            if (discoverer != ''):
                catalog.events[name].add_quantity('discoverer', discoverer,
                                                  source)

    catalog.journal_events()
    return


def do_asiago_spectra(catalog):
    current_task = catalog.current_task
    html = load_cached_url(catalog.args, current_task,
                           ('http://sngroup.oapd.inaf.it./'
                            'cgi-bin/output_class.cgi?sn=1990'),
                           os.path.join(PATH.REPO_EXTERNAL_SPECTRA,
                                        'Asiago/spectra.html'))
    if not html:
        return

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
                name = (td.text.strip()
                        .replace('PSNJ', 'PSN J')
                        .replace('GAIA', 'Gaia'))
                if name.startswith('SN '):
                    name = 'SN' + name[3:]
                if not name:
                    name = alias
                if is_number(name[:4]):
                    name = 'SN' + name
                oldname = name
                name = catalog.add_event(name)
                reference = 'Asiago Supernova Catalogue'
                refurl = 'http://graspa.oapd.inaf.it/cgi-bin/sncat.php'
                secondarysource = catalog.events[name].add_source(
                    srcname=reference, url=refurl, secondary=True)
                catalog.events[name].add_quantity('alias', oldname,
                                                  secondarysource)
                if alias != name:
                    catalog.events[name].add_quantity('alias', alias,
                                                      secondarysource)
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
            #         mjd = (astrotime(epochstr[:4] + '-' + epochstr[4:6] +
            #                '-' +
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
                    source = catalog.events[name].add_source(
                        srcname=reference, url=refurl)
                catalog.events[name].add_quantity('alias', name,
                                                  secondarysource)
                sources = uniq_cdl(
                    list(filter(None, [source, secondarysource])))
            elif tdi == 12:
                pass
                # fitslink = td.find('a')
                # if fitslink:
                #     fitsurl = fitslink['href']
        if name:
            catalog.events[name].add_quantity('claimedtype', claimedtype,
                                              sources)
            catalog.events[name].add_quantity('ra', ra, sources)
            catalog.events[name].add_quantity('dec', dec, sources)
            catalog.events[name].add_quantity('redshift', redshift, sources)
            catalog.events[name].add_quantity('discoverer', discoverer,
                                              sources)
            catalog.events[name].add_quantity('host', host, sources)

            # if fitsurl:
            #    response = urllib.request.urlopen(
            #        'http://sngroup.oapd.inaf.it./' + fitsurl)
            #    compressed = io.BytesIO(response.read())
            #    decompressed = gzip.GzipFile(fileobj=compressed)
            #    hdulist = fits.open(decompressed)
            #    scidata = hdulist[0].data
            #    print(hdulist[0].header)
            #
            #    print(scidata[3])
            #    sys.exit()

    catalog.journal_events()
    return
