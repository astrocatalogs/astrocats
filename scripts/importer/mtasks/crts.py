"""General data import tasks.
"""
import os
import re
import urllib

from bs4 import BeautifulSoup

from cdecimal import Decimal
from scripts import PATH
from scripts.utils import is_number, pbar

from .. import Events
from ..constants import TRAVIS_QUERY_LIMIT
from ..funcs import add_photometry, load_cached_url


def do_crts(catalog):
    crtsnameerrors = ['2011ax']
    current_task = catalog.current_task
    folders = ['catalina', 'MLS', 'SSS']
    for fold in pbar(folders, current_task):
        html = load_cached_url(catalog.args, current_task,
                               'http://nesssi.cacr.caltech.edu/' + fold +
                               '/AllSN.html',
                               os.path.join(PATH.REPO_EXTERNAL, 'CRTS', fold +
                                            '.html'))
        if not html:
            continue
        bs = BeautifulSoup(html, 'html5lib')
        trs = bs.findAll('tr')
        for tri, tr in enumerate(pbar(trs, current_task)):
            tds = tr.findAll('td')
            if not tds:
                continue
            # refs = []
            aliases = []
            crtsname = ''
            ra = ''
            dec = ''
            lclink = ''
            # ttype = ''
            # ctype = ''
            for tdi, td in enumerate(tds):
                if tdi == 0:
                    crtsname = td.contents[0].text.strip()
                elif tdi == 1:
                    ra = td.contents[0]
                elif tdi == 2:
                    dec = td.contents[0]
                elif tdi == 11:
                    lclink = td.find('a')['onclick']
                    lclink = lclink.split("'")[1]
                elif tdi == 13:
                    aliases = re.sub('[()]', '', re.sub(
                        '<[^<]+?>', '', td.contents[-1].strip()))
                    aliases = [xx.strip('; ') for xx in list(
                        filter(None, aliases.split(' ')))]

            name = ''
            hostmag = ''
            hostupper = False
            validaliases = []
            for ai, alias in enumerate(aliases):
                if alias in ['SN', 'SDSS']:
                    continue
                if alias in crtsnameerrors:
                    continue
                if alias == 'mag':
                    if ai < len(aliases) - 1:
                        ind = ai + 1
                        if aliases[ai + 1] in ['SDSS']:
                            ind = ai + 2
                        elif aliases[ai + 1] in ['gal', 'obj', 'object',
                                                 'source']:
                            ind = ai - 1
                        if '>' in aliases[ind]:
                            hostupper = True
                        hostmag = aliases[ind].strip('>~').replace(
                            ',', '.').replace('m', '.')
                    continue
                if (is_number(alias[:4]) and alias[:2] == '20' and
                        len(alias) > 4):
                    name = 'SN' + alias
                if ((('asassn' in alias and len(alias) > 6) or
                     ('ptf' in alias and len(alias) > 3) or
                     ('ps1' in alias and len(alias) > 3) or
                     'snhunt' in alias or
                     ('mls' in alias and len(alias) > 3) or
                     'gaia' in alias or
                     ('lsq' in alias and len(alias) > 3))):
                    alias = alias.replace('SNHunt', 'SNhunt')
                    validaliases.append(alias)

            if not name:
                name = crtsname
            name = catalog.add_event(name)
            source = catalog.events[name].add_source(
                srcname='Catalina Sky Survey', bibcode='2009ApJ...696..870D',
                url='http://nesssi.cacr.caltech.edu/catalina/AllSN.html')
            catalog.events[name].add_quantity('alias', name, source)
            for alias in validaliases:
                catalog.events[name].add_quantity('alias', alias, source)
            catalog.events[name].add_quantity('ra', ra, source, unit='floatdegrees')
            catalog.events[name].add_quantity('dec', dec, source, unit='floatdegrees')

            if hostmag:
                # 1.0 magnitude error based on Drake 2009 assertion that SN are
                # only considered
                #    real if they are 2 mags brighter than host.
                add_photometry(catalog.events, name, band='C', magnitude=hostmag,
                               e_magnitude=1.0, source=source,
                               host=True, telescope='Catalina Schmidt',
                               upperlimit=hostupper)

            fname2 = (PATH.REPO_EXTERNAL + '/' + fold + '/' +
                      lclink.split('.')[-2].rstrip('p').split('/')[-1] +
                      '.html')
            if task_obj.load_archive(args) and os.path.isfile(fname2):
                with open(fname2, 'r') as ff:
                    html2 = ff.read()
            else:
                with open(fname2, 'w') as ff:
                    response2 = urllib.request.urlopen(lclink)
                    html2 = response2.read().decode('utf-8')
                    ff.write(html2)

            lines = html2.splitlines()
            teles = 'Catalina Schmidt'
            for line in lines:
                if 'javascript:showx' in line:
                    mjdstr = re.search("showx\('(.*?)'\)",
                                       line).group(1).split('(')[0].strip()
                    if not is_number(mjdstr):
                        continue
                    mjd = str(Decimal(mjdstr) + Decimal(53249.0))
                else:
                    continue
                if 'javascript:showy' in line:
                    mag = re.search("showy\('(.*?)'\)", line).group(1)
                if 'javascript:showz' in line:
                    err = re.search("showz\('(.*?)'\)", line).group(1)
                e_mag = err if float(err) > 0.0 else ''
                upl = (float(err) == 0.0)
                add_photometry(catalog.events, name, time=mjd, band='C', magnitude=mag,
                               source=source,
                               includeshost=True, telescope=teles,
                               e_magnitude=e_mag, upperlimit=upl)
            if args.update:
                events = Events.journal_events(
                    tasks, args, events, log)

        if args.travis and tri > TRAVIS_QUERY_LIMIT:
            break

    catalog.journal_events()
    return
