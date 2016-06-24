"""Imports from the 'Vizier' catalog.
"""
from cdecimal import Decimal
import csv
import os
from math import isnan

from astroquery.vizier import Vizier
from astropy.time import Time as astrotime

from scripts import PATH
from .. import Events
from .. funcs import add_photometry, convert_aq_output, jd_to_mjd, make_date_string, uniq_cdl
from .. constants import KM, CLIGHT
from ... utils import get_sig_digits, is_number, pbar, pretty_num, round_sig


def do_vizier(events, stubs, args, tasks, task_obj, log):
    """
    """
    current_task = task_obj.current_task(args)

    Vizier.ROW_LIMIT = -1

    # 2012ApJS..200...12H
    result = Vizier.get_catalogs('J/ApJS/200/12/table1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    oldname = ''
    for row in pbar(table, current_task):
        name = row['SN']
        if is_number(name[:4]):
            name = 'SN' + name
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2012ApJS..200...12H')
        events[name].add_quantity('alias', name, source)
        if '[' not in row['Gal']:
            events[name].add_quantity('host', row['Gal'].replace('_', ' '), source)
        events[name].add_quantity('redshift', str(row['z']), source, kind='heliocentric')
        events[name].add_quantity('redshift', str(row['zCMB']), source, kind='cmb')
        events[name].add_quantity('ebv', str(row['E_B-V_']), source, error=str(row['e_E_B-V_']) if row['e_E_B-V_'] else '')
        events[name].add_quantity('ra', row['RAJ2000'], source)
        events[name].add_quantity('dec', row['DEJ2000'], source)

    # 2012ApJ...746...85S
    result = Vizier.get_catalogs('J/ApJ/746/85/table1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    oldname = ''
    for row in pbar(table, current_task):
        name = row['Name'].replace('SCP', 'SCP-')
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2012ApJ...746...85S')
        events[name].add_quantity('alias', name, source)
        if row['f_Name']:
            events[name].add_quantity('claimedtype', 'Ia', source)
        if row['z']:
            events[name].add_quantity('redshift', str(row['z']), source, kind='spectroscopic')
        else:
            events[name].add_quantity('redshift', str(row['zCl']), source, kind='cluster')
        events[name].add_quantity('ebv', str(row['E_B-V_']), source)
        events[name].add_quantity('ra', row['RAJ2000'], source)
        events[name].add_quantity('dec', row['DEJ2000'], source)

    result = Vizier.get_catalogs('J/ApJ/746/85/table2')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    oldname = ''
    for row in pbar(table, current_task):
        name = row['Name'].replace('SCP', 'SCP-')
        flux = Decimal(float(row['Flux']))
        if flux <= 0.0:
            continue
        err = Decimal(float(row['e_Flux']))
        zp = Decimal(float(row['Zero']))
        sig = get_sig_digits(str(row['Flux']))+1
        magnitude = pretty_num(zp-Decimal(2.5)*(flux.log10()), sig=sig)
        e_magnitude = pretty_num(Decimal(2.5)*(Decimal(1.0) + err/flux).log10(), sig=sig)
        if float(e_magnitude) > 5.0:
            continue
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2012ApJ...746...85S')
        events[name].add_quantity('alias', name, source)
        add_photometry(
            events, name, time=str(row['MJD']), band=row['Filter'], instrument=row['Inst'],
            magnitude=magnitude, e_magnitude=e_magnitude, source=source)

    # 2004ApJ...602..571B
    result = Vizier.get_catalogs('J/ApJ/602/571/table8')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    oldname = ''
    for row in pbar(table, current_task):
        name = 'SN'+row['SN']
        flux = Decimal(float(row['Flux']))
        if flux <= 0.0:
            continue
        err = Decimal(float(row['e_Flux']))
        sig = get_sig_digits(str(row['Flux']))+1
        magnitude = pretty_num(Decimal(25.0)-Decimal(2.5)*(flux.log10()), sig=sig)
        e_magnitude = pretty_num(Decimal(2.5)*(Decimal(1.0) + err/flux).log10(), sig=sig)
        if float(e_magnitude) > 5.0:
            continue
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2004ApJ...602..571B')
        events[name].add_quantity('alias', name, source)
        band = row['Filt']
        system = ''
        telescope = ''
        if band in ['R', 'I']:
            system = 'Cousins'
        if band == 'Z':
            telescope = 'Subaru'
        add_photometry(
            events, name, time=str(row['MJD']), band=band, system=system, telescope=telescope,
            magnitude=magnitude, e_magnitude=e_magnitude, source=source)

    # 2014MNRAS.444.3258M
    result = Vizier.get_catalogs('J/MNRAS/444/3258/SNe')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    oldname = ''
    for row in pbar(table, current_task):
        name = row['SN']
        if name == oldname:
            continue
        oldname = name
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2014MNRAS.444.3258M')
        events[name].add_quantity('alias', name, source)
        events[name].add_quantity('redshift', str(row['z']), source, kind='heliocentric', error=str(row['e_z']))
        events[name].add_quantity('ra', str(row['_RA']), source, unit='floatdegrees')
        events[name].add_quantity('dec', str(row['_DE']), source, unit='floatdegrees')
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2014MNRAS.438.1391P
    result = Vizier.get_catalogs('J/MNRAS/438/1391/table2')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        name = row['SN']
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2014MNRAS.438.1391P')
        events[name].add_quantity('alias', name, source)
        events[name].add_quantity('redshift', str(row['zh']), source, kind='heliocentric')
        events[name].add_quantity('ra', row['RAJ2000'], source)
        events[name].add_quantity('dec', row['DEJ2000'], source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2012ApJ...749...18B
    result = Vizier.get_catalogs('J/ApJ/749/18/table1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        name = row['Name'].replace(' ', '')
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2012ApJ...749...18B')
        events[name].add_quantity('alias', name, source)
        mjd = str(astrotime(2450000.+row['JD'], format='jd').mjd)
        band = row['Filt']
        magnitude = str(row['mag'])
        e_magnitude = str(row['e_mag'])
        e_magnitude = '' if e_magnitude == '--' else e_magnitude
        upperlimit = True if row['l_mag'] == '>' else False
        add_photometry(
            events, name, time=mjd, band=band, magnitude=magnitude, e_magnitude=e_magnitude, instrument='UVOT',
            source=source, upperlimit=upperlimit, telescope='Swift', system='Swift')
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2010A&A...523A...7G
    result = Vizier.get_catalogs('J/A+A/523/A7/table9')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        name = 'SNLS-' + row['SNLS']
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2010A&A...523A...7G')
        events[name].add_quantity('alias', name, source)
        astrot = astrotime(2450000.+row['Date1'], format='jd').datetime
        events[name].add_quantity('discoverdate', make_date_string(astrot.year, astrot.month, astrot.day), source)
        events[name].add_quantity('ebv', str(row['E_B-V_']), source)
        events[name].add_quantity('redshift', str(row['z']), source, kind='heliocentric')
        type_str = row['Type'].replace('*', '?').replace('SN', '').replace('(pec)', ' P').replace('Ia? P?', 'Ia P?')
        events[name].add_quantity('claimedtype', type_str, source)
        events[name].add_quantity('ra', row['RAJ2000'], source)
        events[name].add_quantity('dec', row['DEJ2000'], source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2004A&A...415..863G
    result = Vizier.get_catalogs('J/A+A/415/863/table1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        name = 'SN' + row['SN']
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2004A&A...415..863G')
        events[name].add_quantity('alias', name, source)
        datesplit = row['Date'].split('-')
        date_str = make_date_string(datesplit[0], datesplit[1].lstrip('0'), datesplit[2].lstrip('0'))
        events[name].add_quantity('discoverdate', date_str, source)
        events[name].add_quantity('host', 'Abell ' + str(row['Abell']), source)
        events[name].add_quantity('claimedtype', row['Type'], source)
        events[name].add_quantity('ra', row['RAJ2000'], source)
        events[name].add_quantity('dec', row['DEJ2000'], source)
        if row['zSN']:
            events[name].add_quantity('redshift', str(row['zSN']), source, kind='spectroscopic')
        else:
            events[name].add_quantity('redshift', str(row['zCl']), source, kind='cluster')
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2008AJ....136.2306H
    result = Vizier.get_catalogs('J/AJ/136/2306/sources')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        name = 'SDSS-II ' + str(row['SNID'])
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2008AJ....136.2306H')
        events[name].add_quantity('alias', name, source)
        events[name].add_quantity('claimedtype', row['SpType'].replace('SN.', '').strip(':'), source)
        events[name].add_quantity('ra', row['RAJ2000'], source)
        events[name].add_quantity('dec', row['DEJ2000'], source)

    # 2010ApJ...708..661D
    result = Vizier.get_catalogs('J/ApJ/708/661/sn')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        name = row['SN']
        if not name:
            name = 'SDSS-II ' + str(row['SDSS-II'])
        else:
            name = 'SN' + name
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2010ApJ...708..661D')
        events[name].add_quantity('alias', name, source)
        events[name].add_quantity('alias', 'SDSS-II ' + str(row['SDSS-II']), source)
        events[name].add_quantity('claimedtype', 'II P', source)
        events[name].add_quantity('ra', row['RAJ2000'], source)
        events[name].add_quantity('dec', row['DEJ2000'], source)

    result = Vizier.get_catalogs('J/ApJ/708/661/table1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        if row['f_SN'] == 'a':
            name = 'SDSS-II ' + str(row['SN'])
        else:
            name = 'SN' + row['SN']
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2010ApJ...708..661D')
        events[name].add_quantity('alias', name, source)
        events[name].add_quantity('redshift', str(row['z']), source, error=str(row['e_z']))
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2014ApJ...795...44R
    result = Vizier.get_catalogs('J/ApJ/795/44/ps1_snIa')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        name = row['SN']
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2014ApJ...795...44R')
        events[name].add_quantity('alias', name, source)
        astrot = astrotime(row['tdisc'], format='mjd').datetime
        events[name].add_quantity('discoverdate',  make_date_string(astrot.year, astrot.month, astrot.day), source)
        events[name].add_quantity('redshift', str(row['z']), source, error=str(row['e_z']), kind='heliocentric')
        events[name].add_quantity('ra', row['RAJ2000'], source)
        events[name].add_quantity('dec', row['DEJ2000'], source)
        events[name].add_quantity('claimedtype', 'Ia', source)

    result = Vizier.get_catalogs('J/ApJ/795/44/table6')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        name = row['SN']
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2014ApJ...795...44R')
        events[name].add_quantity('alias', name, source)
        if row['mag'] != '--':
            add_photometry(
                events, name, time=str(row['MJD']), band=row['Filt'], magnitude=str(row['mag']),
                e_magnitude=str(row['e_mag']), source=source, system='AB', telescope='PS1', instrument='PS1')
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 1990A&AS...82..145C
    result = Vizier.get_catalogs('II/189/mag')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)

    with open(os.path.join(PATH.REPO_EXTERNAL, 'II_189_refs.csv')) as f:
        tsvin = csv.reader(f, delimiter='\t', skipinitialspace=True)
        ii189bibdict = {}
        ii189refdict = {}
        for r, row in enumerate(tsvin):
            if row[0] != '0':
                ii189bibdict[r+1] = row[1]
            else:
                ii189refdict[r+1] = row[2]

    for row in pbar(table, current_task):
        if row['band'][0] == '(':
            continue
        name = 'SN' + row['SN']
        events, name = Events.add_event(tasks, args, events, name, log)
        source = ''
        secsource = events[name].add_source(bibcode='1990A&AS...82..145C', secondary=True)
        mjd = str(jd_to_mjd(Decimal(row['JD'])))
        mag = str(row['m'])
        band = row['band'].strip("'")
        if row['r_m'] in ii189bibdict:
            source = events[name].add_source(bibcode=ii189bibdict[row['r_m']])
        else:
            source = events[name].add_source(srcname=ii189refdict[row['r_m']])
        events[name].add_quantity('alias', name, source)

        add_photometry(events, name, time=mjd, band=band, magnitude=mag, source=uniq_cdl([source, secsource]))
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2014yCat.7272....0G
    result = Vizier.get_catalogs('VII/272/snrs')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    this_url = 'https://www.mrao.cam.ac.uk/surveys/snrs/snrs.data.html'
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        name = ''
        if row['Names']:
            names = row['Names'].split(',')
            for nam in names:
                if nam.strip()[:2] == 'SN':
                    name = nam.strip()
                    if is_number(name[2:]):
                        name = name + 'A'
            if not name:
                for nam in names:
                    if nam.strip('()') == nam:
                        name = nam.strip()
                        break
        if not name:
            name = row['SNR'].strip()

        events, name = Events.add_event(tasks, args, events, name, log)
        source = (events[name].add_source(bibcode='2014BASI...42...47G') + ',' +
                  events[name].add_source(srcname='Galactic SNRs', url=this_url))
        events[name].add_quantity('alias', name, source)

        events[name].add_quantity('alias', row['SNR'].strip(), source)
        events[name].add_quantity('alias', 'MWSNR '+row['SNR'].strip('G '), source)

        if row['Names']:
            names = row['Names'].split(',')
            for nam in names:
                events[name].add_quantity('alias', nam.replace('Vela (XYZ)', 'Vela').strip('()').strip(), source)
                if nam.strip()[:2] == 'SN':
                    events[name].add_quantity('discoverdate', nam.strip()[2:], source)

        events[name].add_quantity('host', 'Milky Way', source)
        events[name].add_quantity('ra', row['RAJ2000'], source)
        events[name].add_quantity('dec', row['DEJ2000'], source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2014MNRAS.442..844F
    result = Vizier.get_catalogs('J/MNRAS/442/844/table1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        name = 'SN' + row['SN']
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2014MNRAS.442..844F')
        events[name].add_quantity('alias', name, source)
        events[name].add_quantity('redshift', str(row['zhost']), source, kind='host')
        events[name].add_quantity('ebv', str(row['E_B-V_']), source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    result = Vizier.get_catalogs('J/MNRAS/442/844/table2')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    instr = 'KAIT'
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        name = 'SN' + str(row['SN'])
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2014MNRAS.442..844F')
        events[name].add_quantity('alias', name, source)
        for band in ['B', 'V', 'R', 'I']:
            bandtag = band + 'mag'
            if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                add_photometry(
                    events, name, time=row['MJD'], band=band, magnitude=row[bandtag],
                    e_magnitude=row['e_' + bandtag], source=source, telescope=instr,
                    instrument=instr)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2012MNRAS.425.1789S
    result = Vizier.get_catalogs('J/MNRAS/425/1789/table1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        name = ''.join(row['SimbadName'].split(' '))
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2012MNRAS.425.1789S')
        events[name].add_quantity('alias', name, source)
        events[name].add_quantity('alias', 'SN' + row['SN'], source)
        events[name].add_quantity('host', row['Gal'], source)
        if is_number(row['cz']):
            red_str = str(round_sig(float(row['cz'])*KM/CLIGHT, sig=get_sig_digits(str(row['cz']))))
            events[name].add_quantity('redshift', red_str, source, kind='heliocentric')
        events[name].add_quantity('ebv', str(row['E_B-V_']), source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2015ApJS..219...13W
    result = Vizier.get_catalogs('J/ApJS/219/13/table3')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        name = u'LSQ' + str(row['LSQ'])
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2015ApJS..219...13W')
        events[name].add_quantity('alias', name, source)
        events[name].add_quantity('ra', row['RAJ2000'], source)
        events[name].add_quantity('dec', row['DEJ2000'], source)
        events[name].add_quantity('redshift', row['z'], source,
                                  error=row['e_z'], kind='heliocentric')
        events[name].add_quantity('ebv', row['E_B-V_'], source)
        events[name].add_quantity('claimedtype', 'Ia', source)
    result = Vizier.get_catalogs('J/ApJS/219/13/table2')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        name = 'LSQ' + row['LSQ']
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2015ApJS..219...13W')
        events[name].add_quantity('alias', name, source)
        add_photometry(
            events, name, time=str(jd_to_mjd(Decimal(row['JD']))), instrument='QUEST',
            observatory='La Silla', band=row['Filt'], telescope='ESO Schmidt',
            magnitude=row['mag'], e_magnitude=row['e_mag'], system='Swope', source=source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2012Natur.491..228C
    result = Vizier.get_catalogs('J/other/Nat/491.228/tablef1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    name = 'SN2213-1745'
    events, name = Events.add_event(tasks, args, events, name, log)
    source = events[name].add_source(bibcode='2012Natur.491..228C')
    events[name].add_quantity('alias', name, source)
    events[name].add_quantity('claimedtype', 'SLSN-R', source)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        for band in ['g', 'r', 'i']:
            bandtag = band + '_mag'
            if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                add_photometry(
                    events, name, time=row['MJD' + band + '_'], band=band + "'",
                    magnitude=row[bandtag], e_magnitude=row['e_' + bandtag], source=source)

    result = Vizier.get_catalogs('J/other/Nat/491.228/tablef2')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    name = 'SN1000+0216'
    events, name = Events.add_event(tasks, args, events, name, log)
    source = events[name].add_source(bibcode='2012Natur.491..228C')
    events[name].add_quantity('alias', name, source)
    events[name].add_quantity('claimedtype', 'SLSN-II?', source)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        for band in ['g', 'r', 'i']:
            bandtag = band + '_mag'
            if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                add_photometry(
                    events, name, time=row['MJD' + band + '_'], band=band + "'",
                    magnitude=row[bandtag], e_magnitude=row['e_' + bandtag], source=source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2011Natur.474..484Q
    result = Vizier.get_catalogs('J/other/Nat/474.484/tables1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        name = str(row['Name'])
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2011Natur.474..484Q')
        events[name].add_quantity('alias', name, source)
        add_photometry(
            events, name, time=row['MJD'], band=row['Filt'], telescope=row['Tel'],
            magnitude=row['mag'], e_magnitude=row['e_mag'], source=source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2011ApJ...736..159G
    result = Vizier.get_catalogs('J/ApJ/736/159/table1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    name = 'PTF10vdl'
    events, name = Events.add_event(tasks, args, events, name, log)
    source = events[name].add_source(bibcode='2011ApJ...736..159G')
    events[name].add_quantity('alias', name, source)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        add_photometry(
            events, name, time=str(jd_to_mjd(Decimal(row['JD']))), band=row['Filt'],
            telescope=row['Tel'], magnitude=row['mag'],
            e_magnitude=row['e_mag'] if is_number(row['e_mag']) else '',
            upperlimit=(not is_number(row['e_mag'])), source=source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2012ApJ...760L..33B
    result = Vizier.get_catalogs('J/ApJ/760/L33/table1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    name = 'PTF12gzk'
    events, name = Events.add_event(tasks, args, events, name, log)
    source = events[name].add_source(bibcode='2012ApJ...760L..33B')
    events[name].add_quantity('alias', name, source)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        # Fixing a typo in VizieR table
        if str(row['JD']) == '2455151.456':
            row['JD'] = '2456151.456'
        add_photometry(
            events, name, time=str(jd_to_mjd(Decimal(row['JD']))), band=row['Filt'],
            telescope=row['Inst'], magnitude=row['mag'], e_magnitude=row['e_mag'], source=source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2013ApJ...769...39S
    result = Vizier.get_catalogs('J/ApJ/769/39/table1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    name = 'PS1-12sk'
    events, name = Events.add_event(tasks, args, events, name, log)
    source = events[name].add_source(bibcode='2013ApJ...769...39S')
    events[name].add_quantity('alias', name, source)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        instrument = ''
        telescope = ''
        if row['Inst'] == 'RATCam':
            instrument = row['Inst']
        else:
            telescope = row['Inst']
        add_photometry(
            events, name, time=row['MJD'], band=row['Filt'], telescope=telescope,
            instrument=instrument, magnitude=row['mag'],
            e_magnitude=row['e_mag'] if not row['l_mag'] else '',
            upperlimit=(row['l_mag'] == '>'), source=source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2009MNRAS.394.2266P
    # Note: Instrument info available via links in VizieR, can't auto-parse just yet.
    name = 'SN2005cs'
    events, name = Events.add_event(tasks, args, events, name, log)
    source = events[name].add_source(bibcode='2009MNRAS.394.2266P')
    events[name].add_quantity('alias', name, source)
    result = Vizier.get_catalogs('J/MNRAS/394/2266/table2')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        for band in ['U', 'B', 'V', 'R', 'I']:
            bandtag = band + 'mag'
            if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                e_mag = (row['e_' + bandtag] if row['l_' + bandtag] != '>' else '')
                upl = row['l_' + bandtag] == '>'
                add_photometry(
                    events, name, time=str(jd_to_mjd(Decimal(row['JD']))), band=band,
                    magnitude=row[bandtag], e_magnitude=e_mag,
                    source=source, upperlimit=upl)
        if 'zmag' in row and is_number(row['zmag']) and not isnan(float(row['zmag'])):
            add_photometry(
                events, name, time=str(jd_to_mjd(Decimal(row['JD']))), band='z',
                magnitude=row['zmag'], e_magnitude=row['e_zmag'], source=source)

    result = Vizier.get_catalogs('J/MNRAS/394/2266/table3')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        for band in ['B', 'V', 'R']:
            bandtag = band + 'mag'
            if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                time = str(jd_to_mjd(Decimal(row['JD'])))
                e_mag = (row['e_' + bandtag] if row['l_' + bandtag] != '>' else '')
                add_photometry(
                    events, name, time=time, band=band, magnitude=row[bandtag],
                    e_magnitude=e_mag, source=source, upperlimit=(row['l_' + bandtag] == '>'))

    result = Vizier.get_catalogs('J/MNRAS/394/2266/table4')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        for band in ['J', 'H', 'K']:
            bandtag = band + 'mag'
            if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                add_photometry(
                    events, name, time=str(jd_to_mjd(Decimal(row['JD']))), band=band,
                    magnitude=row[bandtag],
                    e_magnitude=row['e_' + bandtag], source=source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2013AJ....145...99A
    result = Vizier.get_catalogs('J/AJ/145/99/table1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    name = 'SN2003ie'
    events, name = Events.add_event(tasks, args, events, name, log)
    source = events[name].add_source(bibcode='2013AJ....145...99A')
    events[name].add_quantity('alias', name, source)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        if 'Bmag' in row and is_number(row['Bmag']) and not isnan(float(row['Bmag'])):
            add_photometry(events, name, time=row['MJD'], band='B', magnitude=row['Bmag'],
                           e_magnitude=row['e_Bmag'] if not row['l_Bmag'] else '',
                           upperlimit=(row['l_Bmag'] == '>'), source=source)
        if 'Vmag' in row and is_number(row['Vmag']) and not isnan(float(row['Vmag'])):
            add_photometry(events, name, time=row['MJD'], band='V', magnitude=row['Vmag'],
                           e_magnitude=row['e_Vmag'] if is_number(row['e_Vmag']) else '',
                           upperlimit=(not is_number(row['e_Vmag'])), source=source)
        if 'Rmag' in row and is_number(row['Rmag']) and not isnan(float(row['Rmag'])):
            add_photometry(events, name, time=row['MJD'], band='R', magnitude=row['Rmag'],
                           e_magnitude=row['e_Rmag'] if not row['l_Rmag'] else '',
                           upperlimit=(row['l_Rmag'] == '>'), source=source)
        if 'Imag' in row and is_number(row['Imag']) and not isnan(float(row['Imag'])):
            add_photometry(events, name, time=row['MJD'], band='I', magnitude=row['Imag'],
                           e_magnitude=row['e_Imag'], source=source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2011ApJ...729..143C
    name = 'SN2008am'
    events, name = Events.add_event(tasks, args, events, name, log)
    source = events[name].add_source(bibcode='2011ApJ...729..143C')
    events[name].add_quantity('alias', name, source)

    result = Vizier.get_catalogs('J/ApJ/729/143/table1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        add_photometry(
            events, name, time=row['MJD'], band='ROTSE', telescope='ROTSE', magnitude=row['mag'],
            e_magnitude=row['e_mag'] if not row['l_mag'] else '', upperlimit=(row['l_mag'] == '<'), source=source)

    result = Vizier.get_catalogs('J/ApJ/729/143/table2')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        if 'Jmag' in row and is_number(row['Jmag']) and not isnan(float(row['Jmag'])):
            add_photometry(events, name, time=row['MJD'], telescope='PAIRITEL', band='J', magnitude=row['Jmag'],
                           e_magnitude=row['e_Jmag'], source=source)
        if 'Hmag' in row and is_number(row['Hmag']) and not isnan(float(row['Hmag'])):
            add_photometry(events, name, time=row['MJD'], telescope='PAIRITEL', band='H', magnitude=row['Hmag'],
                           e_magnitude=row['e_Hmag'], source=source)
        if 'Ksmag' in row and is_number(row['Ksmag']) and not isnan(float(row['Ksmag'])):
            add_photometry(events, name, time=row['MJD'], telescope='PAIRITEL', band='Ks', magnitude=row['Ksmag'],
                           e_magnitude=row['e_Ksmag'], source=source)

    result = Vizier.get_catalogs('J/ApJ/729/143/table4')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        add_photometry(
            events, name, time=row['MJD'], band=row['Filt'], telescope='P60',
            magnitude=row['mag'], e_magnitude=row['e_mag'], source=source)

    result = Vizier.get_catalogs('J/ApJ/729/143/table5')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        add_photometry(
            events, name, time=row['MJD'], band=row['Filt'], instrument='UVOT', telescope='Swift',
            magnitude=row['mag'], e_magnitude=row['e_mag'], source=source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2011ApJ...728...14P
    name = 'SN2009bb'
    events, name = Events.add_event(tasks, args, events, name, log)
    source = events[name].add_source(bibcode='2011ApJ...728...14P')
    events[name].add_quantity('alias', name, source)

    result = Vizier.get_catalogs('J/ApJ/728/14/table1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        if 'Bmag' in row and is_number(row['Bmag']) and not isnan(float(row['Bmag'])):
            add_photometry(
                events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'],
                band='B', magnitude=row['Bmag'], e_magnitude=row['e_Bmag'], source=source)
        if 'Vmag' in row and is_number(row['Vmag']) and not isnan(float(row['Vmag'])):
            add_photometry(
                events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'],
                band='V', magnitude=row['Vmag'], e_magnitude=row['e_Vmag'], source=source)
        if 'Rmag' in row and is_number(row['Rmag']) and not isnan(float(row['Rmag'])):
            add_photometry(
                events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'],
                band='R', magnitude=row['Rmag'], e_magnitude=row['e_Rmag'], source=source)
        if 'Imag' in row and is_number(row['Imag']) and not isnan(float(row['Imag'])):
            add_photometry(
                events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'],
                band='I', magnitude=row['Imag'], e_magnitude=row['e_Imag'], source=source)

    result = Vizier.get_catalogs('J/ApJ/728/14/table2')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        if 'u_mag' in row and is_number(row['u_mag']) and not isnan(float(row['u_mag'])):
            add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'], band='u', magnitude=row['u_mag'],
                           e_magnitude=row['e_u_mag'], source=source)
        if 'g_mag' in row and is_number(row['g_mag']) and not isnan(float(row['g_mag'])):
            add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'], band='g', magnitude=row['g_mag'],
                           e_magnitude=row['e_g_mag'], source=source)
        if 'r_mag' in row and is_number(row['r_mag']) and not isnan(float(row['r_mag'])):
            add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'], band='r', magnitude=row['r_mag'],
                           e_magnitude=row['e_r_mag'], source=source)
        if 'i_mag' in row and is_number(row['i_mag']) and not isnan(float(row['i_mag'])):
            add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'], band='i', magnitude=row['i_mag'],
                           e_magnitude=row['e_i_mag'], source=source)
        if 'z_mag' in row and is_number(row['z_mag']) and not isnan(float(row['z_mag'])):
            add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), telescope=row['Tel'], band='z', magnitude=row['z_mag'],
                           e_magnitude=row['e_z_mag'], source=source)

    result = Vizier.get_catalogs('J/ApJ/728/14/table3')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        if 'Ymag' in row and is_number(row['Ymag']) and not isnan(float(row['Ymag'])):
            add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), instrument=row['Inst'], band='Y', magnitude=row['Ymag'],
                           e_magnitude=row['e_Ymag'], source=source)
        if 'Jmag' in row and is_number(row['Jmag']) and not isnan(float(row['Jmag'])):
            add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), instrument=row['Inst'], band='J', magnitude=row['Jmag'],
                           e_magnitude=row['e_Jmag'], source=source)
        if 'Hmag' in row and is_number(row['Hmag']) and not isnan(float(row['Hmag'])):
            add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), instrument=row['Inst'], band='H', magnitude=row['Hmag'],
                           e_magnitude=row['e_Hmag'], source=source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2011PAZh...37..837T
    name = 'SN2009nr'
    events, name = Events.add_event(tasks, args, events, name, log)
    source = events[name].add_source(bibcode='2011PAZh...37..837T')
    events[name].add_quantity('alias', name, source)

    result = Vizier.get_catalogs('J/PAZh/37/837/table2')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        mjd = str(jd_to_mjd(Decimal(row['JD']) + 2455000))
        if 'Umag' in row and is_number(row['Umag']) and not isnan(float(row['Umag'])):
            add_photometry(events, name, time=mjd, telescope=row['Tel'], band='U', magnitude=row['Umag'],
                           e_magnitude=row['e_Umag'], source=source)
        if 'Bmag' in row and is_number(row['Bmag']) and not isnan(float(row['Bmag'])):
            add_photometry(events, name, time=mjd, telescope=row['Tel'], band='B', magnitude=row['Bmag'],
                           e_magnitude=row['e_Bmag'], source=source)
        if 'Vmag' in row and is_number(row['Vmag']) and not isnan(float(row['Vmag'])):
            add_photometry(events, name, time=mjd, telescope=row['Tel'], band='V', magnitude=row['Vmag'],
                           e_magnitude=row['e_Vmag'], source=source)
        if 'Rmag' in row and is_number(row['Rmag']) and not isnan(float(row['Rmag'])):
            add_photometry(events, name, time=mjd, telescope=row['Tel'], band='R', magnitude=row['Rmag'],
                           e_magnitude=row['e_Rmag'], source=source)
        if 'Imag' in row and is_number(row['Imag']) and not isnan(float(row['Imag'])):
            add_photometry(events, name, time=mjd, telescope=row['Tel'], band='I', magnitude=row['Imag'],
                           e_magnitude=row['e_Imag'], source=source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2013MNRAS.433.1871B
    name = 'SN2012aw'
    events, name = Events.add_event(tasks, args, events, name, log)
    source = events[name].add_source(bibcode='2013MNRAS.433.1871B')
    events[name].add_quantity('alias', name, source)

    result = Vizier.get_catalogs('J/MNRAS/433/1871/table3a')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        mjd = str(jd_to_mjd(Decimal(row['JD']) + 2456000))
        if 'Umag' in row and is_number(row['Umag']) and not isnan(float(row['Umag'])):
            add_photometry(events, name, time=mjd, telescope=row['Tel'], band='U', magnitude=row['Umag'],
                           e_magnitude=row['e_Umag'], source=source)
        if 'Bmag' in row and is_number(row['Bmag']) and not isnan(float(row['Bmag'])):
            add_photometry(events, name, time=mjd, telescope=row['Tel'], band='B', magnitude=row['Bmag'],
                           e_magnitude=row['e_Bmag'], source=source)
        if 'Vmag' in row and is_number(row['Vmag']) and not isnan(float(row['Vmag'])):
            add_photometry(events, name, time=mjd, telescope=row['Tel'], band='V', magnitude=row['Vmag'],
                           e_magnitude=row['e_Vmag'], source=source)
        if 'Rcmag' in row and is_number(row['Rcmag']) and not isnan(float(row['Rcmag'])):
            add_photometry(events, name, time=mjd, telescope=row['Tel'], band='Rc', magnitude=row['Rcmag'],
                           e_magnitude=row['e_Rcmag'], source=source)
        if 'Icmag' in row and is_number(row['Icmag']) and not isnan(float(row['Icmag'])):
            add_photometry(events, name, time=mjd, telescope=row['Tel'], band='Ic', magnitude=row['Icmag'],
                           e_magnitude=row['e_Icmag'], source=source)

    result = Vizier.get_catalogs('J/MNRAS/433/1871/table3b')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        mjd = str(jd_to_mjd(Decimal(row['JD']) + 2456000))
        if 'gmag' in row and is_number(row['gmag']) and not isnan(float(row['gmag'])):
            add_photometry(events, name, time=mjd, telescope=row['Tel'], band='g', magnitude=row['gmag'],
                           e_magnitude=row['e_gmag'], source=source)
        if 'rmag' in row and is_number(row['rmag']) and not isnan(float(row['rmag'])):
            add_photometry(events, name, time=mjd, telescope=row['Tel'], band='r', magnitude=row['rmag'],
                           e_magnitude=row['e_rmag'], source=source)
        if 'imag' in row and is_number(row['imag']) and not isnan(float(row['imag'])):
            add_photometry(events, name, time=mjd, telescope=row['Tel'], band='i', magnitude=row['imag'],
                           e_magnitude=row['e_imag'], source=source)
        if 'zmag' in row and is_number(row['zmag']) and not isnan(float(row['zmag'])):
            add_photometry(events, name, time=mjd, telescope=row['Tel'], band='z', magnitude=row['zmag'],
                           e_magnitude=row['e_zmag'], source=source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2014AJ....148....1Z
    name = 'SN2012fr'
    events, name = Events.add_event(tasks, args, events, name, log)
    source = events[name].add_source(bibcode='2014AJ....148....1Z')
    events[name].add_quantity('alias', name, source)

    result = Vizier.get_catalogs('J/AJ/148/1/table2')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        mjd = row['MJD']
        if 'Bmag' in row and is_number(row['Bmag']) and not isnan(float(row['Bmag'])):
            add_photometry(events, name, time=mjd, telescope='LJT', instrument='YFOSC', band='B', magnitude=row['Bmag'],
                           e_magnitude=row['e_Bmag'], source=source)
        if 'Vmag' in row and is_number(row['Vmag']) and not isnan(float(row['Vmag'])):
            add_photometry(events, name, time=mjd, telescope='LJT', instrument='YFOSC', band='V', magnitude=row['Vmag'],
                           e_magnitude=row['e_Vmag'], source=source)
        if 'Rmag' in row and is_number(row['Rmag']) and not isnan(float(row['Rmag'])):
            add_photometry(events, name, time=mjd, telescope='LJT', instrument='YFOSC', band='R', magnitude=row['Rmag'],
                           e_magnitude=row['e_Rmag'], source=source)
        if 'Imag' in row and is_number(row['Imag']) and not isnan(float(row['Imag'])):
            add_photometry(events, name, time=mjd, telescope='LJT', instrument='YFOSC', band='I', magnitude=row['Imag'],
                           e_magnitude=row['e_Imag'], source=source)

    result = Vizier.get_catalogs('J/AJ/148/1/table3')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        mjd = row['MJD']
        if 'Umag' in row and is_number(row['Umag']) and not isnan(float(row['Umag'])):
            add_photometry(events, name, time=mjd, telescope='Swift', instrument='UVOT', band='U', magnitude=row['Umag'],
                           e_magnitude=row['e_Umag'], source=source)
        if 'Bmag' in row and is_number(row['Bmag']) and not isnan(float(row['Bmag'])):
            add_photometry(events, name, time=mjd, telescope='Swift', instrument='UVOT', band='B', magnitude=row['Bmag'],
                           e_magnitude=row['e_Bmag'], source=source)
        if 'Vmag' in row and is_number(row['Vmag']) and not isnan(float(row['Vmag'])):
            add_photometry(events, name, time=mjd, telescope='Swift', instrument='UVOT', band='V', magnitude=row['Vmag'],
                           e_magnitude=row['e_Vmag'], source=source)
        if 'UVW1' in row and is_number(row['UVW1']) and not isnan(float(row['UVW1'])):
            add_photometry(events, name, time=mjd, telescope='Swift', instrument='UVOT', band='W1', magnitude=row['UVW1'],
                           e_magnitude=row['e_UVW1'], source=source)
        if 'UVW2' in row and is_number(row['UVW2']) and not isnan(float(row['UVW2'])):
            add_photometry(events, name, time=mjd, telescope='Swift', instrument='UVOT', band='W2', magnitude=row['UVW2'],
                           e_magnitude=row['e_UVW2'], source=source)
        if 'UVM2' in row and is_number(row['UVM2']) and not isnan(float(row['UVM2'])):
            add_photometry(events, name, time=mjd, telescope='Swift', instrument='UVOT', band='M2', magnitude=row['UVM2'],
                           e_magnitude=row['e_UVM2'], source=source)

    result = Vizier.get_catalogs('J/AJ/148/1/table5')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        mjd = row['MJD']
        if 'Bmag' in row and is_number(row['Bmag']) and not isnan(float(row['Bmag'])):
            add_photometry(events, name, time=mjd, telescope='LJT', band='B', magnitude=row['Bmag'],
                           e_magnitude=row['e_Bmag'], source=source)
        if 'Vmag' in row and is_number(row['Vmag']) and not isnan(float(row['Vmag'])):
            add_photometry(events, name, time=mjd, telescope='LJT', band='V', magnitude=row['Vmag'],
                           e_magnitude=row['e_Vmag'], source=source)
        if 'Rmag' in row and is_number(row['Rmag']) and not isnan(float(row['Rmag'])):
            add_photometry(events, name, time=mjd, telescope='LJT', band='R', magnitude=row['Rmag'],
                           e_magnitude=row['e_Rmag'], source=source)
        if 'Imag' in row and is_number(row['Imag']) and not isnan(float(row['Imag'])):
            add_photometry(events, name, time=mjd, telescope='LJT', band='I', magnitude=row['Imag'],
                           e_magnitude=row['e_Imag'], source=source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2015ApJ...805...74B
    name = 'SN2014J'
    events, name = Events.add_event(tasks, args, events, name, log)
    source = events[name].add_source(bibcode='2014AJ....148....1Z')
    events[name].add_quantity('alias', name, source)

    result = Vizier.get_catalogs('J/ApJ/805/74/table1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        mjd = row['MJD']
        if 'mag' in row and is_number(row['mag']) and not isnan(float(row['mag'])):
            add_photometry(events, name, time=mjd, telescope='Swift', instrument='UVOT', band=row['Filt'], magnitude=row['mag'],
                           e_magnitude=row['e_mag'], source=source)
        elif 'maglim' in row and is_number(row['maglim']) and not isnan(float(row['maglim'])):
            add_photometry(events, name, time=mjd, telescope='Swift', instrument='UVOT', band=row['Filt'], magnitude=row['maglim'],
                           upperlimit=True, source=source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2011ApJ...741...97D
    result = Vizier.get_catalogs('J/ApJ/741/97/table2')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        name = str(row['SN'])
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2011ApJ...741...97D')
        events[name].add_quantity('alias', name, source)
        add_photometry(events, name, time=str(jd_to_mjd(Decimal(row['JD']))), band=row['Filt'], magnitude=row['mag'],
                       e_magnitude=row['e_mag'] if is_number(row['e_mag']) else '', upperlimit=(not is_number(row['e_mag'])), source=source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2015MNRAS.448.1206M
    # Note: Photometry from two SN can also be added from this source.
    result = Vizier.get_catalogs('J/MNRAS/448/1206/table3')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        name = str(row['Name'])
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2015MNRAS.448.1206M')
        events[name].add_quantity('alias', name, source)
        events[name].add_quantity('discoverdate', '20' + name[4:6], source)
        events[name].add_quantity('ra', row['RAJ2000'], source, unit='floatdegrees')
        events[name].add_quantity('dec', row['DEJ2000'], source, unit='floatdegrees')
        events[name].add_quantity('redshift', row['zsp'], source, kind='spectroscopic')
        events[name].add_quantity('maxappmag', row['rP1mag'], source, error=row['e_rP1mag'])
        events[name].add_quantity('maxband', 'r', source)
        events[name].add_quantity('claimedtype', 'Ia', source)
    result = Vizier.get_catalogs('J/MNRAS/448/1206/table4')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        name = str(row['Name'])
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2015MNRAS.448.1206M')
        events[name].add_quantity('alias', name, source)
        events[name].add_quantity('discoverdate', '20' + name[4:6], source)
        events[name].add_quantity('ra', row['RAJ2000'], source, unit='floatdegrees')
        events[name].add_quantity('dec', row['DEJ2000'], source, unit='floatdegrees')
        events[name].add_quantity('redshift', row['zph'], source, error=row['e_zph'], kind='photometric')
        events[name].add_quantity('maxappmag', row['rP1mag'], source, error=row['e_rP1mag'])
        events[name].add_quantity('maxband', 'r', source)
        events[name].add_quantity('claimedtype', 'Ia?', source)
    result = Vizier.get_catalogs('J/MNRAS/448/1206/table5')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        name = str(row['Name'])
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2015MNRAS.448.1206M')
        events[name].add_quantity('alias', name, source)
        events[name].add_quantity('discoverdate', '20' + name[4:6], source)
        events[name].add_quantity('ra', row['RAJ2000'], source, unit='floatdegrees')
        events[name].add_quantity('dec', row['DEJ2000'], source, unit='floatdegrees')
        events[name].add_quantity('redshift', row['zsp'], source, kind='spectroscopic')
        events[name].add_quantity('maxappmag', row['rP1mag'], source, error=row['e_rP1mag'])
        events[name].add_quantity('maxband', 'r', source)
        events[name].add_quantity('claimedtype', row['Type'], source)
    result = Vizier.get_catalogs('J/MNRAS/448/1206/table6')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        name = str(row['Name'])
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2015MNRAS.448.1206M')
        events[name].add_quantity('alias', name, source)
        events[name].add_quantity('discoverdate', '20' + name[4:6], source)
        events[name].add_quantity('ra', row['RAJ2000'], source, unit='floatdegrees')
        events[name].add_quantity('dec', row['DEJ2000'], source, unit='floatdegrees')
        events[name].add_quantity('maxappmag', row['rP1mag'], source, error=row['e_rP1mag'])
        events[name].add_quantity('maxband', 'r', source)
        events[name].add_quantity('claimedtype', row['Type'], source)
    result = Vizier.get_catalogs('J/MNRAS/448/1206/tablea2')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        name = str(row['Name'])
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2015MNRAS.448.1206M')
        events[name].add_quantity('alias', name, source)
        events[name].add_quantity('discoverdate', '20' + name[4:6], source)
        events[name].add_quantity('ra', row['RAJ2000'], source, unit='floatdegrees')
        events[name].add_quantity('dec', row['DEJ2000'], source, unit='floatdegrees')
        events[name].add_quantity('maxappmag', row['rP1mag'], source, error=row['e_rP1mag'])
        events[name].add_quantity('maxband', 'r', source)
        events[name].add_quantity('claimedtype', row['Typesoft']+'?', source)
        events[name].add_quantity('claimedtype', row['Typepsnid']+'?', source)
    result = Vizier.get_catalogs('J/MNRAS/448/1206/tablea3')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        name = str(row['Name'])
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2015MNRAS.448.1206M')
        events[name].add_quantity('alias', name, source)
        events[name].add_quantity('discoverdate', '20' + name[4:6], source)
        events[name].add_quantity('ra', row['RAJ2000'], source, unit='floatdegrees')
        events[name].add_quantity('dec', row['DEJ2000'], source, unit='floatdegrees')
        events[name].add_quantity('maxappmag', row['rP1mag'], source, error=row['e_rP1mag'])
        events[name].add_quantity('maxband', 'r', source)
        events[name].add_quantity('claimedtype', 'Candidate', source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2012AJ....143..126B
    result = Vizier.get_catalogs('J/AJ/143/126/table4')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        if not row['Wcl'] or row['Wcl'] == 'N':
            continue
        row = convert_aq_output(row)
        name = str(row['SN']).replace(' ', '')
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2012AJ....143..126B')
        events[name].add_quantity('alias', name, source)
        events[name].add_quantity('claimedtype', 'Ia-' + row['Wcl'], source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    # 2015ApJS..220....9F
    for viztab in ['1', '2']:
        result = Vizier.get_catalogs('J/ApJS/220/9/table' + viztab)
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in pbar(table, current_task):
            row = convert_aq_output(row)
            events, name = Events.add_event(tasks, args, events, row['SN'], log)
            source = events[name].add_source(bibcode='2015ApJS..220....9F')
            events[name].add_quantity('alias', name, source)
            events[name].add_quantity('claimedtype', row['Type'], source)
            events[name].add_quantity('ra', row['RAJ2000'], source, unit='floatdegrees')
            events[name].add_quantity('dec', row['DEJ2000'], source, unit='floatdegrees')
            if '?' not in row['Host']:
                events[name].add_quantity('host', row['Host'].replace('_', ' '), source)
            kind = ''
            if 'Host' in row['n_z']:
                kind = 'host'
            elif 'Spectrum' in row['n_z']:
                kind = 'spectroscopic'
            events[name].add_quantity('redshift', row['z'], source, error=row['e_z'], kind=kind)

    result = Vizier.get_catalogs('J/ApJS/220/9/table8')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        events, name = Events.add_event(tasks, args, events, row['SN'], log)
        source = events[name].add_source(bibcode='2015ApJS..220....9F')
        events[name].add_quantity('alias', name, source)
        events[name].add_quantity('claimedtype', row['Type'], source)
        add_photometry(
            events, name, time=row['MJD'], band=row['Band'], magnitude=row['mag'],
            e_magnitude=row['e_mag'], telescope=row['Tel'], source=source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    result = Vizier.get_catalogs('J/ApJ/673/999/table1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        events, name = Events.add_event(tasks, args, events, 'SN'+row['SN'], log)
        source = events[name].add_source(bibcode='2008ApJ...673..999P')
        events[name].add_quantity('alias', name, source)
        events[name].add_quantity('ra', row['RAJ2000'], source, unit='floatdegrees')
        events[name].add_quantity('dec', row['DEJ2000'], source, unit='floatdegrees')
        events[name].add_quantity('redshift', row['z'], source, kind='host')
        events[name].add_quantity('hostra', row['RAGdeg'], source, unit='floatdegrees')
        events[name].add_quantity('hostdec', row['DEGdeg'], source, unit='floatdegrees')
        events[name].add_quantity('claimedtype', row['Type'].strip(':'), source)
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)

    return events


def do_lennarz(events, stubs, args, tasks, task_obj, log):
    """
    """
    current_task = task_obj.current_task(args)
    Vizier.ROW_LIMIT = -1
    result = Vizier.get_catalogs('J/A+A/538/A120/usc')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)

    bibcode = '2012A&A...538A.120L'
    for row in pbar(table, current_task):
        row = convert_aq_output(row)
        name = 'SN' + row['SN']
        events, name = Events.add_event(tasks, args, events, name, log)

        source = events[name].add_source(bibcode=bibcode)
        events[name].add_quantity('alias', name, source)

        if row['RAJ2000']:
            events[name].add_quantity('ra', row['RAJ2000'], source)
        if row['DEJ2000']:
            events[name].add_quantity('dec', row['DEJ2000'], source)
        if row['RAG']:
            events[name].add_quantity('hostra', row['RAG'], source)
        if row['DEG']:
            events[name].add_quantity('hostdec', row['DEG'], source)
        if row['Gal']:
            events[name].add_quantity('host', row['Gal'], source)
        if row['Type']:
            claimedtypes = row['Type'].split('|')
            for claimedtype in claimedtypes:
                events[name].add_quantity('claimedtype', claimedtype.strip(' -'), source)
        if row['z']:
            if name not in ['SN1985D', 'SN2004cq']:
                events[name].add_quantity('redshift', row['z'], source, kind='host')
        if row['Dist']:
            if row['e_Dist']:
                events[name].add_quantity('lumdist', row['Dist'], source, error=row['e_Dist'], kind='host')
            else:
                events[name].add_quantity('lumdist', row['Dist'], source, kind='host')

        if row['Ddate']:
            datestring = row['Ddate'].replace('-', '/')

            events[name].add_quantity('discoverdate', datestring, source)

            if 'photometry' not in events[name]:
                if 'Dmag' in row and is_number(row['Dmag']) and not isnan(float(row['Dmag'])):
                    datesplit = row['Ddate'].strip().split('-')
                    if len(datesplit) == 3:
                        datestr = row['Ddate'].strip()
                    elif len(datesplit) == 2:
                        datestr = row['Ddate'].strip() + '-01'
                    elif len(datesplit) == 1:
                        datestr = row['Ddate'].strip() + '-01-01'
                    mjd = str(astrotime(datestr).mjd)
                    add_photometry(events, name, time=mjd, band=row['Dband'], magnitude=row['Dmag'], source=source)
        if row['Mdate']:
            datestring = row['Mdate'].replace('-', '/')

            events[name].add_quantity('maxdate', datestring, source)

            if 'photometry' not in events[name]:
                if 'MMag' in row and is_number(row['MMag']) and not isnan(float(row['MMag'])):
                    datesplit = row['Mdate'].strip().split('-')
                    if len(datesplit) == 3:
                        datestr = row['Mdate'].strip()
                    elif len(datesplit) == 2:
                        datestr = row['Mdate'].strip() + '-01'
                    elif len(datesplit) == 1:
                        datestr = row['Mdate'].strip() + '-01-01'
                    mjd = str(astrotime(datestr).mjd)
                    add_photometry(events, name, time=mjd, band=row['Mband'], magnitude=row['Mmag'], source=source)

    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events


def do_snls_spectra(events, stubs, args, tasks, task_obj, log):
    """
    """
    from glob import glob
    from .. funcs import add_spectrum, get_preferred_name
    from .. constants import TRAVIS_QUERY_LIMIT
    from ... utils import pbar_strings

    current_task = task_obj.current_task(args)
    result = Vizier.get_catalogs('J/A+A/507/85/table1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    datedict = {}
    for row in table:
        datedict['SNLS-' + row['SN']] = str(astrotime(row['Date']).mjd)

    oldname = ''
    file_names = glob(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'SNLS/*'))
    for fi, fname in enumerate(pbar_strings(file_names, current_task=current_task)):
        filename = os.path.basename(fname)
        fileparts = filename.split('_')
        name = 'SNLS-' + fileparts[1]
        name = get_preferred_name(events, name)
        if oldname and name != oldname:
            events, stubs = Events.journal_events(tasks, args, events, stubs, log)
        oldname = name
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2009A&A...507...85B')
        events[name].add_quantity('alias', name, source)

        events[name].add_quantity('discoverdate', '20' + fileparts[1][:2], source)

        f = open(fname, 'r')
        data = csv.reader(f, delimiter=' ', skipinitialspace=True)
        specdata = []
        for r, row in enumerate(data):
            if row[0] == '@TELESCOPE':
                telescope = row[1].strip()
            elif row[0] == '@REDSHIFT':
                events[name].add_quantity('redshift', row[1].strip(), source)
            if r < 14:
                continue
            specdata.append(list(filter(None, [x.strip(' \t') for x in row])))
        specdata = [list(i) for i in zip(*specdata)]
        wavelengths = specdata[1]

        fluxes = [pretty_num(float(x)*1.e-16, sig=get_sig_digits(x)) for x in specdata[2]]
        # FIX: this isnt being used
        # errors = [pretty_num(float(x)*1.e-16, sig=get_sig_digits(x)) for x in specdata[3]]

        add_spectrum(
            name=name, waveunit='Angstrom', fluxunit='erg/s/cm^2/Angstrom', wavelengths=wavelengths,
            fluxes=fluxes, u_time='MJD' if name in datedict else '',
            time=datedict[name] if name in datedict else '', telescope=telescope, source=source,
            filename=filename)
        if args.travis and fi >= TRAVIS_QUERY_LIMIT:
            break
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events
