#!/usr/local/bin/python3.5

import csv
import glob
import os
import re
import urllib
import requests
import calendar
import sys
import json
import codecs
import resource
import argparse
import gzip
import io
import shutil
import warnings
import statistics
import warnings
from hashlib import md5
from tqdm import tqdm, trange
from html import unescape
from digits import *
from cdecimal import Decimal
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astroquery.irsa_dust import IrsaDust
from copy import deepcopy
from astropy import constants as const
from astropy import units as un
from astropy.io import fits
from astropy.time import Time as astrotime
from astropy.cosmology import Planck15 as cosmo, z_at_value
from astropy.coordinates import SkyCoord as coord
from collections import OrderedDict, Sequence
from math import log10, floor, sqrt, isnan, ceil
from bs4 import BeautifulSoup, Tag, NavigableString
from string import ascii_letters

parser = argparse.ArgumentParser(description='Generate a catalog JSON file and plot HTML files from SNE data.')
parser.add_argument('--update', '-u',       dest='update',      help='Only update catalog using live sources.',    default=False, action='store_true')
parser.add_argument('--verbose', '-v',      dest='verbose',     help='Print more messages to the screen.',         default=False, action='store_true')
parser.add_argument('--refresh', '-r',      dest='refresh',     help='Ignore most task caches.',                   default=False, action='store_true')
parser.add_argument('--full-refresh', '-f', dest='fullrefresh', help='Ignore all task caches.',                    default=False, action='store_true')
parser.add_argument('--travis', '-tr',      dest='travis',      help='Run import script in test mode for Travis.', default=False, action='store_true')
args = parser.parse_args()

tasks = OrderedDict([
    ("deleteoldevents", {"nicename":"Deleting old events",      "update": False}),
    ("internal",        {"nicename":"metadata and photometry",  "update": False}),
    ("radio",           {"nicename":"radio data",               "update": False}),
    ("xray",            {"nicename":"X-ray data",               "update": False}),
    ("simbad",          {"nicename":"SIMBAD",                   "update": False}),
    ("vizier",          {"nicename":"VizieR",                   "update": False}),
    ("donations",       {"nicename":"donations",                "update": False}),
    ("pessto-dr1",      {"nicename":"PESSTO DR1",               "update": False}),
    ("scp",             {"nicename":"SCP",                      "update": False}),
    ("ascii",           {"nicename":"ASCII",                    "update": False}),
    ("cccp",            {"nicename":"CCCP",                     "update": False, "archived": True}),
    ("suspect",         {"nicename":"SUSPECT",                  "update": False}),
    ("cfa",             {"nicename":"CfA archive photometry",   "update": False}),
    ("ucb",             {"nicename":"UCB photometry",           "update": False}),
    ("sdss",            {"nicename":"SDSS photometry",          "update": False}),
    ("csp",             {"nicename":"CSP photometry",           "update": False}),
    ("itep",            {"nicename":"ITEP",                     "update": False}),
    ("asiago",          {"nicename":"Asiago metadata",          "update": False}),
    ("tns",             {"nicename":"TNS metadata",             "update": True,  "archived": True}),
    ("rochester",       {"nicename":"Latest Supernovae",        "update": True,  "archived": False }),
    ("lennarz",         {"nicename":"Lennarz",                  "update": False}),
    ("gaia",            {"nicename":"GAIA",                     "update": True,  "archived": True}),
    ("ogle",            {"nicename":"%pre OGLE",                "update": True,  "archived": True}),
    ("snls",            {"nicename":"%pre SNLS",                "update": False}),
    ("psthreepi",       {"nicename":"%pre Pan-STARRS 3π",       "update": True,  "archived": True}),
    ("psmds",           {"nicename":"%pre Pan-STARRS MDS",      "update": False}),
    ("crts",            {"nicename":"%pre CRTS",                "update": True,  "archived": True}),
    ("snhunt",          {"nicename":"%pre SNhunt",              "update": True,  "archived": True}),
    ("nedd",            {"nicename":"%pre NED-D",               "update": False}),
    ("cpcs",            {"nicename":"%pre CPCS",                "update": True,  "archived": True}),
    ("ptf",             {"nicename":"%pre PTF",                 "update": False, "archived": False}),
    ("asiagospectra",   {"nicename":"%pre Asiago spectra",      "update": True }),
    ("wiserepspectra",  {"nicename":"%pre WISeREP spectra",     "update": False}),
    ("cfaspectra",      {"nicename":"%pre CfA archive spectra", "update": False}),
    ("snlsspectra",     {"nicename":"%pre SNLS spectra",        "update": False}),
    ("cspspectra",      {"nicename":"%pre CSP spectra",         "update": False}),
    ("ucbspectra",      {"nicename":"%pre UCB spectra",         "update": True, "archived": True}),
    ("suspectspectra",  {"nicename":"%pre SUSPECT spectra",     "update": False}),
    ("snfspectra",      {"nicename":"%pre SNH spectra",         "update": False}),
    ("superfitspectra", {"nicename":"%pre Superfit spectra",    "update": False}),
    ("writeevents",     {"nicename":"%pre writing events",      "update": True })
])

oscbibcode = '2016arXiv160501054G'
oscname = 'The Open Supernova Catalog'
oscurl = 'https://sne.space'
clight = const.c.cgs.value
km = (1.0 * un.km).cgs.value
planckh = const.h.cgs.value
keV = (1.0 * un.keV).cgs.value
travislimit = 10

eventnames = []
events = OrderedDict()
currenttask = ''

with open('rep-folders.txt', 'r') as f:
    repofolders = f.read().splitlines()

repoyears = [int(repofolders[x][-4:]) for x in range(len(repofolders))]
repoyears[0] -= 1

typereps = {
    'CC':       ['CCSN'],
    'Ia':       ['Ia-norm', 'Ia- norm'],
    'I P':      ['I pec', 'I-pec', 'I Pec', 'I-Pec'],
    'Ia P':     ['Ia pec', 'Ia-pec', 'Iapec', 'IaPec', 'Ia-Pec'],
    'Ib P':     ['Ib pec', 'Ib-pec'],
    'Ic P':     ['Ic pec', 'Ic-pec'],
    'Ia/c':     ['Ic/Ia', 'Iac'],
    'Ib/c':     ['Ibc'],
    'Ib/c P':   ['Ib/c-pec', 'Ibc pec', 'Ib/c pec'],
    'II P':     ['II pec', 'IIpec', 'II Pec', 'IIPec', 'IIP', 'IIp', 'II p',
                 'II-pec', 'II P pec', 'II-P', 'II-Pec', 'IIP-pec'],
    'II L':     ['IIL'],
    'IIn':      ['II n'],
    'IIn P':    ['IIn pec', 'IIn-pec'],
    'IIb P':    ['IIb-pec', 'IIb: pec', 'IIb pec'],
    'nIa':      ['nIa'],
    'Ia CSM':   ['Ia-CSM', 'Ia-csm', 'Ia-csm'],
    'SLSN-Ic':  ['SLSN Ic', 'SL-Ic'],
    'SLSN-I':   ['SLSN I', 'SL-I'],
    'SLSN-II':  ['SLSN II', 'SL-II'],
    'Ia-91bg':  ['Ia-pec (1991bg)', 'Ia-91bg-like', 'Ia-91bg like'],
    'Ia-91T':   ['Ia-pec 1991T', 'Ia-91T-like', 'Ia-91T like', 'Ia 91T-like'],
    'Ia-02cx':  ['Ia-02cx-like'],
    'Ib-Ca':    ['Ib - Ca-rich'],
    'II P-97D': ['IIP-pec 1997D'],
    'Ic BL':    ['Ic-broad'],
}

repbetterquantity = {
    'redshift',
    'ebv',
    'velocity',
    'lumdist',
    'discoverdate',
    'maxdate'
}

maxbands = [
    ['B', 'b', 'g'], # B-like bands first
    ['V'],           # if not, V-like bands
    ['R', 'r']       # if not, R-like bands
]

def tq(li, leave = True):
    return tqdm(list(li), desc = currenttask, leave = leave)

def tprint(string):
    try:
        tqdm.write(string)
    except:
        print(string)

def uniq_cdl(values):
    return ','.join(list(set(values)))

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

prefkinds = ['heliocentric', 'cmb', 'spectroscopic', 'photometric', 'host', 'cluster', '']
def frame_priority(attr):
    if 'kind' in attr:
        if attr['kind'] in prefkinds:
            return prefkinds.index(attr['kind'])
        else:
            return len(prefkinds)
    return len(prefkinds)

def alias_priority(name, attr):
    if name == attr:
        return 0
    return 1

def ct_priority(name, attr):
    aliases = attr['source'].split(',')
    max_source_year = -10000
    for alias in aliases:
        if alias == 'D':
            continue
        source = get_source_by_alias(name, alias)
        if 'bibcode' in source:
            source_year = get_source_year(source)
            if source_year > max_source_year:
                max_source_year = source_year
    return -max_source_year

def get_source_year(source):
    if 'bibcode' in source:
        if is_number(source['bibcode'][:4]):
            return int(source['bibcode'][:4])
        else:
            return -10000
    raise(ValueError('No bibcode available for source!'))

def name_clean(name):
    newname = name
    if newname.startswith('MASJ'):
        newname = newname.replace('MASJ', 'MASTER OT J', 1)
    if newname.startswith('PSNJ'):
        newname = newname.replace('PSNJ', 'PSN J', 1)
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
    if newname.startswith('SNF') and is_number(newname[3:]) and len(newname) >= 12:
        newname = 'SNF' + newname[3:11] + '-' + newname[11:]
    if newname.startswith('snf'):
        newname = newname.replace('snf', 'SNF', 1)
    if newname.startswith(('MASTER OT J', 'ROTSE3 J')):
        prefix = newname.split('J')[0]
        coords = newname.split('J')[-1]
        decsign = '+' if '+' in coords else '-'
        coordsplit = coords.replace('+','-').split('-')
        if '.' not in coordsplit[0] and len(coordsplit[0]) > 6 and '.' not in coordsplit[1] and len(coordsplit[1]) > 6:
            newname = (prefix + 'J' + coordsplit[0][:6] + '.' + coordsplit[0][6:] +
                decsign + coordsplit[1][:6] + '.' + coordsplit[1][6:])
    if newname.startswith('SN ') and is_number(newname[3:7]) and len(newname) > 7:
        newname = newname.replace('SN ', 'SN', 1)
    if newname.startswith('SN') and is_number(newname[2:6]) and len(newname) == 7 and newname[6].islower():
        newname = 'SN' + newname[2:6] + newname[6].upper()
    elif (newname.startswith('SN') and is_number(newname[2:6]) and
        (len(newname) == 8 or len(newname) == 9) and newname[6:].isupper()):
        newname = 'SN' + newname[2:6] + newname[6:].lower()
    return newname

def get_aliases(name, includename = True):
    if 'alias' in events[name]:
        aliases = [x['value'] for x in events[name]['alias']]
        if includename and name not in aliases:
            return [name] + aliases
        return aliases
    if includename:
        return [name]
    return []

def add_event(name, load = True, delete = True, source = '', loadifempty = True):
    if loadifempty and args.update and not len(events):
        load_stubs()

    newname = name_clean(name)
    if newname not in events or 'stub' in events[newname]:
        match = ''
        if newname not in events:
            for event in events:
                aliases = get_aliases(event)
                if (len(aliases) > 1 and newname in aliases and
                    ('distinctfrom' not in events[event] or newname not in events[event]['distinctfrom'])):
                    match = event
                    break
            if match:
                newname = match

        if load:
            loadedname = load_event_from_file(name = newname, delete = delete)
            if loadedname:
                if 'stub' in events[loadedname]:
                    raise(ValueError('Failed to find event file for stubbed event'))
                return loadedname

        if match:
            return match

        events[newname] = OrderedDict()
        events[newname]['name'] = newname
        if source:
            add_quantity(newname, 'alias', newname, source)
        if args.verbose and 'stub' not in events[newname]:
            tprint('Added new event ' + newname)
        return newname
    else:
        return newname

def event_exists(name):
    if name in events:
        return True
    for ev in events:
        if name in get_aliases(ev):
            return True
    return False

def get_preferred_name(name):
    if name not in events:
        matches = []
        for event in events:
            aliases = get_aliases(event)
            if len(aliases) > 1 and name in aliases:
                return event
        return name
    else:
        return name

def get_event_filename(name):
    return(name.replace('/', '_'))

def snname(string):
    newstring = string.replace(' ', '').upper()
    if (newstring[:2] == "SN"):
        head = newstring[:6]
        tail = newstring[6:]
        if len(tail) >= 2 and tail[1] != '?':
            tail = tail.lower()
        newstring = head + tail

    return newstring

def add_source(name, reference = '', url = '', bibcode = '', secondary = ''):
    nsources = len(events[name]['sources']) if 'sources' in events[name] else 0
    if not reference:
        if not bibcode:
            raise(ValueError('Bibcode must be specified if name is not.'))

        if bibcode and len(bibcode) != 19:
            raise(ValueError('Bibcode "' + bibcode + '" must be exactly 19 characters long'))

        reference = bibcode

    reference = reference.replace('ATEL', 'ATel').replace('Atel', 'ATel').replace('ATel #', 'ATel ').replace('ATel#', 'ATel').replace('ATel', 'ATel ')
    reference = reference.replace('CBET', 'CBET ')
    reference = reference.replace('IAUC', 'IAUC ')
    reference = ' '.join(reference.split())

    if reference.startswith('ATel ') and not bibcode:
        atelnum = reference.split()[-1]
        if is_number(atelnum) and atelnum in atelsdict:
            bibcode = atelsdict[atelnum]

    if reference.startswith('CBET ') and not bibcode:
        cbetnum = reference.split()[-1]
        if is_number(cbetnum) and cbetnum in cbetsdict:
            bibcode = cbetsdict[cbetnum]

    if reference.startswith('IAUC ') and not bibcode:
        cbetnum = reference.split()[-1]
        if is_number(cbetnum) and cbetnum in iaucsdict:
            bibcode = iaucsdict[cbetnum]

    if 'sources' not in events[name] or (reference not in [x['name'] for x in events[name]['sources']] and
        (not bibcode or bibcode not in [x['bibcode'] if 'bibcode' in x else '' for x in events[name]['sources']])):
        source = str(nsources + 1)
        newsource = OrderedDict()
        newsource['name'] = reference
        if url:
            newsource['url'] = url
        if bibcode:
            newsource['bibcode'] = bibcode
        newsource['alias'] =  source
        if secondary:
            newsource['secondary'] = True
        events[name].setdefault('sources',[]).append(newsource)
    else:
        if reference in [x['name'] for x in events[name]['sources']]:
            source = [x['alias'] for x in events[name]['sources']][
                [x['name'] for x in events[name]['sources']].index(reference)]
        elif bibcode and bibcode in [x['bibcode'] if 'bibcode' in x else '' for x in events[name]['sources']]:
            source = [x['alias'] for x in events[name]['sources']][
                [x['bibcode'] if 'bibcode' in x else '' for x in events[name]['sources']].index(bibcode)]
        else:
            raise(ValueError("Couldn't find source that should exist!"))
    return source

def get_source_by_alias(name, alias):
    for source in events[name]['sources']:
        if source['alias'] == alias:
            return source
    raise(ValueError('Source alias not found!'))

def add_photometry(name, time = "", u_time = "MJD", e_time = "", telescope = "", instrument = "", band = "",
                   magnitude = "", e_magnitude = "", source = "", upperlimit = False, system = "",
                   observatory = "", observer = "", host = False, includeshost = False, survey = "",
                   flux = "", fluxdensity = "", e_flux = "", e_fluxdensity = "", u_flux = "", u_fluxdensity = "", frequency = "",
                   u_frequency = "", counts = "", e_counts = "", nhmw = "", photonindex = "", unabsorbedflux = "",
                   e_unabsorbedflux = "", energy = "", u_energy = "", e_lower_magnitude = "", e_upper_magnitude = ""):
    if (not time and not host) or (not magnitude and not flux and not fluxdensity and not counts and not unabsorbedflux):
        warnings.warn('Time or brightness not specified when adding photometry, not adding.\n')
        tprint('Name : "' + name + '", Time: "' + time + '", Band: "' + band + '", AB magnitude: "' + magnitude + '"')
        return

    if (not host and not is_number(time)) or (not is_number(magnitude) and not is_number(flux) and not is_number(fluxdensity) and not is_number(counts)):
        warnings.warn('Time or brightness not numerical, not adding.\n')
        tprint('Name : "' + name + '", Time: "' + time + '", Band: "' + band + '", AB magnitude: "' + magnitude + '"')
        return

    if (not upperlimit and ((e_magnitude and not is_number(e_magnitude)) or (flux and not is_number(e_flux)) or
        (fluxdensity and not is_number(e_fluxdensity)) or (counts and not is_number(e_counts)))):
        warnings.warn('Brightness error not numerical, not adding.\n')
        tprint('Name : "' + name + '", Time: "' + time + '", Band: "' + band + '", AB error: "' + e_magnitude + '"')
        return

    if e_time and not is_number(e_time):
        warnings.warn('Time error not numerical, not adding.\n')
        tprint('Name : "' + name + '", Time: "' + time + '", Time error: "' + e_time + '"')
        return

    if (flux or fluxdensity) and ((not u_flux and not u_fluxdensity) or (not frequency and not band and not energy)):
        warnings.warn('Unit and band/frequency must be set when adding photometry by flux or flux density, not adding.')
        tprint('Name : "' + name + '", Time: "' + time)
        return

    if not source:
        ValueError('Photometry must have source before being added!')

    if is_erroneous(name, 'photometry', source):
        return

    # Look for duplicate data and don't add if duplicate
    if 'photometry' in events[name]:
        for photo in events[name]['photometry']:
            if ((('host' not in photo and not host) or ('host' in photo and host)) and
                (('u_time' not in photo and not u_time) or (photo['u_time'] == u_time)) and
                (('time' not in photo and not time) or
                 (isinstance(photo['time'], str) and isinstance(time, str) and Decimal(photo['time']) == Decimal(time)) or
                 (isinstance(photo['time'], list) and isinstance(time, list) and photo['time'] == time)) and
                (('magnitude' not in photo and not magnitude) or
                 ('magnitude' in photo and Decimal(photo['magnitude']) == Decimal(magnitude)) or
                 ('magnitude' in photo and not magnitude)) and
                (('flux' not in photo and not flux) or
                 ('flux' in photo and Decimal(photo['flux']) == Decimal(flux)) or
                 ('flux' in photo and not flux)) and
                (('unabsorbedflux' not in photo and not unabsorbedflux) or
                 ('unabsorbedflux' in photo and Decimal(photo['unabsorbedflux']) == Decimal(unabsorbedflux)) or
                 ('unabsorbedflux' in photo and not unabsorbedflux)) and
                (('fluxdensity' not in photo and not fluxdensity) or
                 ('fluxdensity' in photo and Decimal(photo['fluxdensity']) == Decimal(fluxdensity)) or
                 ('fluxdensity' in photo and not fluxdensity)) and
                (('counts' not in photo and not counts) or
                 ('counts' in photo and Decimal(photo['counts']) == Decimal(counts)) or
                 ('counts' in photo and not counts)) and
                (('energy' not in photo and not energy) or
                 ('energy' in photo and Decimal(photo['energy']) == Decimal(energy)) or
                 ('energy' in photo and not energy)) and
                (('frequency' not in photo and not frequency) or
                 ('frequency' in photo and Decimal(photo['frequency']) == Decimal(frequency)) or
                 ('frequency' in photo and not frequency)) and
                (('band' not in photo and not band) or
                 ('band' in photo and photo['band'] == band) or
                 ('band' in photo and not band)) and
                (('photonindex' not in photo and not photonindex) or
                 ('photonindex' in photo and photo['photonindex'] == photonindex) or
                 ('photonindex' in photo and not photonindex)) and
                (('e_magnitude' not in photo and not e_magnitude) or
                 ('e_magnitude' in photo and e_magnitude and Decimal(photo['e_magnitude']) == Decimal(e_magnitude)) or
                 ('e_magnitude' in photo and not e_magnitude)) and
                (('e_flux' not in photo and not e_flux) or
                 ('e_flux' in photo and e_flux and Decimal(photo['e_flux']) == Decimal(e_flux)) or
                 ('e_flux' in photo and not e_flux)) and
                (('e_unabsorbedflux' not in photo and not e_unabsorbedflux) or
                 ('e_unabsorbedflux' in photo and e_unabsorbedflux and Decimal(photo['e_unabsorbedflux']) == Decimal(e_unabsorbedflux)) or
                 ('e_unabsorbedflux' in photo and not e_unabsorbedflux)) and
                (('e_fluxdensity' not in photo and not e_fluxdensity) or
                 ('e_fluxdensity' in photo and e_fluxdensity and Decimal(photo['e_fluxdensity']) == Decimal(e_fluxdensity)) or
                 ('e_fluxdensity' in photo and not e_fluxdensity)) and
                (('e_counts' not in photo and not e_counts) or
                 ('e_counts' in photo and e_counts and Decimal(photo['e_counts']) == Decimal(e_counts)) or
                 ('e_counts' in photo and not e_counts)) and
                (('system' not in photo and not system) or
                 ('system' in photo and photo['system'] == system) or
                 ('system' in photo and not system))):
                return

    photoentry = OrderedDict()
    if time:
        photoentry['time'] = time if isinstance(time, list) or isinstance(time, str) else str(time)
    if e_time:
        photoentry['e_time'] = str(e_time)
    if u_time:
        photoentry['u_time'] = u_time
    if band:
        photoentry['band'] = band
    if system:
        photoentry['system'] = system
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
    if telescope:
        photoentry['telescope'] = telescope
    if instrument:
        photoentry['instrument'] = instrument
    if nhmw:
        photoentry['nhmw'] = nhmw
    if source:
        photoentry['source'] = source
    events[name].setdefault('photometry',[]).append(photoentry)

def trim_str_arr(arr, length = 10):
    return [str(round_sig(float(x), length)) if (len(x) > length and len(str(round_sig(float(x), length))) < len(x)) else x for x in arr]

def add_spectrum(name, waveunit, fluxunit, wavelengths = "", fluxes = "", u_time = "", time = "", instrument = "",
    deredshifted = "", dereddened = "", errorunit = "", errors = "", source = "", snr = "", telescope = "",
    observer = "", reducer = "", filename = "", observatory = "", data = ""):

    if is_erroneous(name, 'spectra', source):
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
    events[name].setdefault('spectra',[]).append(spectrumentry)

def is_erroneous(name, field, sources):
    if 'error' in events[name]:
        for alias in sources.split(','):
            source = get_source_by_alias(name, alias)
            if ('bibcode' in source and source['bibcode'] in
                [x['extra'] for x in events[name]['error'] if x['kind'] == 'bibcode' and x['value'] == field]):
                    return True
            if ('name' in source and source['name'] in
                [x['extra'] for x in events[name]['error'] if x['kind'] == 'name' and x['value'] == field]):
                    return True
    return False

def add_quantity(name, quantity, value, sources, forcereplacebetter = False,
    lowerlimit = '', upperlimit = '', error = '', unit = '', kind = '', extra = ''):
    if not quantity:
        raise(ValueError('Quantity must be specified for add_quantity.'))
    if not sources:
        raise(ValueError('Source must be specified for quantity before it is added.'))
    if not isinstance(value, str) and (not isinstance(value, list) or not isinstance(value[0], str)):
        raise(ValueError('Quantity must be a string or an array of strings.'))

    if is_erroneous(name, quantity, sources):
        return

    svalue = value.strip()
    serror = error.strip()
    sunit = ''

    if not svalue or svalue == '--' or svalue == '-':
        return
    if serror and (not is_number(serror) or float(serror) < 0.):
        raise(ValueError('Quanta error value must be a number and positive.'))

    #Set default units
    if not unit and quantity == 'velocity':
        unit = 'km/s' 
    if not unit and quantity == 'ra':
        unit = 'hours'
    if not unit and quantity == 'dec':
        unit = 'degrees'
    if not unit and quantity in ['lumdist', 'comovingdist']:
        unit = 'Mpc'

    #Handle certain quantity
    if quantity == 'alias':
        newalias = name_clean(svalue)
        if 'distinctfrom' in events[name]:
            if newalias in [x['value'] for x in events[name]['distinctfrom']]:
                return
    if quantity in ['velocity', 'redshift', 'ebv', 'lumdist', 'comovingdist']:
        if not is_number(svalue):
            return
    if quantity == 'host':
        if is_number(svalue):
            return
        svalue = svalue.strip("()").replace('  ', ' ')
        svalue = svalue.replace("APMUKS(BJ)", "APMUKS(BJ) ")
        svalue = svalue.replace("ARP", "ARP ")
        svalue = svalue.replace("CGCG", "CGCG ")
        svalue = svalue.replace("HOLM", "HOLM ")
        svalue = svalue.replace("IC", "IC ")
        svalue = svalue.replace("Intergal.", "Intergalactic")
        svalue = svalue.replace("MCG+", "MCG +")
        svalue = svalue.replace("MCG-", "MCG -")
        svalue = svalue.replace("M+", "MCG +")
        svalue = svalue.replace("M-", "MCG -")
        svalue = svalue.replace("MGC ", "MCG ")
        svalue = svalue.replace("Mrk", "MRK")
        svalue = svalue.replace("MRK", "MRK ")
        svalue = svalue.replace("NGC", "NGC ")
        svalue = svalue.replace("PGC", "PGC ")
        svalue = svalue.replace("SDSS", "SDSS ")
        svalue = svalue.replace("UGC", "UGC ")
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
            svalue = svalue[:3] + ':' + svalue[3:5] + ((':' + zpad(svalue[5:])) if len(svalue) > 5 else '')
        else:
            svalue = svalue.replace(' ', ':')
            if 'dec' in quantity:
                valuesplit = svalue.split(':')
                svalue = (('+' if float(valuesplit[0]) > 0.0 else '-') + valuesplit[0].strip('+-').zfill(2) +
                    (':' + valuesplit[1].zfill(2) if len(valuesplit) > 1 else '') +
                    (':' + zpad(valuesplit[2]) if len(valuesplit) > 2 else ''))

        if 'ra' in quantity:
            sunit = 'hours'
        elif 'dec' in quantity:
            sunit = 'degrees'
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
                if 'kind' in ct and kind and ct['kind'] != kind:
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
    if kind:
        quantaentry['kind'] = kind
    if sunit:
        quantaentry['unit'] = sunit
    if lowerlimit:
        quantaentry['lowerlimit'] = lowerlimit
    if upperlimit:
        quantaentry['upperlimit'] = upperlimit
    if extra:
        quantaentry['extra'] = extra
    if (forcereplacebetter or quantity in repbetterquantity) and quantity in events[name]:
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
                    if max(2,get_sig_digits(ctsplit[-1].lstrip('0'))) < max(2,get_sig_digits(svsplit[-1].lstrip('0'))):
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
        events[name].setdefault(quantity,[]).append(quantaentry)

def load_cached_url(url, filepath, timeout = 120):
    filemd5 = ''
    filetxt = ''
    if not args.refresh and os.path.isfile(filepath):
        with open(filepath, 'r') as f:
            filetxt = f.read()
            if args.update:
                filemd5 = md5(filetxt.encode('utf-8')).hexdigest()

    try:
        session = requests.Session()
        response = session.get(url, timeout = timeout)
        txt = response.text
        newmd5 = md5(txt.encode('utf-8')).hexdigest()
        tprint('"' + filemd5 + '", "' + newmd5 + '"')
        if args.update and newmd5 == filemd5:
            tprint('Skipping file in "' + currenttask + '," local and remote copies identical.')
            return False
    except:
        return filetxt
    else:
        with open(filepath, 'w') as f:
            f.write(txt)
    return txt

def make_date_string(year, month = '', day = ''):
    if not year:
        raise ValueError('At least the year must be specified when constructing date string')
    datestring = str(year)
    if month:
        datestring = datestring + '/' + str(month).zfill(2)
    if day:
        datestring = datestring + '/' + str(day).zfill(2)

    return datestring

def get_max_light(name):
    if 'photometry' not in events[name]:
        return (None, None, None, None)

    eventphoto = [(x['u_time'], x['time'], Decimal(x['magnitude']), x['band'] if 'band' in x else '', x['source']) for x in events[name]['photometry'] if
                  ('magnitude' in x and 'time' in x and 'u_time' in x and 'upperlimit' not in x)]
    if not eventphoto:
        return (None, None, None, None)

    mlmag = None
    for mb in maxbands:
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

def get_first_light(name):
    if 'photometry' not in events[name]:
        return (None, None)

    eventphoto = [(Decimal(x['time']) if isinstance(x['time'], str) else Decimal(min(float(y) for y in x['time'])),
        x['source']) for x in events[name]['photometry'] if 'upperlimit' not in x
        and 'time' in x and 'u_time' in x and x['u_time'] == 'MJD']
    if not eventphoto:
        return (None, None)
    flmjd = min([x[0] for x in eventphoto])
    flindex = [x[0] for x in eventphoto].index(flmjd)
    flmjd = float(flmjd)
    flsource = eventphoto[flindex][1]
    return (astrotime(flmjd, format='mjd').datetime, flsource)

def set_first_max_light(name):
    if 'maxappmag' not in events[name]:
        (mldt, mlmag, mlband, mlsource) = get_max_light(name)
        if mldt:
            source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
            add_quantity(name, 'maxdate', make_date_string(mldt.year, mldt.month, mldt.day), uniq_cdl([source,mlsource]))
        if mlmag:
            source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
            add_quantity(name, 'maxappmag', pretty_num(mlmag), uniq_cdl([source,mlsource]))
        if mlband:
            source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
            add_quantity(name, 'maxband', mlband, uniq_cdl([source,mlsource]))

    if 'discoverdate' not in events[name] or max([len(x['value'].split('/')) for x in events[name]['discoverdate']]) < 3:
        (fldt, flsource) = get_first_light(name)
        if fldt:
            source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
            add_quantity(name, 'discoverdate', make_date_string(fldt.year, fldt.month, fldt.day), uniq_cdl([source,flsource]))

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
            source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
            add_quantity(name, 'discoverdate', make_date_string(fldt.year, fldt.month, fldt.day), 'D,' + minspecsource)

def get_best_redshift(name):
    bestsig = -1
    bestkind = 10
    for z in events[name]['redshift']:
        kind = prefkinds.index(z['kind'] if 'kind' in z else '')
        sig = get_sig_digits(z['value'])
        if sig > bestsig and kind <= bestkind:
            bestz = z['value']
            bestkind = kind
            bestsig = sig

    return (bestz, bestkind, bestsig)

def jd_to_mjd(jd):
    return jd - Decimal(2400000.5)

def utf8(x):
    return str(x, 'utf-8')

def convert_aq_output(row):
    return OrderedDict([(x, str(row[x]) if is_number(row[x]) else row[x]) for x in row.colnames])

def set_preferred_names():
    for name in list(sorted(list(events.keys()))):
        if name not in events:
            continue
        newname = ''
        aliases = get_aliases(name)
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
            discoverer = ','.join(x['value'].upper() for x in events[name]['discoverer'])
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
            if load_event_from_file(newname):
                continue
            if load_event_from_file(name, delete = True):
                tprint('Changing event name (' + name + ') to preferred name (' + newname + ').')
                events[newname] = events[name]
                events[newname]['name'] = newname
                del(events[name])
                journal_events()

# Merge and remove duplicate events
def merge_duplicates():
    keys = list(sorted(list(events.keys())))
    for n1, name1 in enumerate(tq(keys[:])):
        if name1 not in events:
            continue
        allnames1 = get_aliases(name1) + (['AT' + name1[2:]] if (name1.startswith('SN') and is_number(name1[2:6])) else [])
        for name2 in keys[n1+1:]:
            if name2 not in events or name1 == name2:
                continue
            allnames2 = get_aliases(name2) + (['AT' + name2[2:]] if (name2.startswith('SN') and is_number(name2[2:6])) else [])
            if any(i in allnames1 for i in allnames2):
                tprint('Found single event with multiple entries (' + name1 + ' and ' + name2 + '), merging.')
                load1 = load_event_from_file(name1, delete = True)
                load2 = load_event_from_file(name2, delete = True)
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
                        copy_to_event(name2, name1)
                        keys.append(name1)
                        del(events[name2])
                    else:
                        copy_to_event(name1, name2)
                        keys.append(name2)
                        del(events[name1])
                else:
                    print ('Duplicate already deleted')
                journal_events()

def derive_and_sanitize():
    biberrordict = {
        "2012Sci..337..942D":"2012Sci...337..942D",
        "2012MNRAS.420.1135":"2012MNRAS.420.1135S",
        "2014MNRAS.438,368":"2014MNRAS.438..368T",
        "2006ApJ...636...400Q":"2006ApJ...636..400Q",
        "0609268":"2007AJ....133...58K",
        "2004MNRAS.tmp..131P":"2004MNRAS.352..457P",
        "2013MNRAS.tmp.1499F":"2013MNRAS.433.1312F",
        "1991MNRAS.247P.410B":"1991A&A...247..410B",
        "2011Sci.333..856S":"2011Sci...333..856S"
    }

    # Calculate some columns based on imported data, sanitize some fields
    for name in events:
        aliases = get_aliases(name, includename = False)
        if name not in aliases:
            if 'sources' in events[name]:
                add_quantity(name, 'alias', name, '1')
            else:
                source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
                add_quantity(name, 'alias', name, source)
        events[name]['alias'] = list(sorted(events[name]['alias'], key=lambda key: alias_priority(name, key)))
        aliases = get_aliases(name)

        set_first_max_light(name)

        if 'claimedtype' in events[name]:
            events[name]['claimedtype'] = list(sorted(events[name]['claimedtype'], key=lambda key: ct_priority(name, key)))
        if 'discoverdate' not in events[name]:
            prefixes = ['MLS', 'SSS', 'CSS']
            for alias in aliases:
                for prefix in prefixes:
                    if alias.startswith(prefix) and is_number(alias.replace(prefix, '')[:2]):
                        discoverdate = '/'.join(['20' + alias.replace(prefix, '')[:2],
                            alias.replace(prefix, '')[2:4], alias.replace(prefix, '')[4:6]])
                        print ('Added discoverdate from name: ' + discoverdate)
                        source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
                        add_quantity(name, 'discoverdate', discoverdate, source)
                        break
                if 'discoverdate' in events[name]:
                    break
        if 'discoverdate' not in events[name]:
            prefixes = ['ASASSN-', 'PS1-', 'PS1', 'PS', 'iPTF', 'PTF', 'SCP-']
            for alias in aliases:
                for prefix in prefixes:
                    if alias.startswith(prefix) and is_number(alias.replace(prefix, '')[:2]):
                        discoverdate = '20' + alias.replace(prefix, '')[:2]
                        print ('Added discoverdate from name: ' + discoverdate)
                        source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
                        add_quantity(name, 'discoverdate', discoverdate, source)
                        break
                if 'discoverdate' in events[name]:
                    break
        if 'discoverdate' not in events[name]:
            prefixes = ['SNF']
            for alias in aliases:
                for prefix in prefixes:
                    if alias.startswith(prefix) and is_number(alias.replace(prefix, '')[:4]):
                        discoverdate = '/'.join([alias.replace(prefix, '')[:4],
                            alias.replace(prefix, '')[4:6], alias.replace(prefix, '')[6:8]])
                        print ('Added discoverdate from name: ' + discoverdate)
                        source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
                        add_quantity(name, 'discoverdate', discoverdate, source)
                        break
                if 'discoverdate' in events[name]:
                    break
        if 'discoverdate' not in events[name]:
            prefixes = ['AT', 'SN']
            for alias in aliases:
                for prefix in prefixes:
                    if alias.startswith(prefix) and is_number(alias.replace(prefix, '')[:4]):
                        discoverdate = alias.replace(prefix, '')[:4]
                        print ('Added discoverdate from name: ' + discoverdate)
                        source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
                        add_quantity(name, 'discoverdate', discoverdate, source)
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
                        noprefix = noprefix.replace('+','|').replace('-','|')
                        nops = noprefix.split('|')
                        if len(nops) < 2:
                            continue
                        rastr = nops[0]
                        decstr = nops[1]
                        ra = ':'.join([rastr[:2], rastr[2:4], rastr[4:6]]) + ('.' + rastr[6:] if len(rastr) > 6 else '') 
                        dec = decsign + ':'.join([decstr[:2], decstr[2:4], decstr[4:6]]) + ('.' + decstr[6:] if len(decstr) > 6 else '')
                        print ('Added ra/dec from name: ' + ra + ' ' + dec)
                        source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
                        add_quantity(name, 'ra', ra, source)
                        add_quantity(name, 'dec', ra, source)
                        break
                if 'ra' in events[name]:
                    break
        if ('ra' in events[name] and 'dec' in events[name] and 
            (not 'host' in events[name] or not any([x['value'] == 'Milky Way' for x in events[name]['host']]))):
            if name not in extinctionsdict:
                result = IrsaDust.get_query_table(events[name]['ra'][0]['value'] + " " + events[name]['dec'][0]['value'], section = 'ebv')
                ebv = result['ext SandF mean'][0]
                ebverr = result['ext SandF std'][0]
                extinctionsdict[name] = [ebv, ebverr]
            source = add_source(name, bibcode = '2011ApJ...737..103S')
            add_quantity(name, 'ebv', str(extinctionsdict[name][0]), source, error = str(extinctionsdict[name][1]))
        if 'claimedtype' in events[name]:
            events[name]['claimedtype'][:] = [ct for ct in events[name]['claimedtype'] if (ct['value'] != '?' and ct['value'] != '-')]
        if 'claimedtype' not in events[name] and name.startswith('AT'):
            source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
            add_quantity(name, 'claimedtype', 'Candidate', source)
        if 'redshift' not in events[name] and 'velocity' in events[name]:
            # Find the "best" velocity to use for this
            bestsig = 0
            for hv in events[name]['velocity']:
                sig = get_sig_digits(hv['value'])
                if sig > bestsig:
                    besthv = hv['value']
                    bestsig = sig
            if bestsig > 0 and is_number(besthv):
                voc = float(besthv)*1.e5/clight
                source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
                add_quantity(name, 'redshift', pretty_num(sqrt((1. + voc)/(1. - voc)) - 1., sig = bestsig), source, kind = 'heliocentric')
        if 'redshift' not in events[name] and has_task('nedd') and 'host' in events[name]:
            reference = "NED-D"
            refurl = "http://ned.ipac.caltech.edu/Library/Distances/"
            for host in events[name]['host']:
                if host['value'] in nedddict:
                    secondarysource = add_source(name, reference = reference, url = refurl, secondary = True)
                    meddist = statistics.median(nedddict[host['value']])
                    redshift = pretty_num(z_at_value(cosmo.comoving_distance, float(meddist) * un.Mpc), sig = get_sig_digits(str(meddist)))
                    add_quantity(name, 'redshift', redshift, secondarysource, kind = 'host')
        if 'maxabsmag' not in events[name] and 'maxappmag' in events[name] and 'lumdist' in events[name]:
            # Find the "best" distance to use for this
            bestsig = 0
            for ld in events[name]['lumdist']:
                sig = get_sig_digits(ld['value'])
                if sig > bestsig:
                    bestld = ld['value']
                    bestsig = sig
            if bestsig > 0 and is_number(bestld) and float(bestld) > 0.:
                source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
                add_quantity(name, 'maxabsmag', pretty_num(float(events[name]['maxappmag'][0]['value']) -
                    5.0*(log10(float(bestld)*1.0e6) - 1.0), sig = bestsig), source)
        if 'redshift' in events[name]:
            # Find the "best" redshift to use for this
            (bestz, bestkind, bestsig) = get_best_redshift(name)
            if bestsig > 0:
                bestz = float(bestz)
                if 'velocity' not in events[name]:
                    source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
                    add_quantity(name, 'velocity', pretty_num(clight/km*((bestz + 1.)**2. - 1.)/
                        ((bestz + 1.)**2. + 1.), sig = bestsig), source, kind = prefkinds[bestkind])
                if bestz > 0.:
                    if 'lumdist' not in events[name]:
                        dl = cosmo.luminosity_distance(bestz)
                        source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
                        add_quantity(name, 'lumdist', pretty_num(dl.value, sig = bestsig), source, kind = prefkinds[bestkind])
                        if 'maxabsmag' not in events[name] and 'maxappmag' in events[name]:
                            source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
                            add_quantity(name, 'maxabsmag', pretty_num(float(events[name]['maxappmag'][0]['value']) -
                                5.0*(log10(dl.to('pc').value) - 1.0), sig = bestsig), source)
                    if 'comovingdist' not in events[name]:
                        dl = cosmo.comoving_distance(bestz)
                        source = add_source(name, bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
                        add_quantity(name, 'comovingdist', pretty_num(dl.value, sig = bestsig), source)
        if 'photometry' in events[name]:
            events[name]['photometry'].sort(key=lambda x: ((float(x['time']) if isinstance(x['time'], str) else
                min([float(y) for y in x['time']])) if 'time' in x else 0.0,
                x['band'] if 'band' in x else '', float(x['magnitude']) if 'magnitude' in x else ''))
        if 'spectra' in events[name] and list(filter(None, ['time' in x for x in events[name]['spectra']])):
            events[name]['spectra'].sort(key=lambda x: (float(x['time']) if 'time' in x else 0.0))
        if 'sources' in events[name]:
            for source in events[name]['sources']:
                if 'bibcode' in source:
                    #First sanitize the bibcode
                    if len(source['bibcode']) != 19:
                        source['bibcode'] = urllib.parse.unquote(unescape(source['bibcode'])).replace('A.A.', 'A&A')
                    if source['bibcode'] in biberrordict:
                        source['bibcode'] = biberrordict[source['bibcode']]

                    if source['bibcode'] not in bibauthordict:
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

                        bibauthordict[bibcode] = unescape(bibcodeauthor).strip()

            for source in events[name]['sources']:
                if 'bibcode' in source and source['bibcode'] in bibauthordict and bibauthordict[source['bibcode']]:
                    source['name'] = bibauthordict[source['bibcode']]
        if 'redshift' in events[name]:
            events[name]['redshift'] = list(sorted(events[name]['redshift'], key=lambda key: frame_priority(key)))
        if 'velocity' in events[name]:
            events[name]['velocity'] = list(sorted(events[name]['velocity'], key=lambda key: frame_priority(key)))
        if 'claimedtype' in events[name]:
            events[name]['claimedtype'] = list(sorted(events[name]['claimedtype'], key=lambda key: ct_priority(name, key)))

        events[name] = OrderedDict(sorted(events[name].items(), key=lambda key: event_attr_priority(key[0])))


def delete_old_event_files():
    # Delete all old event JSON files
    for folder in repofolders:
        filelist = glob.glob("../" + folder + "/*.json") + glob.glob("../" + folder + "/*.json.gz")
        for f in filelist:
            os.remove(f)

def write_all_events(empty = False, gz = False, delete = False):
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
            for r, year in enumerate(repoyears):
                if int(events[name]['discoverdate'][0]['value'].split('/')[0]) <= year:
                    outdir += repofolders[r]
                    break
        else:
            outdir += str(repofolders[0])

        # Delete non-SN events here without IAU designations (those with only banned types)
        if delete:
            deleteevent = False
            nonsnetypes = ['Dwarf Nova', 'Nova', 'QSO', 'AGN', 'CV', 'Galaxy', 'Impostor', 'Imposter', 'Stellar', 'Gal', 'M-star',
                           'AGN / QSO', 'TDE', 'Varstar', 'Star', 'RCrB', 'dK', 'dM', 'SSO', 'YSO', 'LBV', 'BL Lac']
            nonsnetypes = [x.upper() for x in nonsnetypes]
            nonsneprefixes = ('PNVJ', 'PNV J', 'OGLE-2013-NOVA')
            if name.startswith(nonsneprefixes):
                tprint('Deleting ' + name + ', non-SNe prefix.')
                continue
            if 'claimedtype' in events[name] and not (name.startswith('SN') and is_number(name[2:6])):
                for ct in events[name]['claimedtype']:
                    if ct['value'].upper() not in nonsnetypes:
                        deleteevent = False
                        break
                    if ct['value'].upper() in nonsnetypes:
                        deleteevent = True
                if deleteevent:
                    tprint('Deleting ' + name + ' (' + ct['value'] + ').')
                    continue

        jsonstring = json.dumps({name:events[name]}, indent='\t', separators=(',', ':'), ensure_ascii=False)

        path = outdir + '/' + filename + '.json'
        with codecs.open(path, 'w', encoding='utf8') as f:
            f.write(jsonstring)

        if gz:
            if os.path.getsize(path) > 90000000:
                if not args.travis:
                    tprint('Compressing ' + name)
                with open(path, 'rb') as f_in, gzip.open(path + '.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                os.remove(path)
                os.system('cd ' + outdir + '; git rm ' + filename + '.json; git add -f ' + filename + '.json.gz; cd ' + '../scripts')

def null_field(obj, field):
    return obj[field] if field in obj else ''

def copy_to_event(fromname, destname):
    tprint('Copying ' + fromname + ' to event ' + destname)
    newsourcealiases = {}
    keys = list(sorted(events[fromname].keys(), key=lambda key: event_attr_priority(key)))

    if 'sources' in events[fromname]:
        for source in events[fromname]['sources']:
            if 'bibcode' in source:
                newsourcealiases[source['alias']] = add_source(destname, bibcode = source['bibcode'])
            else:
                newsourcealiases[source['alias']] = (add_source(destname,
                    reference = source['name'] if 'name' in source else '',
                    url = source['url'] if 'url' in source else ''))

    for key in keys:
        if key not in ['name', 'sources']:
            for item in events[fromname][key]:
                isd = False
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
                    add_photometry(destname, u_time = null_field(item, "u_time"), time = null_field(item, "time"),
                        e_time = null_field(item, "e_time"), telescope = null_field(item, "telescope"),
                        instrument = null_field(item, "instrument"), band = null_field(item, "band"),
                        magnitude = null_field(item, "magnitude"), e_magnitude = null_field(item, "e_magnitude"),
                        source = sources, upperlimit = null_field(item, "upperlimit"), system = null_field(item, "system"),
                        observatory = null_field(item, "observatory"), observer = null_field(item, "observer"),
                        host = null_field(item, "host"), survey = null_field(item, "survey"))
                elif key == 'spectra':
                    add_spectrum(destname, null_field(item, "waveunit"), null_field(item, "fluxunit"), data = null_field(item, "data"),
                        u_time = null_field(item, "u_time"), time = null_field(item, "time"),
                        instrument = null_field(item, "instrument"), deredshifted = null_field(item, "deredshifted"),
                        dereddened = null_field(item, "dereddened"), errorunit = null_field(item, "errorunit"),
                        source = sources, snr = null_field(item, "snr"),
                        telescope = null_field(item, "telescope"), observer = null_field(item, "observer"),
                        reducer = null_field(item, "reducer"), filename = null_field(item, "filename"),
                        observatory = null_field(item, "observatory"))
                else:
                    add_quantity(destname, key, item['value'], sources, error = null_field(item, "error"),
                        unit = null_field(item, "unit"), kind = null_field(item, "kind"))

def load_event_from_file(name = '', location = '', clean = False, delete = True, append = False):
    if not name and not location:
        raise ValueError('Either event name or location must be specified to load event')

    path = ''
    namepath = ''
    if location:
        path = location
    if name:
        indir = '../'
        for rep in repofolders:
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
                newevent = clean_event(newevent)
            name = next(reversed(newevent))
            if append:
                indir = '../'
                for rep in repofolders:
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
                newevent = OrderedDict([['temp',newevent[name]]])
                copy_to_event('temp', namename)
            else:
                events.update(newevent)

            if args.verbose and not args.travis:
                tprint('Loaded ' + name)

        if 'writeevents' in tasks and delete and namepath:
            os.remove(namepath)
        return name

def clean_event(dirtyevent):
    bibcodes = []
    name = next(reversed(dirtyevent))

    # This is very hacky and is only necessary because we don't have a proper 'Event' object yet.
    events['temp'] = dirtyevent[name]

    if 'name' not in events['temp']:
        events['temp']['name'] = name
    if 'sources' in events['temp']:
        # Rebuild the sources
        newsources = []
        oldsources = events['temp']['sources']
        del(events['temp']['sources'])
        for s, source in enumerate(oldsources):
            if 'bibcode' in source:
                bibcodes.append(source['bibcode'])
                add_source('temp', bibcode = source['bibcode'])
            else:
                add_source('temp', reference = source['name'], url = source['url'])

    # Clean some legacy fields
    if 'aliases' in events['temp'] and isinstance(events['temp']['aliases'], list):
        source = add_source('temp', bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
        for alias in events['temp']['aliases']:
            add_quantity('temp', 'alias', alias, source)
        del(events['temp']['aliases'])

    if ('distinctfrom' in events['temp'] and isinstance(events['temp']['distinctfrom'], list) and
        isinstance(events['temp']['distinctfrom'][0], str)):
        distinctfroms = [x for x in events['temp']['distinctfrom']]
        del(events['temp']['distinctfrom'])
        source = add_source('temp', bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
        for df in distinctfroms:
            add_quantity('temp', 'distinctfrom', df, source)

    if ('errors' in events['temp'] and isinstance(events['temp']['errors'], list) and
        'sourcekind' in events['temp']['errors'][0]):
        source = add_source('temp', bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
        for err in events['temp']['errors']:
            add_quantity('temp', 'error', err['quantity'], source, kind = err['sourcekind'], extra = err['id'])
        del(events['temp']['errors'])

    if not bibcodes:
        add_source('temp', bibcode = oscbibcode, reference = oscname, url = oscurl, secondary = True)
        bibcodes = [oscbibcode]

    for key in list(events['temp'].keys()):
        if key in ['name', 'distinctfrom', 'sources']:
            pass
        elif key == 'photometry':
            for p, photo in enumerate(events['temp']['photometry']):
                if photo['u_time'] == 'JD':
                    events['temp']['photometry'][p]['u_time'] = 'MJD'
                    events['temp']['photometry'][p]['time'] = str(jd_to_mjd(Decimal(photo['time'])))
                if bibcodes and 'source' not in photo:
                    source = add_source('temp', bibcode = bibcodes[0])
                    events['temp']['photometry'][p]['source'] = source
        else:
            for qi, quantity in enumerate(events['temp'][key]):
                if bibcodes and 'source' not in quantity:
                    source = add_source('temp', bibcode = bibcodes[0])
                    events['temp'][key][qi]['source'] = source

    cleanevent = events['temp']
    del (events['temp'])
    return OrderedDict([[name,cleanevent]])

def has_task(task):
    return task in tasks and (not args.update or tasks[task]['update'])

def do_task(checktask, task, quiet = False):
    global currenttask
    dotask = has_task(task) and checktask == task
    if dotask and not quiet:
        currenttask = (tasks[task]['nicename'] if tasks[task]['nicename'] else task).replace('%pre', 'Updating' if args.update else 'Loading')
    return dotask

def journal_events(clear = True):
    if 'writeevents' in tasks:
        write_all_events()
    if clear:
        clear_events()

def clear_events():
    global events
    events = OrderedDict((k, OrderedDict([['name', events[k]['name']]] + ([['alias', events[k]['alias']]] if 'alias' in events[k] else []) + [['stub', True]])) for k in events)

def load_stubs():
    currenttask = 'Loading event stubs'
    files = []

    for rep in repofolders:
        files += glob.glob('../' + rep + "/*.json") + glob.glob('../' + rep + "/*.json.gz")

    #try:
    #    namepath = '../names.min.json'
    #    with open(namepath, 'r') as f:
    #        names = json.loads(f.read(), object_pairs_hook=OrderedDict)
    #    for fi in tq(files):
    #        name = os.path.basename(os.path.splitext(fi)[0])
    #        if name not in names:
    #            name = name.replace("_", "/")
    #        events[name] = OrderedDict(([['name', name], ['alias', [OrderedDict(([['value', x]])) for x in names[name]]], ['stub', True]]))
    #except:
    #    events = OrderedDict()
    for fi in tq(files):
        name = os.path.basename(os.path.splitext(fi)[0])
        name = add_event(name, delete = False, loadifempty = False)
        events[name] = OrderedDict(([['name', events[name]['name']]] + ([['alias', events[name]['alias']]] if 'alias' in events[name] else []) + [['stub', True]]))

path = '../atels.json'
if os.path.isfile(path):
    with open(path, 'r') as f:
        atelsdict = json.loads(f.read(), object_pairs_hook=OrderedDict)
else:
    atelsdict = OrderedDict()
path = '../cbets.json'
if os.path.isfile(path):
    with open(path, 'r') as f:
        cbetsdict = json.loads(f.read(), object_pairs_hook=OrderedDict)
else:
    cbetsdict = OrderedDict()
path = '../iaucs.json'
if os.path.isfile(path):
    with open(path, 'r') as f:
        iaucsdict = json.loads(f.read(), object_pairs_hook=OrderedDict)
else:
    iaucsdict = OrderedDict()

for task in tasks:
    if do_task(task, 'deleteoldevents'):
        delete_old_event_files()

    # Import data provided directly to OSC
    if do_task(task, 'internal'):
        for datafile in tq(sorted(glob.glob("../sne-internal/*.json"), key=lambda s: s.lower())):
            if args.update:
                if not load_event_from_file(location = datafile, clean = True, delete = False, append = True):
                    raise IOError('Failed to find specified file.')
            else:
                if not load_event_from_file(location = datafile, clean = True, delete = False):
                    raise IOError('Failed to find specified file.')
        journal_events()
    
    if do_task(task, 'radio'):
        for datafile in tq(sorted(glob.glob("../sne-external-radio/*.txt"), key=lambda s: s.lower())):
            name = add_event(os.path.basename(datafile).split('.')[0])
            radiosourcedict = OrderedDict()
            with open(datafile, 'r') as f:
                for li, line in enumerate([x.strip() for x in f.read().splitlines()]):
                    if line.startswith('(') and li <= len(radiosourcedict):
                        radiosourcedict[line.split()[0]] = add_source(name, bibcode = line.split()[-1])
                    elif li in [x + len(radiosourcedict) for x in range(3)]:
                        continue
                    else:
                        cols = list(filter(None, line.split()))
                        source = radiosourcedict[cols[6]]
                        add_photometry(name, time = cols[0], frequency = cols[2], u_frequency = 'GHz', fluxdensity = cols[3],
                            e_fluxdensity = cols[4], u_fluxdensity = 'µJy', instrument = cols[5], source = source)
                        add_quantity(name, 'alias', name, source)
        journal_events()
    
    if do_task(task, 'xray'):
        for datafile in tq(sorted(glob.glob("../sne-external-xray/*.txt"), key=lambda s: s.lower())):
            name = add_event(os.path.basename(datafile).split('.')[0])
            with open(datafile, 'r') as f:
                for li, line in enumerate(f.read().splitlines()):
                    if li == 0:
                        source = add_source(name, bibcode = line.split()[-1])
                    elif li in [1,2,3]:
                        continue
                    else:
                        cols = list(filter(None, line.split()))
                        add_photometry(name, time = cols[:2],
                            energy = cols[2:4], u_energy = 'keV', counts = cols[4], flux = cols[6],
                            unabsorbedflux = cols[8], u_flux = 'ergs/s/cm^2',
                            photonindex = cols[15], instrument = cols[17], nhmw = cols[11],
                            upperlimit = (float(cols[5]) < 0), source = source)
                        add_quantity(name, 'alias', name, source)
        journal_events()
    
    #if do_task(task, 'simbad'):
    #    Simbad.list_votable_fields()
    #    customSimbad = Simbad()
    #    customSimbad.add_votable_fields('otype', 'id(opt)')
    #    result = customSimbad.query_object('SN 20[0-9][0-9]*', wildcard=True)
    #    for r, row in enumerate(result):
    #        if row['OTYPE'].decode() != "SN":
    #            continue
    #        name = row["MAIN_ID"].decode()
    #        aliases = Simbad.query_objectids(name)
    #        print(aliases)
    #        if name[:3] == 'SN ':
    #            name = 'SN' + name[3:]
    #        if name[:2] == 'SN' and is_number(name[2:]):
    #            name = name + 'A'
    #        name = add_event(name)
    #    journal_events()
    
    # Import primary data sources from Vizier
    if do_task(task, 'vizier'):
        Vizier.ROW_LIMIT = -1
    
        # 2012ApJS..200...12H
        result = Vizier.get_catalogs("J/ApJS/200/12/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        oldname = ''
        for row in table:
            name = row['SN']
            if is_number(name[:4]):
                name = 'SN' + name
            name = add_event(name)
            source = add_source(name, bibcode = "2012ApJS..200...12H")
            add_quantity(name, 'alias', name, source)
            if '[' not in row['Gal']:
                add_quantity(name, 'host', row['Gal'].replace('_', ' '), source)
            add_quantity(name, 'redshift', str(row['z']), source, kind = 'heliocentric')
            add_quantity(name, 'redshift', str(row['zCMB']), source, kind = 'cmb')
            add_quantity(name, 'ebv', str(row['E_B-V_']), source, error = str(row['e_E_B-V_']) if row['e_E_B-V_'] else '')
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
    
        # 2012ApJ...746...85S
        result = Vizier.get_catalogs("J/ApJ/746/85/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        oldname = ''
        for row in table:
            name = row['Name'].replace('SCP', 'SCP-')
            name = add_event(name)
            source = add_source(name, bibcode = "2012ApJ...746...85S")
            add_quantity(name, 'alias', name, source)
            if row['f_Name']:
                add_quantity(name, 'claimedtype', 'Ia', source)
            if row['z']:
                add_quantity(name, 'redshift', str(row['z']), source, kind = 'spectroscopic')
            else:
                add_quantity(name, 'redshift', str(row['zCl']), source, kind = 'cluster')
            add_quantity(name, 'ebv', str(row['E_B-V_']), source)
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
    
        result = Vizier.get_catalogs("J/ApJ/746/85/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        oldname = ''
        for row in table:
            name = row['Name'].replace('SCP', 'SCP-')
            flux = Decimal(float(row['Flux']))
            if flux <= 0.0:
                continue
            err = Decimal(float(row['e_Flux']))
            zp = Decimal(float(row['Zero']))
            sig = get_sig_digits(str(row['Flux']))+1
            magnitude = pretty_num(zp-Decimal(2.5)*(flux.log10()), sig = sig)
            e_magnitude = pretty_num(Decimal(2.5)*(Decimal(1.0) + err/flux).log10(), sig = sig)
            if float(e_magnitude) > 5.0:
                continue
            name = add_event(name)
            source = add_source(name, bibcode = "2012ApJ...746...85S")
            add_quantity(name, 'alias', name, source)
            add_photometry(name, time = str(row['MJD']), band = row['Filter'], instrument = row['Inst'],
                magnitude = magnitude, e_magnitude = e_magnitude, source = source)
    
        # 2004ApJ...602..571B
        result = Vizier.get_catalogs("J/ApJ/602/571/table8")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        oldname = ''
        for row in table:
            name = 'SN'+row['SN']
            flux = Decimal(float(row['Flux']))
            if flux <= 0.0:
                continue
            err = Decimal(float(row['e_Flux']))
            sig = get_sig_digits(str(row['Flux']))+1
            magnitude = pretty_num(Decimal(25.0)-Decimal(2.5)*(flux.log10()), sig = sig)
            e_magnitude = pretty_num(Decimal(2.5)*(Decimal(1.0) + err/flux).log10(), sig = sig)
            if float(e_magnitude) > 5.0:
                continue
            name = add_event(name)
            source = add_source(name, bibcode = "2004ApJ...602..571B")
            add_quantity(name, 'alias', name, source)
            add_photometry(name, time = str(row['MJD']), band = row['Filt'],
                magnitude = magnitude, e_magnitude = e_magnitude, source = source)
        
        # 2014MNRAS.444.3258M
        result = Vizier.get_catalogs("J/MNRAS/444/3258/SNe")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        oldname = ''
        for row in table:
            name = row['SN']
            if name == oldname:
                continue
            oldname = name
            name = add_event(name)
            source = add_source(name, bibcode = '2014MNRAS.444.3258M')
            add_quantity(name, 'alias', name, source)
            add_quantity(name, 'redshift', str(row['z']), source, kind = 'heliocentric', error = str(row['e_z']))
            add_quantity(name, 'ra', str(row['_RA']), source, unit = 'floatdegrees')
            add_quantity(name, 'dec', str(row['_DE']), source, unit = 'floatdegrees')
        journal_events()
    
        # 2014MNRAS.438.1391P
        result = Vizier.get_catalogs("J/MNRAS/438/1391/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            name = row['SN']
            name = add_event(name)
            source = add_source(name, bibcode = '2014MNRAS.438.1391P')
            add_quantity(name, 'alias', name, source)
            add_quantity(name, 'redshift', str(row['zh']), source, kind = 'heliocentric')
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
        journal_events()
    
        # 2012ApJ...749...18B
        result = Vizier.get_catalogs("J/ApJ/749/18/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            name = row['Name'].replace(' ','')
            name = add_event(name)
            source = add_source(name, bibcode = '2012ApJ...749...18B')
            add_quantity(name, 'alias', name, source)
            mjd = str(astrotime(2450000.+row['JD'], format='jd').mjd)
            band = row['Filt']
            magnitude = str(row['mag'])
            e_magnitude = str(row['e_mag'])
            e_magnitude = '' if e_magnitude == '--' else e_magnitude
            upperlimit = True if row['l_mag'] == '>' else False
            add_photometry(name, time = mjd, band = band, magnitude = magnitude, e_magnitude = e_magnitude, instrument = 'UVOT',
                source = source, upperlimit = upperlimit, telescope = 'Swift')
        journal_events()
    
        # 2010A&A...523A...7G
        result = Vizier.get_catalogs("J/A+A/523/A7/table9")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            name = 'SNLS-' + row['SNLS']
            name = add_event(name)
            source = add_source(name, bibcode = '2010A&A...523A...7G')
            add_quantity(name, 'alias', name, source)
            astrot = astrotime(2450000.+row['Date1'], format='jd').datetime
            add_quantity(name, 'discoverdate', make_date_string(astrot.year, astrot.month, astrot.day), source)
            add_quantity(name, 'ebv', str(row['E_B-V_']), source)
            add_quantity(name, 'redshift', str(row['z']), source, kind = 'heliocentric')
            add_quantity(name, 'claimedtype', (row['Type'].replace('*', '?').replace('SN','')
                .replace('(pec)',' P').replace('Ia? P?', 'Ia P?')), source)
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
        journal_events()
    
        # 2004A&A...415..863G
        result = Vizier.get_catalogs("J/A+A/415/863/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            name = 'SN' + row['SN']
            name = add_event(name)
            source = add_source(name, bibcode = '2004A&A...415..863G')
            add_quantity(name, 'alias', name, source)
            datesplit = row['Date'].split('-')
            add_quantity(name, 'discoverdate', make_date_string(datesplit[0], datesplit[1].lstrip('0'), datesplit[2].lstrip('0')), source)
            add_quantity(name, 'host', 'Abell ' + str(row['Abell']), source)
            add_quantity(name, 'claimedtype', row['Type'], source)
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
            if row['zSN']:
                add_quantity(name, 'redshift', str(row['zSN']), source, kind = 'spectroscopic')
            else:
                add_quantity(name, 'redshift', str(row['zCl']), source, kind = 'cluster')
        journal_events()
    
        # 2008AJ....136.2306H
        result = Vizier.get_catalogs("J/AJ/136/2306/sources")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            name = 'SDSS-II ' + str(row['SNID'])
            name = add_event(name)
            source = add_source(name, bibcode = '2008AJ....136.2306H')
            add_quantity(name, 'alias', name, source)
            add_quantity(name, 'claimedtype', row['SpType'].replace('SN.', '').strip(':'), source)
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
    
        # 2010ApJ...708..661D
        result = Vizier.get_catalogs("J/ApJ/708/661/sn")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            name = row['SN']
            if not name:
                name = 'SDSS-II ' + str(row['SDSS-II'])
            else:
                name = 'SN' + name
            name = add_event(name)
            source = add_source(name, bibcode = '2010ApJ...708..661D')
            add_quantity(name, 'alias', name, source)
            add_quantity(name, 'alias', 'SDSS-II ' + str(row['SDSS-II']), source)
            add_quantity(name, 'claimedtype', 'II P', source)
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
    
        result = Vizier.get_catalogs("J/ApJ/708/661/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            if row['f_SN'] == 'a':
                name = 'SDSS-II ' + str(row['SN'])
            else:
                name = 'SN' + row['SN']
            name = add_event(name)
            source = add_source(name, bibcode = '2010ApJ...708..661D')
            add_quantity(name, 'alias', name, source)
            add_quantity(name, 'redshift', str(row['z']), source, error = str(row['e_z']))
        journal_events()
    
        # 2014ApJ...795...44R
        restpositionerrors = ['SN2010fl']
        result = Vizier.get_catalogs("J/ApJ/795/44/ps1_snIa")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            name = row['SN']
            name = add_event(name)
            source = add_source(name, bibcode = '2014ApJ...795...44R')
            add_quantity(name, 'alias', name, source)
            astrot = astrotime(row['tdisc'], format='mjd').datetime
            add_quantity(name, 'discoverdate',  make_date_string(astrot.year, astrot.month, astrot.day), source)
            add_quantity(name, 'redshift', str(row['z']), source, error = str(row['e_z']), kind = 'heliocentric')
            if name not in restpositionerrors:
                add_quantity(name, 'ra', row['RAJ2000'], source)
                add_quantity(name, 'dec', row['DEJ2000'], source)
            add_quantity(name, 'claimedtype', 'Ia', source)
    
        result = Vizier.get_catalogs("J/ApJ/795/44/table6")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            name = row['SN']
            name = add_event(name)
            source = add_source(name, bibcode = '2014ApJ...795...44R')
            add_quantity(name, 'alias', name, source)
            if row['mag'] != '--':
                add_photometry(name, time = str(row['MJD']), band = row['Filt'], magnitude = str(row['mag']),
                    e_magnitude = str(row['e_mag']), source = source, system = 'AB')
        journal_events()
    
        # 1990A&AS...82..145C
        result = Vizier.get_catalogs("II/189/mag")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
    
        with open('../sne-external/II_189_refs.csv') as f:
            tsvin = csv.reader(f, delimiter='\t', skipinitialspace=True)
            ii189bibdict = {}
            ii189refdict = {}
            for r, row in enumerate(tsvin):
                if row[0] != '0':
                    ii189bibdict[r+1] = row[1]
                else:
                    ii189refdict[r+1] = row[2]
    
        for row in table:
            if row['band'][0] == '(':
                continue
            name = 'SN' + row['SN']
            name = add_event(name)
            source = ''
            secsource = add_source(name, bibcode = '1990A&AS...82..145C', secondary = True)
            mjd = str(jd_to_mjd(Decimal(row['JD'])))
            mag = str(row['m'])
            band = row['band'].strip("'")
            if row['r_m'] in ii189bibdict:
                source = add_source(name, bibcode = ii189bibdict[row['r_m']])
            else:
                source = add_source(name, reference = ii189refdict[row['r_m']])
            add_quantity(name, 'alias', name, source)
    
            add_photometry(name, time = mjd, band = band, magnitude = mag, source = uniq_cdl([source,secsource]))
        journal_events()
    
        # 2014yCat.7272....0G
        result = Vizier.get_catalogs("VII/272/snrs")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
    
        for row in table:
            row = convert_aq_output(row)
            name = ''
            if row["Names"]:
                names = row["Names"].split(',')
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
                name = row["SNR"].strip()
    
            name = add_event(name)
            source = (add_source(name, bibcode = '2014BASI...42...47G') + ',' +
                      add_source(name, reference = 'Galactic SNRs', url = 'https://www.mrao.cam.ac.uk/surveys/snrs/snrs.data.html'))
            add_quantity(name, 'alias', name, source)
    
            add_quantity(name, "alias", row["SNR"].strip(), source)
            add_quantity(name, "alias", 'MWSNR '+row["SNR"].strip('G '), source)
    
            if row["Names"]:
                names = row["Names"].split(',')
                for nam in names:
                    add_quantity(name, "alias", nam.replace('Vela (XYZ)', 'Vela').strip('()').strip(), source)
                    if nam.strip()[:2] == 'SN':
                        add_quantity(name, 'discoverdate', nam.strip()[2:], source)
    
            add_quantity(name, 'host', 'Milky Way', source)
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
        journal_events()
    
        # 2014MNRAS.442..844F
        result = Vizier.get_catalogs("J/MNRAS/442/844/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            name = 'SN' + row['SN']
            name = add_event(name)
            source = add_source(name, bibcode = '2014MNRAS.442..844F')
            add_quantity(name, 'alias', name, source)
            add_quantity(name, 'redshift', str(row['zhost']), source, kind = 'host')
            add_quantity(name, 'ebv', str(row['E_B-V_']), source)
        journal_events()
    
        result = Vizier.get_catalogs("J/MNRAS/442/844/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            name = 'SN' + str(row['SN'])
            name = add_event(name)
            source = add_source(name, bibcode = "2014MNRAS.442..844F")
            add_quantity(name, 'alias', name, source)
            if 'Bmag' in row and is_number(row['Bmag']) and not isnan(float(row['Bmag'])):
                add_photometry(name, time = row['MJD'], band = 'B', magnitude = row['Bmag'], e_magnitude = row['e_Bmag'], source = source)
            if 'Vmag' in row and is_number(row['Vmag']) and not isnan(float(row['Vmag'])):
                add_photometry(name, time = row['MJD'], band = 'V', magnitude = row['Vmag'], e_magnitude = row['e_Vmag'], source = source)
            if 'Rmag' in row and is_number(row['Rmag']) and not isnan(float(row['Rmag'])):
                add_photometry(name, time = row['MJD'], band = 'R', magnitude = row['Rmag'], e_magnitude = row['e_Rmag'], source = source)
            if 'Imag' in row and is_number(row['Imag']) and not isnan(float(row['Imag'])):
                add_photometry(name, time = row['MJD'], band = 'I', magnitude = row['Imag'], e_magnitude = row['e_Imag'], source = source)
        journal_events()
    
        # 2012MNRAS.425.1789S
        result = Vizier.get_catalogs("J/MNRAS/425/1789/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            name = ''.join(row['SimbadName'].split(' '))
            name = add_event(name)
            source = add_source(name, bibcode = '2012MNRAS.425.1789S')
            add_quantity(name, 'alias', name, source)
            add_quantity(name, 'alias', 'SN' + row['SN'], source)
            add_quantity(name, 'host', row['Gal'], source)
            if is_number(row['cz']):
                add_quantity(name, 'redshift', str(round_sig(float(row['cz'])*km/clight, sig = get_sig_digits(str(row['cz'])))), source, kind = 'heliocentric')
            add_quantity(name, 'ebv', str(row['E_B-V_']), source)
        journal_events()
    
        # 2015ApJS..219...13W
        result = Vizier.get_catalogs("J/ApJS/219/13/table3")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            name = u'LSQ' + str(row['LSQ'])
            name = add_event(name)
            source = add_source(name, bibcode = "2015ApJS..219...13W")
            add_quantity(name, 'alias', name, source)
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
            add_quantity(name, 'redshift', row['z'], source, error = row['e_z'], kind = 'heliocentric')
            add_quantity(name, 'ebv', row['E_B-V_'], source)
            add_quantity(name, 'claimedtype', 'Ia', source)
        result = Vizier.get_catalogs("J/ApJS/219/13/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            name = 'LSQ' + row['LSQ']
            name = add_event(name)
            source = add_source(name, bibcode = "2015ApJS..219...13W")
            add_quantity(name, 'alias', name, source)
            add_photometry(name, time = str(jd_to_mjd(Decimal(row['JD']))), instrument = 'QUEST', telescope = 'ESO Schmidt',
                observatory = 'La Silla', band = row['Filt'],
                magnitude = row['mag'], e_magnitude = row['e_mag'], system = "Swope", source = source)
        journal_events()
    
        # 2012Natur.491..228C
        result = Vizier.get_catalogs("J/other/Nat/491.228/tablef1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        name = 'SN2213-1745'
        name = add_event(name)
        source = add_source(name, bibcode = "2012Natur.491..228C")
        add_quantity(name, 'alias', name, source)
        add_quantity(name, 'claimedtype', 'SLSN-R', source)
        for row in table:
            row = convert_aq_output(row)
            if "g_mag" in row and is_number(row["g_mag"]) and not isnan(float(row["g_mag"])):
                add_photometry(name, time = row["MJDg_"], band = "g'", magnitude = row["g_mag"], e_magnitude = row["e_g_mag"], source = source)
            if "r_mag" in row and is_number(row["r_mag"]) and not isnan(float(row["r_mag"])):
                add_photometry(name, time = row["MJDr_"], band = "r'", magnitude = row["r_mag"], e_magnitude = row["e_r_mag"], source = source)
            if "i_mag" in row and is_number(row["i_mag"]) and not isnan(float(row["i_mag"])):
                add_photometry(name, time = row["MJDi_"], band = "i'", magnitude = row["i_mag"], e_magnitude = row["e_i_mag"], source = source)
    
        result = Vizier.get_catalogs("J/other/Nat/491.228/tablef2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        name = 'SN1000+0216'
        name = add_event(name)
        source = add_source(name, bibcode = "2012Natur.491..228C")
        add_quantity(name, 'alias', name, source)
        add_quantity(name, 'claimedtype', 'SLSN-II?', source)
        for row in table:
            row = convert_aq_output(row)
            if "g_mag" in row and is_number(row["g_mag"]) and not isnan(float(row["g_mag"])):
                add_photometry(name, time = row["MJDg_"], band = "g'", magnitude = row["g_mag"], e_magnitude = row["e_g_mag"], source = source)
            if "r_mag" in row and is_number(row["r_mag"]) and not isnan(float(row["r_mag"])):
                add_photometry(name, time = row["MJDr_"], band = "r'", magnitude = row["r_mag"], e_magnitude = row["e_r_mag"], source = source)
            if "i_mag" in row and is_number(row["i_mag"]) and not isnan(float(row["i_mag"])):
                add_photometry(name, time = row["MJDi_"], band = "i'", magnitude = row["i_mag"], e_magnitude = row["e_i_mag"], source = source)
        journal_events()
    
        # 2011Natur.474..484Q
        result = Vizier.get_catalogs("J/other/Nat/474.484/tables1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            name = str(row['Name'])
            name = add_event(name)
            source = add_source(name, bibcode = "2011Natur.474..484Q")
            add_quantity(name, 'alias', name, source)
            add_photometry(name, time = row['MJD'], band = row['Filt'], telescope = row['Tel'], magnitude = row['mag'], e_magnitude = row['e_mag'], source = source)
        journal_events()
    
        # 2011ApJ...736..159G
        result = Vizier.get_catalogs("J/ApJ/736/159/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        name = 'PTF10vdl'
        name = add_event(name)
        source = add_source(name, bibcode = "2011ApJ...736..159G")
        add_quantity(name, 'alias', name, source)
        for row in table:
            row = convert_aq_output(row)
            add_photometry(name, time = str(jd_to_mjd(Decimal(row['JD']))), band = row['Filt'], telescope = row['Tel'], magnitude = row['mag'],
                           e_magnitude = row['e_mag'] if is_number(row['e_mag']) else '', upperlimit = (not is_number(row['e_mag'])), source = source)
        journal_events()
    
        # 2012ApJ...760L..33B
        result = Vizier.get_catalogs("J/ApJ/760/L33/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        name = 'PTF12gzk'
        name = add_event(name)
        source = add_source(name, bibcode = "2012ApJ...760L..33B")
        add_quantity(name, 'alias', name, source)
        for row in table:
            row = convert_aq_output(row)
            # Fixing a typo in VizieR table
            if str(row['JD']) == '2455151.456':
                row['JD'] = '2456151.456'
            add_photometry(name, time = str(jd_to_mjd(Decimal(row['JD']))), band = row['Filt'], telescope = row['Inst'], magnitude = row['mag'],
                           e_magnitude = row['e_mag'], source = source)
        journal_events()
    
        # 2013ApJ...769...39S
        result = Vizier.get_catalogs("J/ApJ/769/39/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        name = 'PS1-12sk'
        name = add_event(name)
        source = add_source(name, bibcode = "2013ApJ...769...39S")
        add_quantity(name, 'alias', name, source)
        for row in table:
            row = convert_aq_output(row)
            instrument = ''
            telescope = ''
            if row['Inst'] == 'RATCam':
                instrument = row['Inst']
            else:
                telescope = row['Inst']
            add_photometry(name, time = row['MJD'], band = row['Filt'], telescope = telescope, instrument = instrument, magnitude = row['mag'],
                           e_magnitude = row['e_mag'] if not row['l_mag'] else '', upperlimit = (row['l_mag'] == '>'), source = source)
        journal_events()
    
        # 2009MNRAS.394.2266P
        # Note: Instrument info available via links in VizieR, can't auto-parse just yet.
        name = 'SN2005cs'
        name = add_event(name)
        source = add_source(name, bibcode = "2009MNRAS.394.2266P")
        add_quantity(name, 'alias', name, source)
        result = Vizier.get_catalogs("J/MNRAS/394/2266/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            if "Umag" in row and is_number(row["Umag"]) and not isnan(float(row["Umag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), band = "U", magnitude = row["Umag"],
                               e_magnitude = (row["e_Umag"] if row['l_Umag'] != '>' else ''), source = source, upperlimit = (row['l_Umag'] == '>'))
            if "Bmag" in row and is_number(row["Bmag"]) and not isnan(float(row["Bmag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), band = "B", magnitude = row["Bmag"],
                               e_magnitude = (row["e_Bmag"] if row['l_Bmag'] != '>' else ''), source = source, upperlimit = (row['l_Bmag'] == '>'))
            if "Vmag" in row and is_number(row["Vmag"]) and not isnan(float(row["Vmag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), band = "V", magnitude = row["Vmag"],
                               e_magnitude = (row["e_Vmag"] if row['l_Vmag'] != '>' else ''), source = source, upperlimit = (row['l_Vmag'] == '>'))
            if "Rmag" in row and is_number(row["Rmag"]) and not isnan(float(row["Rmag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), band = "R", magnitude = row["Rmag"],
                               e_magnitude = (row["e_Rmag"] if row['l_Rmag'] != '>' else ''), source = source, upperlimit = (row['l_Rmag'] == '>'))
            if "Imag" in row and is_number(row["Imag"]) and not isnan(float(row["Imag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), band = "I", magnitude = row["Imag"],
                               e_magnitude = (row["e_Imag"] if row['l_Imag'] != '>' else ''), source = source, upperlimit = (row['l_Imag'] == '>'))
            if "zmag" in row and is_number(row["zmag"]) and not isnan(float(row["zmag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), band = "z", magnitude = row["zmag"],
                               e_magnitude = row["e_zmag"], source = source)
    
        result = Vizier.get_catalogs("J/MNRAS/394/2266/table3")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            if "Bmag" in row and is_number(row["Bmag"]) and not isnan(float(row["Bmag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), band = "B", magnitude = row["Bmag"],
                               e_magnitude = (row["e_Bmag"] if row['l_Bmag'] != '>' else ''), source = source, upperlimit = (row['l_Bmag'] == '>'))
            if "Vmag" in row and is_number(row["Vmag"]) and not isnan(float(row["Vmag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), band = "V", magnitude = row["Vmag"],
                               e_magnitude = (row["e_Vmag"] if row['l_Vmag'] != '>' else ''), source = source, upperlimit = (row['l_Vmag'] == '>'))
            if "Rmag" in row and is_number(row["Rmag"]) and not isnan(float(row["Rmag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), band = "R", magnitude = row["Rmag"],
                               e_magnitude = (row["e_Rmag"] if row['l_Rmag'] != '>' else ''), source = source, upperlimit = (row['l_Rmag'] == '>'))
    
        result = Vizier.get_catalogs("J/MNRAS/394/2266/table4")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            if "Jmag" in row and is_number(row["Jmag"]) and not isnan(float(row["Jmag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), band = "J", magnitude = row["Jmag"],
                               e_magnitude = row["e_Jmag"], source = source)
            if "Hmag" in row and is_number(row["Hmag"]) and not isnan(float(row["Hmag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), band = "H", magnitude = row["Hmag"],
                               e_magnitude = row["e_Hmag"], source = source)
            if "Kmag" in row and is_number(row["Kmag"]) and not isnan(float(row["Kmag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), band = "K", magnitude = row["Kmag"],
                               e_magnitude = row["e_Kmag"], source = source)
        journal_events()
    
        # 2013AJ....145...99A
        result = Vizier.get_catalogs("J/AJ/145/99/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        name = 'SN2003ie'
        name = add_event(name)
        source = add_source(name, bibcode = "2013AJ....145...99A")
        add_quantity(name, 'alias', name, source)
        for row in table:
            row = convert_aq_output(row)
            if "Bmag" in row and is_number(row["Bmag"]) and not isnan(float(row["Bmag"])):
                add_photometry(name, time = row["MJD"], band = "B", magnitude = row["Bmag"],
                               e_magnitude = row["e_Bmag"] if not row["l_Bmag"] else '',
                               upperlimit = (row['l_Bmag'] == '>'), source = source)
            if "Vmag" in row and is_number(row["Vmag"]) and not isnan(float(row["Vmag"])):
                add_photometry(name, time = row["MJD"], band = "V", magnitude = row["Vmag"],
                               e_magnitude = row["e_Vmag"] if is_number(row["e_Vmag"]) else '',
                               upperlimit = (not is_number(row["e_Vmag"])), source = source)
            if "Rmag" in row and is_number(row["Rmag"]) and not isnan(float(row["Rmag"])):
                add_photometry(name, time = row["MJD"], band = "R", magnitude = row["Rmag"],
                               e_magnitude = row["e_Rmag"] if not row["l_Rmag"] else '',
                               upperlimit = (row['l_Rmag'] == '>'), source = source)
            if "Imag" in row and is_number(row["Imag"]) and not isnan(float(row["Imag"])):
                add_photometry(name, time = row["MJD"], band = "I", magnitude = row["Imag"],
                               e_magnitude = row["e_Imag"], source = source)
        journal_events()
    
        # 2011ApJ...729..143C
        name = 'SN2008am'
        name = add_event(name)
        source = add_source(name, bibcode = "2011ApJ...729..143C")
        add_quantity(name, 'alias', name, source)
    
        result = Vizier.get_catalogs("J/ApJ/729/143/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            add_photometry(name, time = row['MJD'], band = 'ROTSE', telescope = 'ROTSE', magnitude = row['mag'],
                           e_magnitude = row['e_mag'] if not row['l_mag'] else '', upperlimit = (row['l_mag'] == '<'), source = source)
    
        result = Vizier.get_catalogs("J/ApJ/729/143/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            if "Jmag" in row and is_number(row["Jmag"]) and not isnan(float(row["Jmag"])):
                add_photometry(name, time = row["MJD"], telescope = "PAIRITEL", band = "J", magnitude = row["Jmag"],
                               e_magnitude = row["e_Jmag"], source = source)
            if "Hmag" in row and is_number(row["Hmag"]) and not isnan(float(row["Hmag"])):
                add_photometry(name, time = row["MJD"], telescope = "PAIRITEL", band = "H", magnitude = row["Hmag"],
                               e_magnitude = row["e_Hmag"], source = source)
            if "Ksmag" in row and is_number(row["Ksmag"]) and not isnan(float(row["Ksmag"])):
                add_photometry(name, time = row["MJD"], telescope = "PAIRITEL", band = "Ks", magnitude = row["Ksmag"],
                               e_magnitude = row["e_Ksmag"], source = source)
    
        result = Vizier.get_catalogs("J/ApJ/729/143/table4")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            add_photometry(name, time = row['MJD'], band = row['Filt'], telescope = 'P60', magnitude = row['mag'],
                           e_magnitude = row['e_mag'], source = source)
    
        result = Vizier.get_catalogs("J/ApJ/729/143/table5")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            add_photometry(name, time = row['MJD'], band = row['Filt'], instrument = 'UVOT', telescope = 'Swift', magnitude = row['mag'],
                           e_magnitude = row['e_mag'], source = source)
        journal_events()
    
        # 2011ApJ...728...14P
        name = 'SN2009bb'
        name = add_event(name)
        source = add_source(name, bibcode = "2011ApJ...728...14P")
        add_quantity(name, 'alias', name, source)
    
        result = Vizier.get_catalogs("J/ApJ/728/14/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            if "Bmag" in row and is_number(row["Bmag"]) and not isnan(float(row["Bmag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), telescope = row['Tel'], band = "B", magnitude = row["Bmag"],
                               e_magnitude = row["e_Bmag"], source = source)
            if "Vmag" in row and is_number(row["Vmag"]) and not isnan(float(row["Vmag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), telescope = row['Tel'], band = "V", magnitude = row["Vmag"],
                               e_magnitude = row["e_Vmag"], source = source)
            if "Rmag" in row and is_number(row["Rmag"]) and not isnan(float(row["Rmag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), telescope = row['Tel'], band = "R", magnitude = row["Rmag"],
                               e_magnitude = row["e_Rmag"], source = source)
            if "Imag" in row and is_number(row["Imag"]) and not isnan(float(row["Imag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), telescope = row['Tel'], band = "I", magnitude = row["Imag"],
                               e_magnitude = row["e_Imag"], source = source)
    
        result = Vizier.get_catalogs("J/ApJ/728/14/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            if "u_mag" in row and is_number(row["u_mag"]) and not isnan(float(row["u_mag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), telescope = row['Tel'], band = "u'", magnitude = row["u_mag"],
                               e_magnitude = row["e_u_mag"], source = source)
            if "g_mag" in row and is_number(row["g_mag"]) and not isnan(float(row["g_mag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), telescope = row['Tel'], band = "g'", magnitude = row["g_mag"],
                               e_magnitude = row["e_g_mag"], source = source)
            if "r_mag" in row and is_number(row["r_mag"]) and not isnan(float(row["r_mag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), telescope = row['Tel'], band = "r'", magnitude = row["r_mag"],
                               e_magnitude = row["e_r_mag"], source = source)
            if "i_mag" in row and is_number(row["i_mag"]) and not isnan(float(row["i_mag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), telescope = row['Tel'], band = "i'", magnitude = row["i_mag"],
                               e_magnitude = row["e_i_mag"], source = source)
            if "z_mag" in row and is_number(row["z_mag"]) and not isnan(float(row["z_mag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), telescope = row['Tel'], band = "z'", magnitude = row["z_mag"],
                               e_magnitude = row["e_z_mag"], source = source)
    
        result = Vizier.get_catalogs("J/ApJ/728/14/table3")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            if "Ymag" in row and is_number(row["Ymag"]) and not isnan(float(row["Ymag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), instrument = row['Inst'], band = "Y", magnitude = row["Ymag"],
                               e_magnitude = row["e_Ymag"], source = source)
            if "Jmag" in row and is_number(row["Jmag"]) and not isnan(float(row["Jmag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), instrument = row['Inst'], band = "J", magnitude = row["Jmag"],
                               e_magnitude = row["e_Jmag"], source = source)
            if "Hmag" in row and is_number(row["Hmag"]) and not isnan(float(row["Hmag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), instrument = row['Inst'], band = "H", magnitude = row["Hmag"],
                               e_magnitude = row["e_Hmag"], source = source)
        journal_events()
    
        # 2011PAZh...37..837T
        name = 'SN2009nr'
        name = add_event(name)
        source = add_source(name, bibcode = "2011PAZh...37..837T")
        add_quantity(name, 'alias', name, source)
    
        result = Vizier.get_catalogs("J/PAZh/37/837/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            mjd = str(jd_to_mjd(Decimal(row["JD"]) + 2455000))
            if "Umag" in row and is_number(row["Umag"]) and not isnan(float(row["Umag"])):
                add_photometry(name, time = mjd, telescope = row['Tel'], band = "U", magnitude = row["Umag"],
                               e_magnitude = row["e_Umag"], source = source)
            if "Bmag" in row and is_number(row["Bmag"]) and not isnan(float(row["Bmag"])):
                add_photometry(name, time = mjd, telescope = row['Tel'], band = "B", magnitude = row["Bmag"],
                               e_magnitude = row["e_Bmag"], source = source)
            if "Vmag" in row and is_number(row["Vmag"]) and not isnan(float(row["Vmag"])):
                add_photometry(name, time = mjd, telescope = row['Tel'], band = "V", magnitude = row["Vmag"],
                               e_magnitude = row["e_Vmag"], source = source)
            if "Rmag" in row and is_number(row["Rmag"]) and not isnan(float(row["Rmag"])):
                add_photometry(name, time = mjd, telescope = row['Tel'], band = "R", magnitude = row["Rmag"],
                               e_magnitude = row["e_Rmag"], source = source)
            if "Imag" in row and is_number(row["Imag"]) and not isnan(float(row["Imag"])):
                add_photometry(name, time = mjd, telescope = row['Tel'], band = "I", magnitude = row["Imag"],
                               e_magnitude = row["e_Imag"], source = source)
        journal_events()
    
        # 2013MNRAS.433.1871B
        name = 'SN2012aw'
        name = add_event(name)
        source = add_source(name, bibcode = "2013MNRAS.433.1871B")
        add_quantity(name, 'alias', name, source)
    
        result = Vizier.get_catalogs("J/MNRAS/433/1871/table3a")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            mjd = str(jd_to_mjd(Decimal(row["JD"]) + 2456000))
            if "Umag" in row and is_number(row["Umag"]) and not isnan(float(row["Umag"])):
                add_photometry(name, time = mjd, telescope = row['Tel'], band = "U", magnitude = row["Umag"],
                               e_magnitude = row["e_Umag"], source = source)
            if "Bmag" in row and is_number(row["Bmag"]) and not isnan(float(row["Bmag"])):
                add_photometry(name, time = mjd, telescope = row['Tel'], band = "B", magnitude = row["Bmag"],
                               e_magnitude = row["e_Bmag"], source = source)
            if "Vmag" in row and is_number(row["Vmag"]) and not isnan(float(row["Vmag"])):
                add_photometry(name, time = mjd, telescope = row['Tel'], band = "V", magnitude = row["Vmag"],
                               e_magnitude = row["e_Vmag"], source = source)
            if "Rcmag" in row and is_number(row["Rcmag"]) and not isnan(float(row["Rcmag"])):
                add_photometry(name, time = mjd, telescope = row['Tel'], band = "Rc", magnitude = row["Rcmag"],
                               e_magnitude = row["e_Rcmag"], source = source)
            if "Icmag" in row and is_number(row["Icmag"]) and not isnan(float(row["Icmag"])):
                add_photometry(name, time = mjd, telescope = row['Tel'], band = "Ic", magnitude = row["Icmag"],
                               e_magnitude = row["e_Icmag"], source = source)
    
        result = Vizier.get_catalogs("J/MNRAS/433/1871/table3b")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            mjd = str(jd_to_mjd(Decimal(row["JD"]) + 2456000))
            if "gmag" in row and is_number(row["gmag"]) and not isnan(float(row["gmag"])):
                add_photometry(name, time = mjd, telescope = row['Tel'], band = "g", magnitude = row["gmag"],
                               e_magnitude = row["e_gmag"], source = source)
            if "rmag" in row and is_number(row["rmag"]) and not isnan(float(row["rmag"])):
                add_photometry(name, time = mjd, telescope = row['Tel'], band = "r", magnitude = row["rmag"],
                               e_magnitude = row["e_rmag"], source = source)
            if "imag" in row and is_number(row["imag"]) and not isnan(float(row["imag"])):
                add_photometry(name, time = mjd, telescope = row['Tel'], band = "i", magnitude = row["imag"],
                               e_magnitude = row["e_imag"], source = source)
            if "zmag" in row and is_number(row["zmag"]) and not isnan(float(row["zmag"])):
                add_photometry(name, time = mjd, telescope = row['Tel'], band = "z", magnitude = row["zmag"],
                               e_magnitude = row["e_zmag"], source = source)
        journal_events()
    
        # 2014AJ....148....1Z
        name = 'SN2012fr'
        name = add_event(name)
        source = add_source(name, bibcode = "2014AJ....148....1Z")
        add_quantity(name, 'alias', name, source)
    
        result = Vizier.get_catalogs("J/AJ/148/1/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            mjd = row['MJD']
            if "Bmag" in row and is_number(row["Bmag"]) and not isnan(float(row["Bmag"])):
                add_photometry(name, time = mjd, telescope = "LJT", instrument = "YFOSC", band = "B", magnitude = row["Bmag"],
                               e_magnitude = row["e_Bmag"], source = source)
            if "Vmag" in row and is_number(row["Vmag"]) and not isnan(float(row["Vmag"])):
                add_photometry(name, time = mjd, telescope = "LJT", instrument = "YFOSC", band = "V", magnitude = row["Vmag"],
                               e_magnitude = row["e_Vmag"], source = source)
            if "Rmag" in row and is_number(row["Rmag"]) and not isnan(float(row["Rmag"])):
                add_photometry(name, time = mjd, telescope = "LJT", instrument = "YFOSC", band = "R", magnitude = row["Rmag"],
                               e_magnitude = row["e_Rmag"], source = source)
            if "Imag" in row and is_number(row["Imag"]) and not isnan(float(row["Imag"])):
                add_photometry(name, time = mjd, telescope = "LJT", instrument = "YFOSC", band = "I", magnitude = row["Imag"],
                               e_magnitude = row["e_Imag"], source = source)
    
        result = Vizier.get_catalogs("J/AJ/148/1/table3")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            mjd = row['MJD']
            if "Umag" in row and is_number(row["Umag"]) and not isnan(float(row["Umag"])):
                add_photometry(name, time = mjd, telescope = "Swift", instrument = "UVOT", band = "U", magnitude = row["Umag"],
                               e_magnitude = row["e_Umag"], source = source)
            if "Bmag" in row and is_number(row["Bmag"]) and not isnan(float(row["Bmag"])):
                add_photometry(name, time = mjd, telescope = "Swift", instrument = "UVOT", band = "B", magnitude = row["Bmag"],
                               e_magnitude = row["e_Bmag"], source = source)
            if "Vmag" in row and is_number(row["Vmag"]) and not isnan(float(row["Vmag"])):
                add_photometry(name, time = mjd, telescope = "Swift", instrument = "UVOT", band = "V", magnitude = row["Vmag"],
                               e_magnitude = row["e_Vmag"], source = source)
            if "UVW1" in row and is_number(row["UVW1"]) and not isnan(float(row["UVW1"])):
                add_photometry(name, time = mjd, telescope = "Swift", instrument = "UVOT", band = "W1", magnitude = row["UVW1"],
                               e_magnitude = row["e_UVW1"], source = source)
            if "UVW2" in row and is_number(row["UVW2"]) and not isnan(float(row["UVW2"])):
                add_photometry(name, time = mjd, telescope = "Swift", instrument = "UVOT", band = "W2", magnitude = row["UVW2"],
                               e_magnitude = row["e_UVW2"], source = source)
            if "UVM2" in row and is_number(row["UVM2"]) and not isnan(float(row["UVM2"])):
                add_photometry(name, time = mjd, telescope = "Swift", instrument = "UVOT", band = "M2", magnitude = row["UVM2"],
                               e_magnitude = row["e_UVM2"], source = source)
    
        result = Vizier.get_catalogs("J/AJ/148/1/table5")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            mjd = row['MJD']
            if "Bmag" in row and is_number(row["Bmag"]) and not isnan(float(row["Bmag"])):
                add_photometry(name, time = mjd, telescope = "LJT", band = "B", magnitude = row["Bmag"],
                               e_magnitude = row["e_Bmag"], source = source)
            if "Vmag" in row and is_number(row["Vmag"]) and not isnan(float(row["Vmag"])):
                add_photometry(name, time = mjd, telescope = "LJT", band = "V", magnitude = row["Vmag"],
                               e_magnitude = row["e_Vmag"], source = source)
            if "Rmag" in row and is_number(row["Rmag"]) and not isnan(float(row["Rmag"])):
                add_photometry(name, time = mjd, telescope = "LJT", band = "R", magnitude = row["Rmag"],
                               e_magnitude = row["e_Rmag"], source = source)
            if "Imag" in row and is_number(row["Imag"]) and not isnan(float(row["Imag"])):
                add_photometry(name, time = mjd, telescope = "LJT", band = "I", magnitude = row["Imag"],
                               e_magnitude = row["e_Imag"], source = source)
        journal_events()
    
        # 2015ApJ...805...74B
        name = 'SN2014J'
        name = add_event(name)
        source = add_source(name, bibcode = "2014AJ....148....1Z")
        add_quantity(name, 'alias', name, source)
    
        result = Vizier.get_catalogs("J/ApJ/805/74/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            mjd = row['MJD']
            if "mag" in row and is_number(row["mag"]) and not isnan(float(row["mag"])):
                add_photometry(name, time = mjd, telescope = "Swift", instrument = "UVOT", band = row["Filt"], magnitude = row["mag"],
                               e_magnitude = row["e_mag"], source = source)
            elif "maglim" in row and is_number(row["maglim"]) and not isnan(float(row["maglim"])):
                add_photometry(name, time = mjd, telescope = "Swift", instrument = "UVOT", band = row["Filt"], magnitude = row["maglim"],
                               upperlimit = True, source = source)
        journal_events()
    
        # 2011ApJ...741...97D
        result = Vizier.get_catalogs("J/ApJ/741/97/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            name = str(row['SN'])
            name = add_event(name)
            source = add_source(name, bibcode = "2011ApJ...741...97D")
            add_quantity(name, 'alias', name, source)
            add_photometry(name, time = str(jd_to_mjd(Decimal(row['JD']))), band = row['Filt'], magnitude = row['mag'],
                           e_magnitude = row['e_mag'] if is_number(row['e_mag']) else '', upperlimit = (not is_number(row['e_mag'])), source = source)
        journal_events()
    
        # 2015MNRAS.448.1206M
        # Note: Photometry from two SN can also be added from this source.
        result = Vizier.get_catalogs("J/MNRAS/448/1206/table3")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            name = str(row['Name'])
            name = add_event(name)
            source = add_source(name, bibcode = "2015MNRAS.448.1206M")
            add_quantity(name, 'alias', name, source)
            add_quantity(name, 'discoverdate', '20' + name[4:6], source)
            add_quantity(name, 'ra', row['RAJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'dec', row['DEJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'redshift', row['zsp'], source, kind = 'spectroscopic')
            add_quantity(name, 'maxappmag', row['rP1mag'], source, error = row['e_rP1mag'])
            add_quantity(name, 'maxband', 'r', source)
            add_quantity(name, 'claimedtype', 'Ia', source)
        result = Vizier.get_catalogs("J/MNRAS/448/1206/table4")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            name = str(row['Name'])
            name = add_event(name)
            source = add_source(name, bibcode = "2015MNRAS.448.1206M")
            add_quantity(name, 'alias', name, source)
            add_quantity(name, 'discoverdate', '20' + name[4:6], source)
            add_quantity(name, 'ra', row['RAJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'dec', row['DEJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'redshift', row['zph'], source, error = row['e_zph'], kind = 'photometric')
            add_quantity(name, 'maxappmag', row['rP1mag'], source, error = row['e_rP1mag'])
            add_quantity(name, 'maxband', 'r', source)
            add_quantity(name, 'claimedtype', 'Ia?', source)
        result = Vizier.get_catalogs("J/MNRAS/448/1206/table5")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            name = str(row['Name'])
            name = add_event(name)
            source = add_source(name, bibcode = "2015MNRAS.448.1206M")
            add_quantity(name, 'alias', name, source)
            add_quantity(name, 'discoverdate', '20' + name[4:6], source)
            add_quantity(name, 'ra', row['RAJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'dec', row['DEJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'redshift', row['zsp'], source, kind = 'spectroscopic')
            add_quantity(name, 'maxappmag', row['rP1mag'], source, error = row['e_rP1mag'])
            add_quantity(name, 'maxband', 'r', source)
            add_quantity(name, 'claimedtype', row['Type'], source)
        result = Vizier.get_catalogs("J/MNRAS/448/1206/table6")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            name = str(row['Name'])
            name = add_event(name)
            source = add_source(name, bibcode = "2015MNRAS.448.1206M")
            add_quantity(name, 'alias', name, source)
            add_quantity(name, 'discoverdate', '20' + name[4:6], source)
            add_quantity(name, 'ra', row['RAJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'dec', row['DEJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'maxappmag', row['rP1mag'], source, error = row['e_rP1mag'])
            add_quantity(name, 'maxband', 'r', source)
            add_quantity(name, 'claimedtype', row['Type'], source)
        result = Vizier.get_catalogs("J/MNRAS/448/1206/tablea2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            name = str(row['Name'])
            name = add_event(name)
            source = add_source(name, bibcode = "2015MNRAS.448.1206M")
            add_quantity(name, 'alias', name, source)
            add_quantity(name, 'discoverdate', '20' + name[4:6], source)
            add_quantity(name, 'ra', row['RAJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'dec', row['DEJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'maxappmag', row['rP1mag'], source, error = row['e_rP1mag'])
            add_quantity(name, 'maxband', 'r', source)
            add_quantity(name, 'claimedtype', row['Typesoft']+'?', source)
            add_quantity(name, 'claimedtype', row['Typepsnid']+'?', source)
        result = Vizier.get_catalogs("J/MNRAS/448/1206/tablea3")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            row = convert_aq_output(row)
            name = str(row['Name'])
            name = add_event(name)
            source = add_source(name, bibcode = "2015MNRAS.448.1206M")
            add_quantity(name, 'alias', name, source)
            add_quantity(name, 'discoverdate', '20' + name[4:6], source)
            add_quantity(name, 'ra', row['RAJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'dec', row['DEJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'maxappmag', row['rP1mag'], source, error = row['e_rP1mag'])
            add_quantity(name, 'maxband', 'r', source)
            add_quantity(name, 'claimedtype', 'Candidate', source)
        journal_events()
    
        # 2012AJ....143..126B
        result = Vizier.get_catalogs("J/AJ/143/126/table4")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in table:
            if not row['Wcl'] or row['Wcl'] == 'N':
                continue
            row = convert_aq_output(row)
            name = str(row['SN']).replace(' ', '')
            name = add_event(name)
            source = add_source(name, bibcode = "2012AJ....143..126B")
            add_quantity(name, 'alias', name, source)
            add_quantity(name, 'claimedtype', 'Ia-' + row['Wcl'], source)
        journal_events()
    
    if do_task(task, 'donations'):
        # Nicholl 04-01-16 donation
        with open("../sne-external/Nicholl-04-01-16/bibcodes.json", 'r') as f:
            bcs = json.loads(f.read())
    
        for datafile in sorted(glob.glob("../sne-external/Nicholl-04-01-16/*.txt"), key=lambda s: s.lower()):
            name = os.path.basename(datafile).split('_')[0]
            name = add_event(name)
            bibcode = ''
            for bc in bcs:
                if name in bcs[bc]:
                    bibcode = bc
            if not bibcode:
                raise(ValueError('Bibcode not found!'))
            source = add_source(name, bibcode = bibcode)
            add_quantity(name, 'alias', name, source)
            with open(datafile,'r') as f:
                tsvin = csv.reader(f, delimiter='\t', skipinitialspace=True)
                for r, rrow in enumerate(tsvin):
                    row = list(filter(None, rrow))
                    if not row:
                        continue
                    if row[0][0] == '#' and row[0] != '#MJD':
                        continue
                    if row[0] == '#MJD':
                        bands = [x for x in row[1:] if x and 'err' not in x]
                        continue
                    mjd = row[0]
                    if not is_number(mjd):
                        continue
                    for v, val in enumerate(row[1::2]):
                        upperlimit = ''
                        if '>' in val:
                            upperlimit = True
                        mag = val.strip('>')
                        if not is_number(mag) or isnan(float(mag)) or float(mag) > 90.0:
                            continue
                        err = ''
                        if is_number(row[2*v+2]) and not isnan(float(row[2*v+2])):
                            err = row[2*v+2]
                        add_photometry(name, time = mjd, band = bands[v], magnitude = mag,
                            e_magnitude = err, upperlimit = upperlimit, source = source)
        journal_events()
    
        # Maggi 04-11-16 donation (MC SNRs)
        with open('../sne-external/Maggi-04-11-16/LMCSNRs_OpenSNe.csv') as f:
            tsvin = csv.reader(f, delimiter=',')
            for row in tsvin:
                name = 'MCSNR ' + row[0]
                name = add_event(name)
                source = add_source(name, bibcode = '2016A&A...585A.162M')
                add_quantity(name, 'alias', name, source)
                if row[1] != 'noname':
                    add_quantity(name, "alias", row[1], source)
                add_quantity(name, 'ra', row[2], source)
                add_quantity(name, 'dec', row[3], source)
                add_quantity(name, 'host', 'LMC', source)
                if row[4] == '1':
                    add_quantity(name, 'claimedtype', 'Ia', source)
                elif row[4] == '2':
                    add_quantity(name, 'claimedtype', 'CC', source)
        with open('../sne-external/Maggi-04-11-16/SMCSNRs_OpenSNe.csv') as f:
            tsvin = csv.reader(f, delimiter=',')
            for row in tsvin:
                name = 'MCSNR ' + row[0]
                name = add_event(name)
                source = add_source(name, reference = 'Pierre Maggi')
                add_quantity(name, 'alias', name, source)
                add_quantity(name, "alias", row[1], source)
                add_quantity(name, "alias", row[2], source)
                add_quantity(name, 'ra', row[3], source)
                add_quantity(name, 'dec', row[4], source)
                add_quantity(name, 'host', 'SMC', source)
        journal_events()
    
        # Galbany 04-18-16 donation
        folders = next(os.walk('../sne-external/galbany-04-18-16/'))[1]
        bibcode = '2016AJ....151...33G'
        for folder in folders:
            infofiles = glob.glob("../sne-external/galbany-04-18-16/" + folder + "/*.info")
            photfiles = glob.glob("../sne-external/galbany-04-18-16/" + folder + "/*.out*")
    
            zhel = ''
            zcmb = ''
            zerr = ''
            for path in infofiles:
                with open(path, 'r') as f:
                    lines = f.read().splitlines()
                    for line in lines:
                        splitline = line.split(':')
                        field = splitline[0].strip().lower()
                        value = splitline[1].strip()
                        if field == 'name':
                            name = value[:6].upper() + (value[6].upper() if len(value) == 7 else value[6:])
                            name = add_event(name)
                            source = add_source(name, bibcode = bibcode)
                            add_quantity(name, 'alias', name, source)
                        elif field == 'type':
                            claimedtype = value.replace('SN', '')
                            add_quantity(name, 'claimedtype', claimedtype, source)
                        elif field == 'zhel':
                            zhel = value
                        elif field == 'redshift_error':
                            zerr = value
                        elif field == 'zcmb':
                            zcmb = value
                        elif field == 'ra':
                            add_quantity(name, 'ra', value, source, unit = 'floatdegrees')
                        elif field == 'dec':
                            add_quantity(name, 'dec', value, source, unit = 'floatdegrees')
                        elif field == 'host':
                            add_quantity(name, 'host', value.replace('- ', '-').replace('G ', 'G'), source)
                        elif field == 'e(b-v)_mw':
                            add_quantity(name, 'ebv', value, source)
    
            add_quantity(name, 'redshift', zhel, source, error = zerr, kind = 'heliocentric')
            add_quantity(name, 'redshift', zcmb, source, error = zerr, kind = 'cmb')
    
            for path in photfiles:
                with open(path, 'r') as f:
                    band = ''
                    lines = f.read().splitlines()
                    for li, line in enumerate(lines):
                        if li in [0, 2, 3]:
                            continue
                        if li == 1:
                            band = line.split(':')[-1].strip()
                        else:
                            cols = list(filter(None, line.split()))
                            if not cols:
                                continue
                            add_photometry(name, time = cols[0], magnitude = cols[1], e_magnitude = cols[2],
                                           band = band, system = cols[3], telescope = cols[4], source = source)
        journal_events()
    
    if do_task(task, 'pessto-dr1'):
        with open("../sne-external/PESSTO_MPHOT.csv", 'r') as f:
            tsvin = csv.reader(f, delimiter=',')
            for ri, row in enumerate(tsvin):
                if ri == 0:
                    bands = [x.split('_')[0] for x in row[3::2]]
                    systems = [x.split('_')[1].capitalize().replace('Ab', 'AB') for x in row[3::2]]
                    continue
                name = row[1]
                name = add_event(name)
                source = add_source(name, bibcode = "2015A&A...579A..40S")
                add_quantity(name, 'alias', name, source)
                for hi, ci in enumerate(range(3,len(row)-1,2)):
                    if not row[ci]:
                        continue
                    add_photometry(name, time = row[2], magnitude = row[ci], e_magnitude = row[ci+1],
                        band = bands[hi], system = systems[hi], telescope = 'Swift' if systems[hi] == 'Swift' else '',
                        source = source)
        journal_events()
    
    if do_task(task, 'scp'):
        with open("../sne-external/SCP09.csv", 'r') as f:
            tsvin = csv.reader(f, delimiter=',')
            for ri, row in enumerate(tq(tsvin)):
                if ri == 0:
                    continue
                name = row[0].replace('SCP', 'SCP-')
                name = add_event(name)
                source = add_source(name, reference = 'Supernova Cosmology Project', url = 'http://supernova.lbl.gov/2009ClusterSurvey/')
                add_quantity(name, 'alias', name, source)
                if row[1]:
                    add_quantity(name, 'alias', row[1], source)
                if row[2]:
                    add_quantity(name, 'redshift', row[2], source, kind = 'spectroscopic' if row[3] == 'sn' else 'host')
                if row[4]:
                    add_quantity(name, 'redshift', row[2], source, kind = 'cluster')
                if row[6]:
                    claimedtype = row[6].replace('SN ', '')
                    kind = ('spectroscopic/light curve' if 'a' in row[7] and 'c' in row[7] else
                        'spectroscopic' if 'a' in row[7] else 'light curve' if 'c' in row[7] else '')
                    if claimedtype != '?':
                        add_quantity(name, 'claimedtype', claimedtype, source, kind = kind)
        journal_events()
    
    if do_task(task, 'ascii'):
        # 2006ApJ...645..841N
        with open("../sne-external/2006ApJ...645..841N-table3.csv", 'r') as f:
            tsvin = csv.reader(f, delimiter=',')
            for ri, row in enumerate(tsvin):
                name = 'SNLS-' + row[0]
                name = add_event(name)
                source = add_source(name, bibcode = '2006ApJ...645..841N')
                add_quantity(name, 'alias', name, source)
                add_quantity(name, 'redshift', row[1], source, kind = 'spectroscopic')
                astrot = astrotime(float(row[4]) + 2450000., format = 'jd').datetime
                add_quantity(name, 'discoverdate', make_date_string(astrot.year, astrot.month, astrot.day), source)
        journal_events()
    
        # Anderson 2014
        for datafile in sorted(glob.glob("../sne-external/SNII_anderson2014/*.dat"), key=lambda s: s.lower()):
            basename = os.path.basename(datafile)
            if not is_number(basename[:2]):
                continue
            if basename == '0210_V.dat':
                name = 'SN0210'
            else:
                name = ('SN20' if int(basename[:2]) < 50 else 'SN19') + basename.split('_')[0]
            name = add_event(name)
            source = add_source(name, bibcode = '2014ApJ...786...67A')
            add_quantity(name, 'alias', name, source)
    
            if name in ['SN1999ca','SN2003dq','SN2008aw']:
                system = 'Swope'
            else:
                system = 'Landolt'
    
            with open(datafile,'r') as f:
                tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)
                for row in tsvin:
                    if not row[0]:
                        continue
                    add_photometry(name, time = str(jd_to_mjd(Decimal(row[0]))), band = 'V', magnitude = row[1], e_magnitude = row[2], system = system, source = source)
        journal_events()
    
        # stromlo
        stromlobands = ['B','V','R','I','VM','RM']
        with open('../sne-external/J_A+A_415_863-1/photometry.csv', 'r') as f:
            tsvin = csv.reader(f, delimiter=',')
            for row in tsvin:
                name = row[0]
                name = add_event(name)
                source = add_source(name, bibcode = "2004A&A...415..863G")
                add_quantity(name, 'alias', name, source)
                mjd = str(jd_to_mjd(Decimal(row[1])))
                for ri, ci in enumerate(range(2,len(row),3)):
                    if not row[ci]:
                        continue
                    band = stromlobands[ri]
                    upperlimit = True if (not row[ci+1] and row[ci+2]) else False
                    e_upper_magnitude = str(abs(Decimal(row[ci+1]))) if row[ci+1] else ''
                    e_lower_magnitude = str(abs(Decimal(row[ci+2]))) if row[ci+2] else ''
                    add_photometry(name, time = mjd, band = band, magnitude = row[ci],
                        e_upper_magnitude = e_upper_magnitude, e_lower_magnitude = e_lower_magnitude,
                        upperlimit = upperlimit, telescope = 'MSSSO 1.3m' if band in ['VM', 'RM'] else 'CTIO',
                        instrument = 'MaCHO' if band in ['VM', 'RM'] else '', source = source)
        journal_events()
    
        # 2015MNRAS.449..451W
        with open("../sne-external/2015MNRAS.449..451W.dat", 'r') as f:
            data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace = True)
            for r, row in enumerate(data):
                if r == 0:
                    continue
                namesplit = row[0].split('/')
                name = namesplit[-1]
                if name.startswith('SN'):
                    name = name.replace(' ', '')
                name = add_event(name)
                source = add_source(name, bibcode = '2015MNRAS.449..451W')
                add_quantity(name, 'alias', name, source)
                if len(namesplit) > 1:
                    add_quantity(name, 'alias', namesplit[0], source)
                add_quantity(name, 'claimedtype', row[1], source)
                add_photometry(name, time = row[2], band = row[4], magnitude = row[3], source = source)
        journal_events()
    
    # CCCP
    if do_task(task, 'cccp'):
        cccpbands = ['B', 'V', 'R', 'I']
        for datafile in sorted(glob.glob("../sne-external/CCCP/apj407397*.txt"), key=lambda s: s.lower()):
            with open(datafile,'r') as f:
                tsvin = csv.reader(f, delimiter='\t', skipinitialspace=True)
                for r, row in enumerate(tsvin):
                    if r == 0:
                        continue
                    elif r == 1:
                        name = 'SN' + row[0].split('SN ')[-1]
                        name = add_event(name)
                        source = add_source(name, bibcode = '2012ApJ...744...10K')
                        add_quantity(name, 'alias', name, source)
                    elif r >= 5:
                        mjd = str(Decimal(row[0]) + 53000)
                        for b, band in enumerate(cccpbands):
                            if row[2*b + 1]:
                                if not row[2*b + 2]:
                                    upplim = True
                                add_photometry(name, time = mjd, band = band, magnitude = row[2*b + 1].strip('>'),
                                    e_magnitude = row[2*b + 2], upperlimit = (not row[2*b + 2]), source = source)
    
        if tasks['cccp']['archived']:
            with open('../sne-external/CCCP/sc_cccp.html', 'r') as f:
                html = f.read()
        else:
            session = requests.Session()
            response = session.get("https://webhome.weizmann.ac.il/home/iair/sc_cccp.html")
            html = response.text
            with open('../sne-external/CCCP/sc_cccp.html', 'w') as f:
                f.write(html)
    
        soup = BeautifulSoup(html, "html5lib")
        links = soup.body.findAll("a")
        for link in tq(links):
            if 'sc_sn' in link['href']:
                name = add_event(link.text.replace(' ', ''))
                source = add_source(name, reference = 'CCCP', url = 'https://webhome.weizmann.ac.il/home/iair/sc_cccp.html')
                add_quantity(name, 'alias', name, source)
    
                if tasks['cccp']['archived']:
                    with open('../sne-external/CCCP/' + link['href'].split('/')[-1], 'r') as f:
                        html2 = f.read()
                else:
                    response2 = session.get("https://webhome.weizmann.ac.il/home/iair/" + link['href'])
                    html2 = response2.text
                    with open('../sne-external/CCCP/' + link['href'].split('/')[-1], 'w') as f:
                        f.write(html2)
    
                soup2 = BeautifulSoup(html2, "html5lib")
                links2 = soup2.body.findAll("a")
                for link2 in links2:
                    if ".txt" in link2['href'] and '_' in link2['href']:
                        band = link2['href'].split('_')[1].split('.')[0].upper()
                        if tasks['cccp']['archived']:
                            fname = '../sne-external/CCCP/' + link2['href'].split('/')[-1]
                            if not os.path.isfile(fname):
                                continue
                            with open(fname, 'r') as f:
                                html3 = f.read()
                        else:
                            response3 = session.get("https://webhome.weizmann.ac.il/home/iair/cccp/" + link2['href'])
                            if response3.status_code == 404:
                                continue
                            html3 = response3.text
                            with open('../sne-external/CCCP/' + link2['href'].split('/')[-1], 'w') as f:
                                f.write(html3)
                        table = [[str(Decimal(y.strip())).rstrip('0') for y in x.split(",")] for x in list(filter(None, html3.split("\n")))]
                        for row in table:
                            add_photometry(name, time = str(Decimal(row[0]) + 53000), band = band, magnitude = row[1], e_magnitude = row[2], source = source)
        journal_events()
    
    # Suspect catalog
    if do_task(task, 'suspect'): 
        with open('../sne-external/suspectreferences.csv','r') as f:
            tsvin = csv.reader(f, delimiter=',', skipinitialspace=True)
            suspectrefdict = {}
            for row in tsvin:
                suspectrefdict[row[0]] = row[1]
    
        for datafile in tq(sorted(glob.glob("../sne-external/SUSPECT/*.html"), key=lambda s: s.lower())):
            basename = os.path.basename(datafile)
            basesplit = basename.split('-')
            name = basesplit[1]
            name = add_event(name)
            if name.startswith('SN') and is_number(name[2:]):
                name = name + 'A'
            band = basesplit[3].split('.')[0]
            ei = int(basesplit[2])
            bandlink = 'file://' + os.path.abspath(datafile)
            bandresp = urllib.request.urlopen(bandlink)
            bandsoup = BeautifulSoup(bandresp, "html5lib")
            bandtable = bandsoup.find('table')
    
            names = bandsoup.body.findAll(text=re.compile("Name"))
            reference = ''
            for link in bandsoup.body.findAll('a'):
                if 'adsabs' in link['href']:
                    reference = str(link).replace('"', "'")
    
            bibcode = unescape(suspectrefdict[reference])
            source = add_source(name, bibcode = bibcode)
    
            secondaryreference = "SUSPECT"
            secondaryrefurl = "https://www.nhn.ou.edu/~suspect/"
            secondarysource = add_source(name, reference = secondaryreference, url = secondaryrefurl, secondary = True)
            add_quantity(name, 'alias', name, secondarysource)
    
            if ei == 1:
                year = re.findall(r'\d+', name)[0]
                add_quantity(name, 'discoverdate', year, secondarysource)
                add_quantity(name, 'host', names[1].split(':')[1].strip(), secondarysource)
    
                redshifts = bandsoup.body.findAll(text=re.compile("Redshift"))
                if redshifts:
                    add_quantity(name, 'redshift', redshifts[0].split(':')[1].strip(), secondarysource, kind = 'heliocentric')
                hvels = bandsoup.body.findAll(text=re.compile("Heliocentric Velocity"))
                #if hvels:
                #    add_quantity(name, 'velocity', hvels[0].split(':')[1].strip().split(' ')[0],
                #        secondarysource, kind = 'heliocentric')
                types = bandsoup.body.findAll(text=re.compile("Type"))
    
                add_quantity(name, 'claimedtype', types[0].split(':')[1].strip().split(' ')[0], secondarysource)
    
            for r, row in enumerate(bandtable.findAll('tr')):
                if r == 0:
                    continue
                col = row.findAll('td')
                mjd = str(jd_to_mjd(Decimal(col[0].contents[0])))
                mag = col[3].contents[0]
                if mag.isspace():
                    mag = ''
                else:
                    mag = str(mag)
                e_magnitude = col[4].contents[0]
                if e_magnitude.isspace():
                    e_magnitude = ''
                else:
                    e_magnitude = str(e_magnitude)
                add_photometry(name, time = mjd, band = band, magnitude = mag, e_magnitude = e_magnitude, source = secondarysource + ',' + source)
        journal_events()
    
    # CfA data
    if do_task(task, 'cfa'): 
        for fname in tq(sorted(glob.glob("../sne-external/cfa-input/*.dat"), key=lambda s: s.lower())):
            f = open(fname,'r')
            tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)
            csv_data = []
            for r, row in enumerate(tsvin):
                new = []
                for item in row:
                    new.extend(item.split("\t"))
                csv_data.append(new)
    
            for r, row in enumerate(csv_data):
                for c, col in enumerate(row):
                    csv_data[r][c] = col.strip()
                csv_data[r] = [_f for _f in csv_data[r] if _f]
    
            eventname = os.path.basename(os.path.splitext(fname)[0])
    
            eventparts = eventname.split('_')
    
            name = snname(eventparts[0])
            name = add_event(name)
            secondaryname = 'CfA Supernova Archive'
            secondaryurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
            secondarysource = add_source(name, reference = secondaryname, url = secondaryurl, secondary = True)
            add_quantity(name, 'alias', name, secondarysource)
    
            year = re.findall(r'\d+', name)[0]
            add_quantity(name, 'discoverdate', year, secondarysource)
    
            eventbands = list(eventparts[1])
    
            tu = 'MJD'
            jdoffset = Decimal(0.)
            for rc, row in enumerate(csv_data):
                if len(row) > 0 and row[0][0] == "#":
                    if len(row[0]) > 2 and row[0][:3] == "#JD":
                        tu = 'JD'
                        rowparts = row[0].split('-')
                        jdoffset = Decimal(rowparts[1])
                    elif len(row[0]) > 6 and row[0][:7] == "#Julian":
                        tu = 'JD'
                        jdoffset = Decimal(0.)
                    elif len(row) > 1 and row[1].lower() == "photometry":
                        for ci, col in enumerate(row[2:]):
                            if col[0] == "(":
                                refstr = ' '.join(row[2+ci:])
                                refstr = refstr.replace('(','').replace(')','')
                                bibcode = unescape(refstr)
                                source = add_source(name, bibcode = bibcode)
    
                    elif len(row) > 1 and row[1] == "HJD":
                        tu = "HJD"
    
                    continue
                elif len(row) > 0:
                    mjd = row[0]
                    for v, val in enumerate(row):
                        if v == 0:
                            if tu == 'JD':
                                mjd = str(jd_to_mjd(Decimal(val) + jdoffset))
                                tuout = 'MJD'
                            elif tu == 'HJD':
                                mjd = str(jd_to_mjd(Decimal(val)))
                                tuout = 'MJD'
                            else:
                                mjd = val
                                tuout = tu
                        elif v % 2 != 0:
                            if float(row[v]) < 90.0:
                                add_photometry(name, u_time = tuout, time = mjd, band = eventbands[(v-1)//2], magnitude = row[v], e_magnitude = row[v+1], source = secondarysource + ',' + source)
            f.close()
    
        # Hicken 2012
        f = open("../sne-external/hicken-2012-standard.dat", 'r')
        tsvin = csv.reader(f, delimiter='|', skipinitialspace=True)
        for r, row in enumerate(tq(tsvin)):
            if r <= 47:
                continue
    
            if row[0][:2] != 'sn':
                name = 'SN' + row[0].strip()
            else:
                name = row[0].strip()
    
            name = add_event(name)
    
            source = add_source(name, bibcode = '2012ApJS..200...12H')
            add_quantity(name, 'alias', name, source)
            add_quantity(name, 'claimedtype', 'Ia', source)
            add_photometry(name, u_time = 'MJD', time = row[2].strip(), band = row[1].strip(),
                magnitude = row[6].strip(), e_magnitude = row[7].strip(), source = source)
        
        # Bianco 2014
        tsvin = open("../sne-external/bianco-2014-standard.dat", 'r')
        tsvin = csv.reader(tsvin, delimiter=' ', skipinitialspace=True)
        for row in tq(tsvin):
            name = 'SN' + row[0]
            name = add_event(name)
    
            source = add_source(name, bibcode = '2014ApJS..213...19B')
            add_quantity(name, 'alias', name, source)
            add_photometry(name, u_time = 'MJD', time = row[2], band = row[1], magnitude = row[3],
                e_magnitude = row[4], telescope = row[5], system = "Standard", source = source)
        f.close()
        journal_events()
    
    # Now import the UCB SNDB
    if do_task(task, 'ucb'): 
        for fname in tq(sorted(glob.glob("../sne-external/SNDB/*.dat"), key=lambda s: s.lower())):
            f = open(fname,'r')
            tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)
    
            eventname = os.path.basename(os.path.splitext(fname)[0])
    
            eventparts = eventname.split('.')
    
            name = snname(eventparts[0])
            name = add_event(name)
    
            reference = "UCB Filippenko Group's Supernova Database (SNDB)"
            refurl = "http://heracles.astro.berkeley.edu/sndb/info"
            refbib = "2012MNRAS.425.1789S"
            source = add_source(name, reference = reference, url = refurl, bibcode = refbib, secondary = True)
            add_quantity(name, 'alias', name, source)
    
            year = re.findall(r'\d+', name)[0]
            add_quantity(name, 'discoverdate', year, source)
    
            for r, row in enumerate(tsvin):
                if len(row) > 0 and row[0] == "#":
                    continue
                mjd = row[0]
                magnitude = row[1]
                if magnitude and float(magnitude) > 99.0:
                    continue
                e_magnitude = row[2]
                band = row[4]
                telescope = row[5]
                add_photometry(name, time = mjd, telescope = telescope, band = band, magnitude = magnitude,
                    e_magnitude = e_magnitude, source = source)
            f.close()
        journal_events()
        
    # Import SDSS
    if do_task(task, 'sdss'): 
        with open('../sne-external/SDSS/2010ApJ...708..661D.txt', 'r') as f:
            bibcodes2010 = f.read().split("\n")
        sdssbands = ['u', 'g', 'r', 'i', 'z']
        for fname in tq(sorted(glob.glob("../sne-external/SDSS/*.sum"), key=lambda s: s.lower())):
            f = open(fname,'r')
            tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)
    
            basename = os.path.basename(fname)
            if basename in bibcodes2010:
                bibcode = '2010ApJ...708..661D'
            else:
                bibcode = '2008AJ....136.2306H'
    
            for r, row in enumerate(tsvin):
                if r == 0:
                    if row[5] == "RA:":
                        name = "SDSS-II " + row[3]
                    else:
                        name = "SN" + row[5]
                    name = add_event(name)
                    source = add_source(name, bibcode = bibcode)
                    add_quantity(name, 'alias', name, source)
                    add_quantity(name, 'alias', "SDSS-II " + row[3], source)
    
                    if row[5] != "RA:":
                        year = re.findall(r'\d+', name)[0]
                        add_quantity(name, 'discoverdate', year, source)
    
                    add_quantity(name, 'ra', row[-4], source, unit = 'floatdegrees')
                    add_quantity(name, 'dec', row[-2], source, unit = 'floatdegrees')
                if r == 1:
                    error = row[4] if float(row[4]) >= 0.0 else ''
                    add_quantity(name, 'redshift', row[2], source, error = error, kind = 'heliocentric')
                if r >= 19:
                    # Skip bad measurements
                    if int(row[0]) > 1024:
                        continue
    
                    mjd = row[1]
                    band = sdssbands[int(row[2])]
                    magnitude = row[3]
                    e_magnitude = row[4]
                    telescope = "SDSS"
                    add_photometry(name, time = mjd, telescope = telescope, band = band, magnitude = magnitude,
                        e_magnitude = e_magnitude, source = source, system = "SDSS")
            f.close()
        journal_events()
    
    #Import GAIA
    if do_task(task, 'gaia'): 
        fname = '../sne-external/GAIA/alerts.csv'
        csvtxt = load_cached_url('http://gsaweb.ast.cam.ac.uk/alerts/alerts.csv', fname)
        if not csvtxt:
            continue
        tsvin = csv.reader(csvtxt.splitlines(), delimiter=',', skipinitialspace=True)
        reference = "Gaia Photometric Science Alerts"
        refurl = "http://gsaweb.ast.cam.ac.uk/alerts/alertsindex"
        for ri, row in enumerate(tq(tsvin)):
            if ri == 0 or not row:
                continue
            name = add_event(row[0])
            source = add_source(name, reference = reference, url = refurl)
            add_quantity(name, 'alias', name, source)
            year = '20' + re.findall(r'\d+', row[0])[0]
            add_quantity(name, 'discoverdate', year, source)
            add_quantity(name, 'ra', row[2], source, unit = 'floatdegrees')
            add_quantity(name, 'dec', row[3], source, unit = 'floatdegrees')
            if row[7] and row[7] != 'unknown':
                add_quantity(name, 'claimedtype', row[7].replace('SNe', '').replace('SN', '').strip(), source)
            elif (True in [x in row[9].upper() for x in ['SN CANDIATE', 'CANDIDATE SN', 'HOSTLESS SN']]):
                add_quantity(name, 'claimedtype', 'Candidate', source)
    
            if 'aka' in row[9].replace('gakaxy','galaxy').lower() and 'AKARI' not in row[9]:
                commentsplit = (row[9].replace('_', ' ').replace('MLS ', 'MLS').replace('CSS ', 'CSS').
                    replace('SN iPTF', 'iPTF').replace('SN ', 'SN').replace('AT ', 'AT').split())
                for csi, cs in enumerate(commentsplit):
                    if 'aka' in cs.lower() and csi < len(commentsplit) - 1:
                        alias = commentsplit[csi+1].strip('(),:.').replace('PSNJ', 'PSN J')
                        if alias[:6] == 'ASASSN' and alias[6] != '-':
                            alias = 'ASASSN-' + alias[6:]
                        add_quantity(name, 'alias', alias, source)
                        break
    
            fname = '../sne-external/GAIA/' + row[0] + '.csv'
            if not args.fullrefresh and tasks['gaia']['archived'] and os.path.isfile(fname):
                with open(fname, 'r') as f:
                    csvtxt = f.read()
            else:
                response = urllib.request.urlopen("http://gsaweb.ast.cam.ac.uk/alerts/alert/" + row[0] + "/lightcurve.csv")
                with open(fname, 'w') as f:
                    csvtxt = response.read().decode('utf-8')
                    f.write(csvtxt)
    
            tsvin2 = csv.reader(csvtxt.splitlines())
            for ri2, row2 in enumerate(tsvin2):
                if ri2 <= 1 or not row2:
                    continue
                mjd = str(jd_to_mjd(Decimal(row2[1].strip())))
                magnitude = row2[2].strip()
                if magnitude == 'null':
                    continue
                e_magnitude = 0.
                telescope = 'GAIA'
                band = 'G'
                add_photometry(name, time = mjd, telescope = telescope, band = band, magnitude = magnitude, e_magnitude = e_magnitude, source = source)
            if args.update:
                journal_events()
        journal_events()
    
    # Import CSP
    if do_task(task, 'csp'): 
        cspbands = ['u', 'B', 'V', 'g', 'r', 'i', 'Y', 'J', 'H', 'K']
        for fname in ta(sorted(glob.glob("../sne-external/CSP/*.dat"), key=lambda s: s.lower())):
            f = open(fname,'r')
            tsvin = csv.reader(f, delimiter='\t', skipinitialspace=True)
    
            eventname = os.path.basename(os.path.splitext(fname)[0])
    
            eventparts = eventname.split('opt+')
    
            name = snname(eventparts[0])
            name = add_event(name)
    
            reference = "Carnegie Supernova Project"
            refurl = "http://csp.obs.carnegiescience.edu/data"
            source = add_source(name, reference = reference, url = refurl)
            add_quantity(name, 'alias', name, source)
    
            year = re.findall(r'\d+', name)[0]
            add_quantity(name, 'discoverdate', year, source)
    
            for r, row in enumerate(tsvin):
                if len(row) > 0 and row[0][0] == "#":
                    if r == 2:
                        add_quantity(name, 'redshift', row[0].split(' ')[-1], source, kind = 'cmb')
                        add_quantity(name, 'ra', row[1].split(' ')[-1], source)
                        add_quantity(name, 'dec', row[2].split(' ')[-1], source)
                    continue
                for v, val in enumerate(row):
                    if v == 0:
                        mjd = val
                    elif v % 2 != 0:
                        if float(row[v]) < 90.0:
                            add_photometry(name, time = mjd, observatory = 'LCO', band = cspbands[(v-1)//2],
                                system = 'CSP', magnitude = row[v], e_magnitude = row[v+1], source = source)
            f.close()
        journal_events()
    
    # Import ITEP
    if do_task(task, 'itep'): 
        itepphotometryerrors = ['SN1995N']
        itepbadsources = ['2004ApJ...602..571B']
    
        needsbib = []
        with open("../sne-external/itep-refs.txt",'r') as f:
            refrep = f.read().splitlines()
        refrepf = dict(list(zip(refrep[1::2], refrep[::2])))
        f = open("../sne-external/itep-lc-cat-28dec2015.txt",'r')
        tsvin = csv.reader(f, delimiter='|', skipinitialspace=True)
        curname = ''
        for r, row in enumerate(tq(tsvin)):
            if r <= 1 or len(row) < 7:
                continue
            name = 'SN' + row[0].strip()
            mjd = str(jd_to_mjd(Decimal(row[1].strip())))
            band = row[2].strip()
            magnitude = row[3].strip()
            e_magnitude = row[4].strip()
            reference = row[6].strip().strip(',')
    
            if curname != name:
                curname = name
                name = add_event(name)
    
                secondaryreference = "Sternberg Astronomical Institute Supernova Light Curve Catalogue"
                secondaryrefurl = "http://dau.itep.ru/sn/node/72"
                secondarysource = add_source(name, reference = secondaryreference, url = secondaryrefurl, secondary = True)
                add_quantity(name, 'alias', name, secondarysource)
    
                year = re.findall(r'\d+', name)[0]
                add_quantity(name, 'discoverdate', year, secondarysource)
            if reference in refrepf:
                bibcode = unescape(refrepf[reference])
                source = add_source(name, bibcode = bibcode)
            else:
                needsbib.append(reference)
                source = add_source(name, reference = reference) if reference else ''
    
            if name not in itepphotometryerrors and bibcode not in itepbadsources:
                add_photometry(name, time = mjd, band = band, magnitude = magnitude, e_magnitude = e_magnitude, source = secondarysource + ',' + source)
        f.close()
        
        # Write out references that could use a bibcode
        needsbib = list(OrderedDict.fromkeys(needsbib))
        with open('../itep-needsbib.txt', 'w') as f:
            f.writelines(["%s\n" % i for i in needsbib])
        journal_events()
    
    # Now import the Asiago catalog
    if do_task(task, 'asiago'): 
        asiagopositionerrors = ['SN2011in', 'SN2012ac', 'SN2012at', 'SN2013bz']
        #response = urllib.request.urlopen('http://graspa.oapd.inaf.it/cgi-bin/sncat.php')
        path = os.path.abspath('../sne-external/asiago-cat.php')
        response = urllib.request.urlopen('file://' + path)
        html = response.read().decode('utf-8')
        html = html.replace("\r", "")
    
        soup = BeautifulSoup(html, "html5lib")
        table = soup.find("table")
    
        records = []
        for r, row in enumerate(table.findAll('tr')):
            if r == 0:
                continue
    
            col = row.findAll('td')
            records.append([utf8(x.renderContents()) for x in col])
    
        for record in tq(records):
            if len(record) > 1 and record[1] != '':
                name = snname("SN" + record[1]).strip('?')
                name = add_event(name)
    
                reference = 'Asiago Supernova Catalogue'
                refurl = 'http://graspa.oapd.inaf.it/cgi-bin/sncat.php'
                refbib = '1989A&AS...81..421B'
                source = add_source(name, reference = reference, url = refurl, bibcode = refbib, secondary = True)
                add_quantity(name, 'alias', name, source)
    
                year = re.findall(r'\d+', name)[0]
                add_quantity(name, 'discoverdate', year, source)
    
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
                    monthstr = ''.join(re.findall("[a-zA-Z]+", monthday))
                    monthstr = str(list(calendar.month_abbr).index(monthstr))
                    datestring = datestring + '/' + monthstr
    
                    dayarr = re.findall(r'\d+', monthday)
                    if dayarr:
                        daystr = dayarr[0]
                        datestring = datestring + '/' + daystr
    
                add_quantity(name, datekey + 'date', datestring, source)
    
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
                    add_quantity(name, 'host', hostname, source)
                if (claimedtype != ''):
                    add_quantity(name, 'claimedtype', claimedtype, source)
                if (redshift != ''):
                    add_quantity(name, 'redshift', redshift, source, kind = 'host')
                if (velocity != ''):
                    add_quantity(name, 'velocity', velocity, source, kind = 'host')
                if (hostra != ''):
                    add_quantity(name, 'hostra', hostra, source, unit = 'nospace')
                if (hostdec != ''):
                    add_quantity(name, 'hostdec', hostdec, source, unit = 'nospace')
                if (ra != '' and name not in asiagopositionerrors):
                    add_quantity(name, 'ra', ra, source, unit = 'nospace')
                if (dec != '' and name not in asiagopositionerrors):
                    add_quantity(name, 'dec', dec, source, unit = 'nospace')
                if (discoverer != ''):
                    add_quantity(name, 'discoverer', discoverer, source)
        journal_events()
    
    if do_task(task, 'lennarz'): 
        Vizier.ROW_LIMIT = -1
        result = Vizier.get_catalogs("J/A+A/538/A120/usc")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
    
        bibcode = "2012A&A...538A.120L"
        for row in tq(table):
            row = convert_aq_output(row)
            name = 'SN' + row['SN']
            name = add_event(name)
    
            source = add_source(name, bibcode = bibcode)
            add_quantity(name, 'alias', name, source)
    
            if row['RAJ2000']:
                add_quantity(name, 'ra', row['RAJ2000'], source)
            if row['DEJ2000']:
                add_quantity(name, 'dec', row['DEJ2000'], source)
            if row['RAG']:
                add_quantity(name, 'hostra', row['RAG'], source)
            if row['DEG']:
                add_quantity(name, 'hostdec', row['DEG'], source)
            if row['Gal']:
                add_quantity(name, 'host', row['Gal'], source)
            if row['Type']:
                claimedtypes = row['Type'].split('|')
                for claimedtype in claimedtypes:
                    add_quantity(name, 'claimedtype', claimedtype.strip(' -'), source)
            if row['z']:
                if name not in ['SN1985D', 'SN2004cq']:
                    add_quantity(name, 'redshift', row['z'], source, kind = 'host')
            if row['Dist']:
                if row['e_Dist']:
                    add_quantity(name, 'lumdist', row['Dist'], source, error = row['e_Dist'], kind = 'host')
                else:
                    add_quantity(name, 'lumdist', row['Dist'], source, kind = 'host')
    
            if row['Ddate']:
                datestring = row['Ddate'].replace('-', '/')
    
                add_quantity(name, 'discoverdate', datestring, source)
    
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
                        add_photometry(name, time = mjd, band = row['Dband'], magnitude = row['Dmag'], source = source)
            if row['Mdate']:
                datestring = row['Mdate'].replace('-', '/')
    
                add_quantity(name, 'maxdate', datestring, source)
    
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
                        add_photometry(name, time = mjd, band = row['Mband'], magnitude = row['Mmag'], source = source)
        f.close()
        journal_events()
    
    if do_task(task, 'tns'):
        session = requests.Session()
        csvtxt = load_cached_url("https://wis-tns.weizmann.ac.il/search?&isTNS_AT=yes&num_page=1&format=html&sort=desc&order=id&format=csv&page=0",
            "../sne-external/TNS/index.csv")
        if not csvtxt:
            continue
        maxid = csvtxt.splitlines()[1].split(",")[0].strip('"')
        maxpages = ceil((int(maxid) - 10000)/1000.)
    
        for page in tq(range(maxpages)):
            fname = '../sne-external/TNS/page-' + str(page).zfill(2) + '.csv'
            if tasks['tns']['archived'] and os.path.isfile(fname) and page != maxpages:
                with open(fname, 'r') as f:
                    csvtxt = f.read()
            else:
                with open(fname, 'w') as f:
                    session = requests.Session()
                    response = session.get("https://wis-tns.weizmann.ac.il/search?&isTNS_AT=yes&num_page=1000&format=html&edit[type]=&edit[objname]=&edit[id]=&sort=asc&order=id&display[redshift]=1&display[hostname]=1&display[host_redshift]=1&display[source_group_name]=1&display[programs_name]=1&display[internal_name]=1&display[isTNS_AT]=1&display[public]=1&display[end_pop_period]=0&display[spectra_count]=1&display[discoverymag]=1&display[discmagfilter]=1&display[discoverydate]=1&display[discoverer]=1&display[sources]=1&display[bibcode]=1&format=csv&page=" + str(page))
                    csvtxt = response.text
                    f.write(csvtxt)
    
            tsvin = csv.reader(csvtxt.splitlines(), delimiter=',')
            for ri, row in enumerate(tq(tsvin, leave = False)):
                if ri == 0:
                    continue
                if row[4] and 'SN' not in row[4]:
                    continue
                name = row[1].replace(' ', '')
                name = add_event(name)
                source = add_source(name, reference = 'Transient Name Server', url = 'https://wis-tns.weizmann.ac.il')
                add_quantity(name, 'alias', name, source)
                if row[2]:
                    add_quantity(name, 'ra', row[2], source)
                if row[3]:
                    add_quantity(name, 'dec', row[3], source)
                if row[4]:
                    add_quantity(name, 'claimedtype', row[4].replace('SN', '').strip(), source)
                if row[5]:
                    add_quantity(name, 'redshift', row[5], source, kind = 'spectroscopic')
                if row[6]:
                    add_quantity(name, 'host', row[6], source)
                if row[7]:
                    add_quantity(name, 'redshift', row[7], source, kind = 'host')
                if row[8]:
                    add_quantity(name, 'discoverer', row[8], source)
                if row[9]:
                    observers = row[9].split(',')
                    for observer in observers:
                        add_quantity(name, 'observer', observer.strip(), source)
                if row[10]:
                    add_quantity(name, 'alias', row[10], source)
                if row[8] and row[14] and row[15] and row[16]:
                    survey = row[8]
                    magnitude = row[14]
                    band = row[15].split('-')[0]
                    mjd = astrotime(row[16]).mjd
                    add_photometry(name, time = mjd, magnitude = magnitude, band = band, survey = survey, source = source)
                if args.update:
                    journal_events()
        journal_events()
    
    if do_task(task, 'rochester'): 
        rochesterpaths = ['http://www.rochesterastronomy.org/snimages/snredshiftall.html', 'http://www.rochesterastronomy.org/sn2016/snredshift.html']
        rochesterupdate = [False, True]
    
        # These are known to be in error on the Rochester page, so ignore them.
        rochesterredshifterrors = ['LSQ12bgl','LSQ12axx']
        rochesterphotometryerrors = ['SNF20080514-002','SN1998ev']
        rochestertypeerrors = ['SN1054A']
        rochestercoordinateerrors = ['MASTER OT J095321.02+202721.2', 'SN1996D', 'SN1998ew', 'SN2003an']
    
        for p, path in enumerate(tq(rochesterpaths)):
            if args.update and not rochesterupdate[p]:
                continue
    
            filepath = '../sne-external/rochester/' + os.path.basename(path)
            html = load_cached_url(path, filepath)
            if not html:
                continue
    
            soup = BeautifulSoup(html, "html5lib")
            rows = soup.findAll('tr')
            secondaryreference = "Latest Supernovae"
            secondaryrefurl = "http://www.rochesterastronomy.org/snimages/snredshiftall.html"
            for r, row in enumerate(tq(rows)):
                if r == 0:
                    continue
                cols = row.findAll('td')
                if not len(cols):
                    continue
    
                name = ''
                if cols[14].contents:
                    aka = str(cols[14].contents[0]).strip()
                    if is_number(aka.strip('?')):
                        aka = 'SN' + aka.strip('?') + 'A'
                        name = add_event(aka)
                    elif len(aka) >= 4 and is_number(aka[:4]):
                        aka = 'SN' + aka
                        name = add_event(aka)
    
                sn = re.sub('<[^<]+?>', '', str(cols[0].contents[0])).strip()
                if is_number(sn.strip('?')):
                    sn = 'SN' + sn.strip('?') + 'A'
                elif len(sn) >= 4 and is_number(sn[:4]):
                    sn = 'SN' + sn
                if not name:
                    if not sn:
                        continue
                    if sn[:8] == 'MASTER J':
                        sn = sn.replace('MASTER J', 'MASTER OT J').replace('SNHunt', 'SNhunt')
                    name = add_event(sn)
    
                reference = cols[12].findAll('a')[0].contents[0].strip()
                refurl = cols[12].findAll('a')[0]['href'].strip()
                source = add_source(name, reference = reference, url = refurl)
                secondarysource = add_source(name, reference = secondaryreference, url = secondaryrefurl, secondary = True)
                add_quantity(name, 'alias', name, secondarysource)
                sources = uniq_cdl(list(filter(None, [source, secondarysource])))
    
                add_quantity(name, 'alias', sn, source)
    
                if cols[14].contents:
                    if aka == 'SNR G1.9+0.3':
                        aka = 'G001.9+00.3'
                    if aka[:4] == 'PS1 ':
                        aka = 'PS1-' + aka[4:]
                    if aka[:8] == 'MASTER J':
                        aka = aka.replace('MASTER J', 'MASTER OT J').replace('SNHunt', 'SNhunt')
                    add_quantity(name, 'alias', aka, source)
    
                if str(cols[1].contents[0]).strip() != 'unk' and name not in rochestertypeerrors:
                    add_quantity(name, 'claimedtype', str(cols[1].contents[0]).strip(' :,'), sources)
                if str(cols[2].contents[0]).strip() != 'anonymous':
                    add_quantity(name, 'host', str(cols[2].contents[0]).strip(), sources)
                if name not in rochestercoordinateerrors:
                    add_quantity(name, 'ra', str(cols[3].contents[0]).strip(), sources)
                    add_quantity(name, 'dec', str(cols[4].contents[0]).strip(), sources)
                if str(cols[6].contents[0]).strip() not in ['2440587', '2440587.292']:
                    astrot = astrotime(float(str(cols[6].contents[0]).strip()), format='jd').datetime
                    add_quantity(name, 'discoverdate', make_date_string(astrot.year, astrot.month, astrot.day), sources)
                if str(cols[7].contents[0]).strip() not in ['2440587', '2440587.292']:
                    astrot = astrotime(float(str(cols[7].contents[0]).strip()), format='jd')
                    if (float(str(cols[8].contents[0]).strip()) <= 90.0 and name not in rochesterphotometryerrors and
                        not any('GRB' in x for x in get_aliases(name))):
                        add_photometry(name, time = str(astrot.mjd), magnitude = str(cols[8].contents[0]).strip(), source = sources)
                if cols[11].contents[0] != 'n/a' and name not in rochesterredshifterrors:
                    add_quantity(name, 'redshift', str(cols[11].contents[0]).strip(), sources)
                add_quantity(name, 'discoverer', str(cols[13].contents[0]).strip(), sources)
                if args.update:
                    journal_events()
    
        if not args.update:
            vsnetfiles = ["latestsne.dat"]
            for vsnetfile in vsnetfiles:
                f = open("../sne-external/" + vsnetfile,'r',encoding='latin1')
                tsvin = csv.reader(f, delimiter=' ', skipinitialspace=True)
                for r, row in enumerate(tsvin):
                    if not row or row[0][:4] in ['http', 'www.'] or len(row) < 3:
                        continue
                    name = row[0].strip()
                    if name[:4].isdigit():
                        name = 'SN' + name
                    if name.startswith('PSNJ'):
                        name = 'PSN J' + name[4:]
                    if name.startswith('MASTEROTJ'):
                        name = name.replace('MASTEROTJ', 'MASTER OT J')
                    name = add_event(name)
                    secondarysource = add_source(name, reference = secondaryreference, url = secondaryrefurl, secondary = True)
                    add_quantity(name, 'alias', name, secondarysource)
    
                    if not is_number(row[1]):
                        continue
                    year = row[1][:4]
                    month = row[1][4:6]
                    day = row[1][6:]
                    if '.' not in day:
                        day = day[:2] + '.' + day[2:]
                    mjd = astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day))
                    magnitude = row[2].rstrip(ascii_letters)
                    if not is_number(magnitude):
                        continue
                    if magnitude.isdigit():
                        if int(magnitude) > 100:
                            magnitude = magnitude[:2] + '.' + magnitude[2:]
    
                    if float(str(cols[8].contents[0]).strip()) >= 90.0:
                        continue
    
                    if len(row) >= 4:
                        if is_number(row[3]):
                            e_magnitude = row[3]
                            refind = 4
                        else:
                            e_magnitude = ''
                            refind = 3
    
                        if refind >= len(row):
                            sources = secondarysource
                        else:
                            reference = ' '.join(row[refind:])
                            source = add_source(name, reference = reference)
                            add_quantity(name, 'alias', name, secondarysource)
                            sources = uniq_cdl([source,secondarysource])
                    else:
                        sources = secondarysource
    
                    band = row[2].lstrip('1234567890.')
    
                    add_photometry(name, time = mjd, band = band, magnitude = magnitude, e_magnitude = e_magnitude, source = sources)
                f.close()
        journal_events()
    
    if do_task(task, 'ogle'): 
        basenames = ['transients', 'transients/2014b', 'transients/2014', 'transients/2013', 'transients/2012']
        oglenames = []
        ogleupdate = [True, False, False, False, False]
        oglepositionerrors = ['OGLE-2015-SN-078']
        for b, bn in enumerate(tq(basenames)):
            if args.update and not ogleupdate[b]:
                continue
    
            filepath = '../sne-external/OGLE/' + bn + '-' + 'transients.html'
            htmltxt = load_cached_url('http://ogle.astrouw.edu.pl/ogle4/' + bn + '/transients.html', filepath)
            if not htmltxt:
                continue
    
            soup = BeautifulSoup(htmltxt, "html5lib")
            links = soup.findAll('a')
            breaks = soup.findAll('br')
            datalinks = []
            datafnames = []
            for a in links:
                if a.has_attr('href'):
                    if '.dat' in a['href']:
                        datalinks.append('http://ogle.astrouw.edu.pl/ogle4/' + bn + '/' + a['href'])
                        datafnames.append(bn + '-' + a['href'].replace('/', '-'))
    
            ec = -1
            reference = 'OGLE-IV Transient Detection System'
            refurl = 'http://ogle.astrouw.edu.pl/ogle4/transients/transients.html'
            for br in tq(breaks):
                sibling = br.nextSibling
                if 'Ra,Dec=' in sibling:
                    line = sibling.replace("\n", '').split('Ra,Dec=')
                    name = line[0].strip()
                    ec += 1
    
                    if 'NOVA' in name or 'dupl' in name:
                        continue
    
                    if name in oglenames:
                        continue
                    oglenames.append(name)
    
                    name = add_event(name)
    
                    mySibling = sibling.nextSibling
                    atelref = ''
                    claimedtype = ''
                    while 'Ra,Dec=' not in mySibling:
                        if isinstance(mySibling, NavigableString):
                            if 'Phot.class=' in str(mySibling):
                                claimedtype = re.sub(r'\([^)]*\)', '', str(mySibling).split('=')[-1]).replace('SN','').strip()
                        if isinstance(mySibling, Tag):
                            atela = mySibling
                            if atela and atela.has_attr('href') and 'astronomerstelegram' in atela['href']:
                                atelref = atela.contents[0].strip()
                                atelurl = atela['href']
                        mySibling = mySibling.nextSibling
                        if mySibling is None:
                            break
    
                    nextSibling = sibling.nextSibling
                    if isinstance(nextSibling, Tag) and nextSibling.has_attr('alt') and nextSibling.contents[0].strip() != 'NED':
                        radec = nextSibling.contents[0].strip().split()
                    else:
                        radec = line[-1].split()
                    ra = radec[0]
                    dec = radec[1]
    
                    fname = '../sne-external/OGLE/' + datafnames[ec]
                    if not args.fullrefresh and tasks['ogle']['archived'] and os.path.isfile(fname):
                        with open(fname, 'r') as f:
                            csvtxt = f.read()
                    else:
                        response = urllib.request.urlopen(datalinks[ec])
                        with open(fname, 'w') as f:
                            csvtxt = response.read().decode('utf-8')
                            f.write(csvtxt)
    
                    lcdat = csvtxt.splitlines()
                    sources = [add_source(name, reference = reference, url = refurl)]
                    add_quantity(name, 'alias', name, sources[0])
                    if atelref and atelref != 'ATel#----':
                        sources.append(add_source(name, reference = atelref, url = atelurl))
                    sources = uniq_cdl(sources)
    
                    if name.startswith('OGLE'):
                        if name[4] == '-':
                            if is_number(name[5:9]):
                                add_quantity(name, 'discoverdate', name[5:9], sources)
                        else:
                            if is_number(name[4:6]):
                                add_quantity(name, 'discoverdate', '20' + name[4:6], sources)
    
                    if name not in oglepositionerrors:
                        add_quantity(name, 'ra', ra, sources)
                        add_quantity(name, 'dec', dec, sources)
                    if claimedtype and claimedtype != '-':
                        add_quantity(name, 'claimedtype', claimedtype, sources)
                    elif 'SN' not in name and 'claimedtype' not in events[name]:
                        add_quantity(name, 'claimedtype', 'Candidate', sources)
                    for row in lcdat:
                        row = row.split()
                        mjd = str(jd_to_mjd(Decimal(row[0])))
                        magnitude = row[1]
                        if float(magnitude) > 90.0:
                            continue
                        e_magnitude = row[2]
                        upperlimit = False
                        if e_magnitude == '-1' or float(e_magnitude) > 10.0:
                            e_magnitude = ''
                            upperlimit = True
                        add_photometry(name, time = mjd, band = 'I', magnitude = magnitude, e_magnitude = e_magnitude, source = sources, upperlimit = upperlimit)
                    if args.update:
                        journal_events()
            journal_events()
    
    if do_task(task, 'snls'): 
        with open("../sne-external/SNLS-ugriz.dat", 'r') as f:
            data = csv.reader(f, delimiter=' ', quotechar='"', skipinitialspace = True)
            for row in data:
                flux = row[3]
                err = row[4]
                # Being extra strict here with the flux constraint, see note below.
                if float(flux) < 3.0*float(err):
                    continue
                name = 'SNLS-' + row[0]
                name = add_event(name)
                source = add_source(name, bibcode = '2010A&A...523A...7G')
                add_quantity(name, 'alias', name, source)
                band = row[1]
                mjd = row[2]
                sig = get_sig_digits(flux.split('E')[0])+1
                # Conversion comes from SNLS-Readme
                # NOTE: Datafiles available for download suggest different zeropoints than 30, need to inquire.
                magnitude = pretty_num(30.0-2.5*log10(float(flux)), sig = sig)
                e_magnitude = pretty_num(2.5*log10(1.0 + float(err)/float(flux)), sig = sig)
                #e_magnitude = pretty_num(2.5*(log10(float(flux) + float(err)) - log10(float(flux))), sig = sig)
                add_photometry(name, time = mjd, band = band, magnitude = magnitude, e_magnitude = e_magnitude, source = source)
        journal_events()
    
    if do_task(task, 'psthreepi'):
        fname = '../sne-external/3pi/page01.html'
        html = load_cached_url("http://psweb.mp.qub.ac.uk/ps1threepi/psdb/public/?page=1&sort=followup_flag_date", fname)
        if not html:
            continue

        bs = BeautifulSoup(html, "html5lib")
        div = bs.find('div', {"class":"pagination"})
        offline = False
        if not div:
            offline = True
        else:
            links = div.findAll('a')
            if not links:
                offline = True
    
        if offline:
            warnings.warn("Pan-STARRS 3pi offline, using local files only.")
            with open(fname, 'r') as f:
                html = f.read()
            bs = BeautifulSoup(html, "html5lib")
            div = bs.find('div', {"class":"pagination"})
            links = div.findAll('a')
    
        numpages = int(links[-2].contents[0])
        oldnumpages = len(glob.glob('../sne-external/3pi/page*'))
        for page in tq(range(1,numpages)):
            fname = '../sne-external/3pi/page' + str(page).zfill(2) + '.html'
            if not args.fullrefresh and tasks['psthreepi']['archived'] and os.path.isfile(fname) and page < oldnumpages:
                with open(fname, 'r') as f:
                    html = f.read()
            elif not offline:
                response = urllib.request.urlopen("http://psweb.mp.qub.ac.uk/ps1threepi/psdb/public/?page=" + str(page) + "&sort=followup_flag_date")
                with open(fname, 'w') as f:
                    html = response.read().decode('utf-8')
                    f.write(html)
            else:
                continue
    
            bs = BeautifulSoup(html, "html5lib")
            trs = bs.findAll('tr')
            for tr in tq(trs):
                tds = tr.findAll('td')
                if not tds:
                    continue
                refs = []
                aliases = []
                ttype = ''
                ctype = ''
                for tdi, td in enumerate(tds):
                    if tdi == 0:
                        psname = td.contents[0]
                        pslink = psname['href']
                        psname = psname.text
                    elif tdi == 1:
                        ra = td.contents[0]
                    elif tdi == 2:
                        dec = td.contents[0]
                    elif tdi == 3:
                        ttype = td.contents[0]
                        if ttype != 'sn' and ttype != 'orphan':
                            break
                    elif tdi == 5:
                        if not td.contents:
                            continue
                        ctype = td.contents[0]
                        if ctype == 'Observed':
                            ctype = ''
                    elif tdi == 16:
                        if td.contents:
                            crossrefs = td.findAll('a')
                            for cref in crossrefs:
                                if 'atel' in cref.contents[0].lower():
                                    refs.append([cref.contents[0], cref['href']])
                                elif is_number(cref.contents[0][:4]):
                                    continue
                                else:
                                    aliases.append(cref.contents[0])
    
                if ttype != 'sn' and ttype != 'orphan':
                    break
    
                name = ''
                for alias in aliases:
                    if alias[:2] == 'SN':
                        name = alias
                if not name:
                    name = psname
                name = add_event(name)
                sources = [add_source(name, reference = 'Pan-STARRS 3Pi', url = 'http://psweb.mp.qub.ac.uk/ps1threepi/psdb/')]
                add_quantity(name, 'alias', name, sources[0])
                for ref in refs:
                    sources.append(add_source(name, reference = ref[0], url = ref[1]))
                source = uniq_cdl(sources)
                for alias in aliases:
                    newalias = alias
                    if alias[:3] in ['CSS', 'SSS', 'MLS']:
                        newalias = alias.replace('-', ':', 1)
                    newalias = newalias.replace('PSNJ', 'PSN J')
                    add_quantity(name, 'alias', newalias, source)
                add_quantity(name, 'ra', ra, source)
                add_quantity(name, 'dec', dec, source)
                add_quantity(name, 'claimedtype', ctype, source)
    
                fname2 = '../sne-external/3pi/candidate-' + pslink.rstrip('/').split('/')[-1] + '.html'
                if tasks['psthreepi']['archived'] and os.path.isfile(fname2):
                    with open(fname2, 'r') as f:
                        html2 = f.read()
                elif not offline:
                    pslink = 'http://psweb.mp.qub.ac.uk/ps1threepi/psdb/public/' + pslink
                    with open(fname2, 'w') as f:
                        response2 = urllib.request.urlopen(pslink)
                        html2 = response2.read().decode('utf-8')
                        f.write(html2)
                else:
                    continue
    
                bs2 = BeautifulSoup(html2, "html5lib")
                scripts = bs2.findAll('script')
                nslines = []
                nslabels = []
                for script in scripts:
                    if 'jslcdata.push' not in script.text:
                        continue
                    slines = script.text.splitlines()
                    for line in slines:
                        if 'jslcdata.push' in line:
                            nslines.append(json.loads(line.strip().replace('jslcdata.push(','').replace(');','')))
                        if 'jslabels.push' in line and 'blanks' not in line and 'non det' not in line:
                            nslabels.append(json.loads(line.strip().replace('jslabels.push(','').replace(');',''))['label'])
                for li, line in enumerate(nslines[:len(nslabels)]):
                    if not line:
                        continue
                    for obs in line:
                        add_photometry(name, time = obs[0], band = nslabels[li], magnitude = obs[1], e_magnitude = obs[2], source = source,
                            telescope = 'Pan-STARRS1')
                for li, line in enumerate(nslines[2*len(nslabels):]):
                    if not line:
                        continue
                    for obs in line:
                        add_photometry(name, time = obs[0], band = nslabels[li], magnitude = obs[1], upperlimit = True, source = source,
                            telescope = 'Pan-STARRS1')
                assoctab = bs2.find('table', {"class":"generictable"})
                hostname = ''
                redshift = ''
                if assoctab:
                    trs = assoctab.findAll('tr')
                    headertds = [x.contents[0] for x in trs[1].findAll('td')]
                    tds = trs[1].findAll('td')
                    for tdi, td in enumerate(tds):
                        if tdi == 1:
                            hostname = td.contents[0].strip()
                        elif tdi == 4:
                            if 'z' in headertds:
                                redshift = td.contents[0].strip()
                # Skip galaxies with just SDSS id
                if is_number(hostname):
                    continue
                add_quantity(name, 'host', hostname, source)
                if redshift:
                    add_quantity(name, 'redshift', redshift, source, kind = 'host')
                if args.update:
                    journal_events()
            journal_events()
    
    if do_task(task, 'psmds'):
        with open('../sne-external/MDS/apj506838t1_mrt.txt') as f:
            for ri, row in enumerate(f.read().splitlines()):
                if ri < 35:
                    continue
                cols = [x.strip() for x in row.split(',')]
                name = add_event(cols[0])
                source = add_source(name, bibcode = '2015ApJ...799..208S')
                add_quantity(name, 'alias', name, source)
                add_quantity(name, 'ra', cols[2], source)
                add_quantity(name, 'dec', cols[3], source)
                astrot = astrotime(float(cols[4]), format='mjd').datetime
                add_quantity(name, 'discoverdate', make_date_string(astrot.year, astrot.month, astrot.day), source)
                add_quantity(name, 'redshift', cols[5], source, kind = 'spectroscopic')
                add_quantity(name, 'claimedtype', 'II P', source)
        journal_events()
    
    if do_task(task, 'crts'):
        crtsnameerrors = ['2011ax']
    
        folders = ["catalina", "MLS", "SSS"]
        for fold in tq(folders):
            html = load_cached_url("http://nesssi.cacr.caltech.edu/" + fold + "/AllSN.html", '../sne-external/CRTS/' + fold + '.html')
            if not html:
                continue
            bs = BeautifulSoup(html, "html5lib")
            trs = bs.findAll('tr')
            for tr in tq(trs):
                tds = tr.findAll('td')
                if not tds:
                    continue
                refs = []
                aliases = []
                ttype = ''
                ctype = ''
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
                        aliases = re.sub('[()]', '', re.sub('<[^<]+?>', '', td.contents[-1].strip()))
                        aliases = [x.strip('; ') for x in list(filter(None, aliases.split(' ')))]
    
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
                            if '>' in aliases[ai+1]:
                                hostupper = True
                            hostmag = aliases[ai+1].strip('>~').replace(',', '.')
                        continue
                    if is_number(alias[:4]) and alias[:2] == '20' and len(alias) > 4:
                        name = 'SN' + alias
                    lalias = alias.lower()
                    if (('asassn' in alias and len(alias) > 6) or ('ptf' in alias and len(alias) > 3) or
                        ('ps1' in alias and len(alias) > 3) or 'snhunt' in alias or
                        ('mls' in alias and len(alias) > 3) or 'gaia' in alias or ('lsq' in alias and len(alias) > 3)):
                        alias = alias.replace('SNHunt', 'SNhunt')
                        validaliases.append(alias)
                if not name:
                    name = crtsname
                name = add_event(name)
                source = add_source(name, reference = 'Catalina Sky Survey', bibcode = '2009ApJ...696..870D',
                    url = 'http://nesssi.cacr.caltech.edu/catalina/AllSN.html')
                add_quantity(name, 'alias', name, source)
                for alias in validaliases:
                    add_quantity(name, 'alias', alias, source)
                add_quantity(name, 'ra', ra, source, unit = 'floatdegrees')
                add_quantity(name, 'dec', dec, source, unit = 'floatdegrees')
    
                if hostmag:
                    # 1.0 magnitude error based on Drake 2009 assertion that SN are only considered real if they are 2 mags brighter than host.
                    add_photometry(name, band = 'C', magnitude = hostmag, e_magnitude = 1.0, source = source, host = True,
                        telescope = 'Catalina Schmidt', upperlimit = hostupper)
    
                fname2 = '../sne-external/' + fold + '/' + lclink.split('.')[-2].rstrip('p').split('/')[-1] + '.html'
                if not args.fullrefresh and tasks['crts']['archived'] and os.path.isfile(fname2):
                    with open(fname2, 'r') as f:
                        html2 = f.read()
                else:
                    with open(fname2, 'w') as f:
                        response2 = urllib.request.urlopen(lclink)
                        html2 = response2.read().decode('utf-8')
                        f.write(html2)
    
                lines = html2.splitlines()
                for line in lines:
                    if 'javascript:showx' in line:
                        mjdstr = re.search("showx\('(.*?)'\)", line).group(1).split('(')[0].strip()
                        if not is_number(mjdstr):
                            continue
                        mjd = str(Decimal(mjdstr) + Decimal(53249.0))
                    else:
                        continue
                    if 'javascript:showy' in line:
                        mag = re.search("showy\('(.*?)'\)", line).group(1)
                    if 'javascript:showz' in line:
                        err = re.search("showz\('(.*?)'\)", line).group(1)
                    add_photometry(name, time = mjd, band = 'C', magnitude = mag, source = source, includeshost = (float(err) > 0.0),
                        telescope = 'Catalina Schmidt', e_magnitude = err if float(err) > 0.0 else '', upperlimit = (float(err) == 0.0))
                if args.update:
                    journal_events()
        journal_events()
    
    if do_task(task, 'snhunt'):
        html = load_cached_url('http://nesssi.cacr.caltech.edu/catalina/current.html', '../sne-external/SNhunt/current.html')
        if not html:
            continue
        text = html.splitlines()
        findtable = False
        for ri, row in enumerate(text):
            if 'Supernova Discoveries' in row:
                findtable = True
            if findtable and '<table' in row:
                tstart = ri+1
            if findtable and '</table>' in row:
                tend = ri-1
        tablestr = '<html><body><table>'
        for row in text[tstart:tend]:
            if row[:3] == 'tr>':
                tablestr = tablestr + '<tr>' + row[3:]
            else:
                tablestr = tablestr + row
        tablestr = tablestr + '</table></body></html>'
        bs = BeautifulSoup(tablestr, 'html5lib')
        trs = bs.find('table').findAll('tr')
        for tr in tq(trs):
            cols = [str(x.text) for x in tr.findAll('td')]
            if not cols:
                continue
            name = re.sub('<[^<]+?>', '', cols[4]).strip().replace(' ', '').replace('SNHunt', 'SNhunt')
            name = add_event(name)
            source = add_source(name, reference = 'Supernova Hunt', url = 'http://nesssi.cacr.caltech.edu/catalina/current.html')
            add_quantity(name, 'alias', name, source)
            host = re.sub('<[^<]+?>', '', cols[1]).strip().replace('_', ' ')
            add_quantity(name, 'host', host, source)
            add_quantity(name, 'ra', cols[2], source, unit = 'floatdegrees')
            add_quantity(name, 'dec', cols[3], source, unit = 'floatdegrees')
            dd = cols[0]
            discoverdate = dd[:4] + '/' + dd[4:6] + '/' + dd[6:8]
            add_quantity(name, 'discoverdate', discoverdate, source)
            discoverers = cols[5].split('/')
            for discoverer in discoverers:
                add_quantity(name, 'discoverer', 'CRTS', source)
                add_quantity(name, 'discoverer', discoverer, source)
            if args.update:
                journal_events()
        journal_events()
    
    if do_task(task, 'nedd'): 
        f = open("../sne-external/NED25.12.1-D-10.4.0-20151123.csv", 'r')
        data = csv.reader(f, delimiter=',', quotechar='"')
        reference = "NED-D"
        refurl = "http://ned.ipac.caltech.edu/Library/Distances/"
        nedddict = {}
        oldhostname = ''
        for r, row in enumerate(data):
            if r <= 12:
                continue
            hostname = row[3]
            if args.update and oldhostname != hostname:
                journal_events()
            distmod = row[4]
            moderr = row[5]
            dist = row[6]
            bibcode = unescape(row[8])
            name = ''
            if hostname.startswith('SN '):
                if is_number(hostname[3:7]):
                    name = 'SN' + hostname[3:]
                else:
                    name = hostname[3:]
            elif hostname.startswith('SNLS '):
                name = 'SNLS-' + hostname[5:].split()[0]
            else:
                cleanhost = hostname.replace('MESSIER 0', 'M').replace('MESSIER ', 'M').strip()
                if True in [x in cleanhost for x in ['UGC', 'PGC', 'IC']]:
                    cleanhost = ' '.join([x.lstrip('0') for x in cleanhost.split()])
                if 'ESO' in cleanhost:
                    cleanhost = cleanhost.replace(' ', '').replace('ESO', 'ESO ')
                nedddict.setdefault(cleanhost,[]).append(Decimal(dist))
    
            if name:
                name = add_event(name)
                secondarysource = add_source(name, reference = reference, url = refurl, secondary = True)
                add_quantity(name, 'alias', name, secondarysource)
                if bibcode:
                    source = add_source(name, bibcode = bibcode)
                    sources = uniq_cdl([source, secondarysource])
                else:
                    sources = secondarysource
                add_quantity(name, 'comovingdist', dist, sources)
            oldhostname = hostname
        journal_events()
    
    # Import CPCS
    if do_task(task, 'cpcs'):
        jsontxt = load_cached_url("http://gsaweb.ast.cam.ac.uk/followup/list_of_alerts?format=json&num=100000&published=1&observed_only=1&hashtag=JG_530ad9462a0b8785bfb385614bf178c6",
            "../sne-external/CPCS/index.json")
        if not jsontxt:
            continue
        alertindex = json.loads(jsontxt, object_pairs_hook=OrderedDict)
        ids = [x["id"] for x in alertindex]
    
        for i, ai in enumerate(tq(ids)):
            name = alertindex[i]['ivorn'].split('/')[-1].strip()
            # Skip a few weird entries
            if name == 'ASASSNli':
                continue
            # Just use a whitelist for now since naming seems inconsistent
            if True in [x in name.upper() for x in ['GAIA', 'OGLE', 'ASASSN', 'MASTER', 'OTJ', 'PS1', 'IPTF']]:
                name = name.replace('Verif', '').replace('_', ' ')
                if 'ASASSN' in name and name[6] != '-':
                    name = 'ASASSN-' + name[6:]
                if 'MASTEROTJ' in name:
                    name = name.replace('MASTEROTJ', 'MASTER OT J')
                if 'OTJ' in name:
                    name = name.replace('OTJ', 'MASTER OT J')
                if name.upper().startswith('IPTF'):
                    name = 'iPTF' + name[4:]
                # Only add events that are classified as SN.
                if event_exists(name):
                    continue
                name = add_event(name)
            else:
                continue
    
            secondarysource = add_source(name, reference = 'Cambridge Photometric Calibration Server', url = 'http://gsaweb.ast.cam.ac.uk/followup/', secondary = True)
            add_quantity(name, 'alias', name, secondarysource)
            add_quantity(name, 'ra', str(alertindex[i]['ra']), secondarysource, unit = 'floatdegrees')
            add_quantity(name, 'dec', str(alertindex[i]['dec']), secondarysource, unit = 'floatdegrees')
    
            alerturl = "http://gsaweb.ast.cam.ac.uk/followup/get_alert_lc_data?alert_id=" + str(ai)
            source = add_source(name, reference = 'CPCS Alert ' + str(ai), url = alerturl)
            fname = '../sne-external/CPCS/alert-' + str(ai).zfill(2) + '.json'
            if tasks['cpcs']['archived'] and os.path.isfile(fname):
                with open(fname, 'r') as f:
                    jsonstr = f.read()
            else:
                session = requests.Session()
                response = session.get(alerturl + "&hashtag=JG_530ad9462a0b8785bfb385614bf178c6")
                with open(fname, 'w') as f:
                    jsonstr = response.text
                    f.write(jsonstr)
    
            try:
                cpcsalert = json.loads(jsonstr)
            except:
                continue
    
            mjds = [round_sig(x, sig=9) for x in cpcsalert['mjd']]
            mags = [round_sig(x, sig=6) for x in cpcsalert['mag']]
            errs = [round_sig(x, sig=6) if (is_number(x) and float(x) > 0.0) else '' for x in cpcsalert['magerr']]
            bnds = cpcsalert['filter']
            obs  = cpcsalert['observatory']
            for mi, mjd in enumerate(mjds):
                add_photometry(name, time = mjd, magnitude = mags[mi], e_magnitude = errs[mi],
                    band = bnds[mi], observatory = obs[mi], source = uniq_cdl([source,secondarysource]))
            if args.update:
                journal_events()
        journal_events()
    
    if do_task(task, 'ptf'):
        #response = urllib.request.urlopen("http://wiserep.weizmann.ac.il/objects/list")
        #bs = BeautifulSoup(response, "html5lib")
        #select = bs.find('select', {"name":"objid"})
        #options = select.findAll('option')
        #for option in options:
        #    print(option.text)
        #    name = option.text
        #    if ((name.startswith('PTF') and is_number(name[3:5])) or
        #        name.startswith('PTFS') or name.startswith('iPTF')):
        #        name = add_event(name)
    
        if tasks['ptf']['archived']:
            with open('../sne-external/PTF/update.html', 'r') as f:
                html = f.read()
        else:
            session = requests.Session()
            response = session.get("http://wiserep.weizmann.ac.il/spectra/update")
            html = response.text
            with open('../sne-external/PTF/update.html', 'w') as f:
                f.write(html)
    
        bs = BeautifulSoup(html, "html5lib")
        select = bs.find('select', {"name":"objid"})
        options = select.findAll('option')
        for option in options:
            name = option.text
            if ((name.startswith('PTF') and is_number(name[3:5])) or
                name.startswith('PTFS') or name.startswith('iPTF')):
                if '(' in name:
                    alias = name.split('(')[0].strip(' ')
                    name = name.split('(')[-1].strip(') ').replace('sn', 'SN')
                    name = add_event(name)
                    add_quantity(name, 'alias', alias, source)
                else:
                    name = add_event(name)
        
        with open('../sne-external/PTF/old-ptf-events.csv') as f:
            for suffix in f.read().splitlines():
                name = add_event('PTF' + suffix)
        with open('../sne-external/PTF/perly-2016.csv') as f:
            for row in f.read().splitlines():
                cols = [x.strip() for x in row.split(',')]
                alias = ''
                if cols[8]:
                    name = cols[8]
                    alias = 'PTF' + cols[0]
                else:
                    name = 'PTF' + cols[0]
                name = add_event(name)
                source = add_source(name, bibcode = '2016arXiv160408207D')
                add_quantity(name, 'alias', name, source)
                if alias:
                    add_quantity(name, 'alias', alias, source)
                add_quantity(name, 'ra', cols[1], source)
                add_quantity(name, 'dec', cols[2], source)
                add_quantity(name, 'claimedtype', 'SLSN-' + cols[3], source)
                add_quantity(name, 'redshift', cols[4], source, kind = 'spectroscopic')
                maxdate = cols[6].replace('-', '/')
                add_quantity(name, 'maxdate', maxdate.lstrip('<'), source, upperlimit = maxdate.startswith('<'))
                add_quantity(name, 'ebv', cols[7], source, kind = 'spectroscopic')
                name = add_event('PTF' + suffix)
        journal_events()
    
    if do_task(task, 'asiagospectra'):
        html = load_cached_url("http://sngroup.oapd.inaf.it./cgi-bin/output_class.cgi?sn=1990", "../sne-external-spectra/Asiago/spectra.html")
        if not html:
            continue
        bs = BeautifulSoup(html, "html5lib")
        trs = bs.findAll('tr')
        for tr in tq(trs):
            tds = tr.findAll('td')
            name = ''
            host = ''
            fitsurl = ''
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
                    name = add_event(name)
                    reference = 'Asiago Supernova Catalogue'
                    refurl = 'http://graspa.oapd.inaf.it/cgi-bin/sncat.php'
                    secondarysource = add_source(name, reference = reference, url = refurl, secondary = True)
                    add_quantity(name, 'alias', name, secondarysource)
                    if alias != name:
                        add_quantity(name, 'alias', alias, secondarysource)
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
                elif tdi == 9:
                    epochstr = td.text.strip()
                    if epochstr:
                        mjd = (astrotime(epochstr[:4] + '-' + epochstr[4:6] + '-' + str(floor(float(epochstr[6:]))).zfill(2)).mjd +
                            float(epochstr[6:]) - floor(float(epochstr[6:])))
                    else:
                        mjd = ''
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
                        source = add_source(name, reference = reference, url = refurl)
                    add_quantity(name, 'alias', name, secondarysource)
                    sources = uniq_cdl(list(filter(None, [source, secondarysource])))
                elif tdi == 12:
                    fitslink = td.find('a')
                    if fitslink:
                        fitsurl = fitslink['href']
            if name:
                add_quantity(name, 'claimedtype', claimedtype, sources)
                add_quantity(name, 'ra', ra, sources)
                add_quantity(name, 'dec', dec, sources)
                add_quantity(name, 'redshift', redshift, sources)
                add_quantity(name, 'discoverer', discoverer, sources)
                add_quantity(name, 'host', host, sources)
    
                #if fitsurl:
                #    response = urllib.request.urlopen("http://sngroup.oapd.inaf.it./" + fitsurl)
                #    compressed = io.BytesIO(response.read())
                #    decompressed = gzip.GzipFile(fileobj=compressed)
                #    hdulist = fits.open(decompressed)
                #    scidata = hdulist[0].data
                #    print(hdulist[0].header)
    
                #    print(scidata[3])
                #    sys.exit()
        journal_events()
    
    if do_task(task, 'wiserepspectra'):
        secondaryreference = 'WISeREP'
        secondaryrefurl = 'http://wiserep.weizmann.ac.il/'
        secondarybibcode = '2012PASP..124..668Y'
        wiserepcnt = 0
    
        # These are known to be in error on the WISeREP page, either fix or ignore them.
        wiserepbibcorrectdict = {'2000AJ....120..367G]':'2000AJ....120..367G',
                                 'Harutyunyan+et+al.+2008':'2008A&A...488..383H',
                                 '0609268':'2007AJ....133...58K',
                                 '2006ApJ...636...400Q':'2006ApJ...636..400Q',
                                 '2011ApJ...741...76':'2011ApJ...741...76C',
                                 '2016PASP...128...961':'2016PASP..128...961',
                                 '2002AJ....1124..417H':'2002AJ....1124.417H',
                                 '2013ApJ…774…58D':'2013ApJ...774...58D',
                                 '2011Sci.333..856S':'2011Sci...333..856S',
                                 '2014MNRAS.438,368':'2014MNRAS.438..368T',
                                 '2012MNRAS.420.1135':'2012MNRAS.420.1135S',
                                 '2012Sci..337..942D':'2012Sci...337..942D',
                                 'stt1839':''}
    
        oldname = ''
        for folder in tq(sorted(next(os.walk("../sne-external-WISEREP"))[1], key=lambda s: s.lower())):
            files = glob.glob("../sne-external-WISEREP/" + folder + '/*')
            for fname in tq(files):
                if '.html' in fname:
                    lfiles = deepcopy(files)
                    with open(fname, 'r') as f:
                        path = os.path.abspath(fname)
                        response = urllib.request.urlopen('file://' + path)
                        bs = BeautifulSoup(response, "html5lib")
                        trs = bs.findAll('tr', {'valign': 'top'})
                        for tri, tr in enumerate(trs):
                            if "Click to show/update object" in str(tr.contents):
                                claimedtype = ''
                                instrument = ''
                                epoch = ''
                                observer = ''
                                reducer = ''
                                specfile = ''
                                produceoutput = True
                                specpath = ''
                                tds = tr.findAll('td')
                                for tdi, td in enumerate(tds):
                                    if td.contents:
                                        if tdi == 3:
                                            name = re.sub('<[^<]+?>', '', str(td.contents[0])).strip()
                                        elif tdi == 5:
                                            claimedtype = re.sub('<[^<]+?>', '', str(td.contents[0])).strip()
                                            if claimedtype == 'SN':
                                                claimedtype = ''
                                                continue
                                            if claimedtype[:3] == 'SN ':
                                                claimedtype = claimedtype[3:].strip()
                                            claimedtype = claimedtype.replace('-like', '').strip()
                                        elif tdi == 9:
                                            instrument = re.sub('<[^<]+?>', '', str(td.contents[0])).strip()
                                        elif tdi == 11:
                                            epoch = re.sub('<[^<]+?>', '', str(td.contents[0])).strip()
                                        elif tdi == 13:
                                            observer = re.sub('<[^<]+?>', '', str(td.contents[0])).strip()
                                            if observer == 'Unknown' or observer == 'Other':
                                                observer = ''
                                        elif tdi == 17:
                                            reducer = re.sub('<[^<]+?>', '', str(td.contents[0])).strip()
                                            if reducer == 'Unknown' or reducer == 'Other':
                                                reducer = ''
                                        elif tdi == 25:
                                            speclinks = td.findAll('a')
                                            try:
                                                for link in speclinks:
                                                    if 'Ascii' in link['href']:
                                                        specfile = link.contents[0].strip()
                                                        tfiles = deepcopy(lfiles)
                                                        for fi, fname in enumerate(lfiles):
                                                            if specfile in fname:
                                                                specpath = fname
                                                                del(tfiles[fi])
                                                                lfiles = deepcopy(tfiles)
                                                                raise(StopIteration)
                                            except StopIteration:
                                                pass
                                            if not specpath:
                                                warnings.warn('Spectrum file not found, "' + specfile + '"')
                                        else:
                                            continue
                            if "Spec Type:</span>" in str(tr.contents) and produceoutput:
                                produceoutput = False
    
                                trstr = str(tr)
                                result = re.search('redshift=(.*?)&amp;', trstr)
                                redshift = ''
                                if result:
                                    redshift = result.group(1)
                                    if not is_number(redshift) or float(redshift) > 100.:
                                        redshift = ''
    
                                result = re.search('publish=(.*?)&amp;', trstr)
                                bibcode = ''
                                if result:
                                    bibcode = unescape(urllib.parse.unquote(urllib.parse.unquote(result.group(1))).split('/')[-1])
    
                                if not bibcode:
                                    biblink = tr.find('a', {'title': 'Link to NASA ADS'})
                                    if biblink:
                                        bibcode = biblink.contents[0]
    
                                if name.startswith('sn'):
                                    name = 'SN' + name[2:]
                                if name.startswith(('CSS', 'SSS', 'MLS')) and ':' not in name:
                                    name = name.replace('-', ':', 1)
                                if name.startswith('MASTERJ'):
                                    name = name.replace('MASTERJ', 'MASTER OT J')
                                if name.startswith('PSNJ'):
                                    name = name.replace('PSNJ', 'PSN J')
                                name = get_preferred_name(name)
                                if oldname and name != oldname:
                                    journal_events()
                                oldname = name
                                name = add_event(name)
    
                                #print(name + " " + claimedtype + " " + epoch + " " + observer + " " + reducer + " " + specfile + " " + bibcode + " " + redshift)
    
                                secondarysource = add_source(name, reference = secondaryreference, url = secondaryrefurl, bibcode = secondarybibcode, secondary = True)
                                add_quantity(name, 'alias', name, secondarysource)
                                if bibcode:
                                    newbibcode = bibcode
                                    if bibcode in wiserepbibcorrectdict:
                                        newbibcode = wiserepbibcorrectdict[bibcode]
                                    if newbibcode:
                                        source = add_source(name, bibcode = unescape(newbibcode))
                                    else:
                                        source = add_source(name, reference = unescape(bibcode))
                                    sources = uniq_cdl([source, secondarysource])
                                else:
                                    sources = secondarysource
    
                                if claimedtype not in ['Other']:
                                    add_quantity(name, 'claimedtype', claimedtype, secondarysource)
                                add_quantity(name, 'redshift', redshift, secondarysource)
    
                                if not specpath:
                                    continue
    
                                with open(specpath,'r') as f:
                                    data = [x.split() for x in f]
    
                                    skipspec = False
                                    newdata = []
                                    oldval = ''
                                    for row in data:
                                        if row and '#' not in row[0]:
                                            if len(row) >= 2 and is_number(row[0]) and is_number(row[1]) and row[1] != oldval:
                                                newdata.append(row)
                                                oldval = row[1]
    
                                    if skipspec or not newdata:
                                        warnings.warn('Skipped adding spectrum file ' + specfile)
                                        continue
    
                                    data = [list(i) for i in zip(*newdata)]
                                    wavelengths = data[0]
                                    fluxes = data[1]
                                    errors = ''
                                    if len(data) == 3:
                                        errors = data[1]
                                    time = str(astrotime(epoch).mjd)
    
                                    if max([float(x) for x in fluxes]) < 1.0e-5:
                                        fluxunit = 'erg/s/cm^2/Angstrom'
                                    else:
                                        fluxunit = 'Uncalibrated'
    
                                    add_spectrum(name = name, waveunit = 'Angstrom', fluxunit = fluxunit, errors = errors, errorunit = fluxunit, wavelengths = wavelengths,
                                        fluxes = fluxes, u_time = 'MJD', time = time, instrument = instrument, source = sources, observer = observer, reducer = reducer,
                                        filename = specfile)
                                    wiserepcnt = wiserepcnt + 1
    
                                    if args.travis and wiserepcnt % travislimit == 0:
                                        break
    
                    tprint('Unadded files: ' + str(len(lfiles) - 1) + "/" + str(len(files)-1))
                    tprint('WISeREP spectrum count: ' + str(wiserepcnt))
        journal_events()
    
    if do_task(task, 'cfaspectra'): 
        # Ia spectra
        oldname = ''
        for name in sorted(next(os.walk("../sne-external-spectra/CfA_SNIa"))[1], key=lambda s: s.lower()):
            fullpath = "../sne-external-spectra/CfA_SNIa/" + name
            origname = name
            if name.startswith('sn') and is_number(name[2:6]):
                name = 'SN' + name[2:]
            if name.startswith('snf') and is_number(name[3:7]):
                name = 'SNF' + name[3:]
            name = get_preferred_name(name)
            if oldname and name != oldname:
                journal_events()
            oldname = name
            name = add_event(name)
            reference = 'CfA Supernova Archive'
            refurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
            source = add_source(name, reference = reference, url = refurl, secondary = True)
            add_quantity(name, 'alias', name, source)
            for fi, fname in enumerate(sorted(glob.glob(fullpath + '/*'), key=lambda s: s.lower())):
                filename = os.path.basename(fname)
                fileparts = filename.split('-')
                if origname.startswith("sn") and is_number(origname[2:6]):
                    year = fileparts[1][:4]
                    month = fileparts[1][4:6]
                    day = fileparts[1][6:]
                    instrument = fileparts[2].split('.')[0]
                else:
                    year = fileparts[2][:4]
                    month = fileparts[2][4:6]
                    day = fileparts[2][6:]
                    instrument = fileparts[3].split('.')[0]
                time = astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day))
                f = open(fname,'r')
                data = csv.reader(f, delimiter=' ', skipinitialspace=True)
                data = [list(i) for i in zip(*data)]
                wavelengths = data[0]
                fluxes = data[1]
                errors = data[2]
                sources = uniq_cdl([source, add_source(name, bibcode = '2012AJ....143..126B'), add_source(name, bibcode = '2008AJ....135.1598M')])
                add_spectrum(name = name, waveunit = 'Angstrom', fluxunit = 'erg/s/cm^2/Angstrom', filename = filename,
                    wavelengths = wavelengths, fluxes = fluxes, u_time = 'MJD', time = time, instrument = instrument,
                    errorunit = "ergs/s/cm^2/Angstrom", errors = errors, source = sources, dereddened = False, deredshifted = False)
                if args.travis and fi >= travislimit:
                    break
        journal_events()
    
        # Ibc spectra
        oldname = ''
        for name in sorted(next(os.walk("../sne-external-spectra/CfA_SNIbc"))[1], key=lambda s: s.lower()):
            fullpath = "../sne-external-spectra/CfA_SNIbc/" + name
            if name.startswith('sn') and is_number(name[2:6]):
                name = 'SN' + name[2:]
            name = get_preferred_name(name)
            if oldname and name != oldname:
                journal_events()
            oldname = name
            name = add_event(name)
            reference = 'CfA Supernova Archive'
            refurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
            source = add_source(name, reference = reference, url = refurl, secondary = True)
            add_quantity(name, 'alias', name, source)
            for fi, fname in enumerate(sorted(glob.glob(fullpath + '/*'), key=lambda s: s.lower())):
                filename = os.path.basename(fname)
                fileparts = filename.split('-')
                instrument = ''
                year = fileparts[1][:4]
                month = fileparts[1][4:6]
                day = fileparts[1][6:].split('.')[0]
                if len(fileparts) > 2:
                    instrument = fileparts[-1].split('.')[0]
                time = astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day))
                f = open(fname,'r')
                data = csv.reader(f, delimiter=' ', skipinitialspace=True)
                data = [list(i) for i in zip(*data)]
                wavelengths = data[0]
                fluxes = data[1]
                sources = uniq_cdl([source, add_source(name, bibcode = '2014AJ....147...99M')])
                add_spectrum(name = name, waveunit = 'Angstrom', fluxunit = 'erg/s/cm^2/Angstrom', wavelengths = wavelengths, filename = filename,
                    fluxes = fluxes, u_time = 'MJD', time = time, instrument = instrument, source = sources,
                    dereddened = False, deredshifted = False)
                if args.travis and fi >= travislimit:
                    break
        journal_events()
    
    if do_task(task, 'snlsspectra'): 
        result = Vizier.get_catalogs("J/A+A/507/85/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        datedict = {}
        for row in table:
            datedict['SNLS-' + row['SN']] = str(astrotime(row['Date']).mjd)
    
        oldname = ''
        for fi, fname in enumerate(sorted(glob.glob('../sne-external-spectra/SNLS/*'), key=lambda s: s.lower())):
            filename = os.path.basename(fname)
            fileparts = filename.split('_')
            name = 'SNLS-' + fileparts[1]
            name = get_preferred_name(name)
            if oldname and name != oldname:
                journal_events()
            oldname = name
            name = add_event(name)
            source = add_source(name, bibcode = "2009A&A...507...85B")
            add_quantity(name, 'alias', name, source)
    
            add_quantity(name, 'discoverdate', '20' + fileparts[1][:2], source)
    
            f = open(fname,'r')
            data = csv.reader(f, delimiter=' ', skipinitialspace=True)
            specdata = []
            for r, row in enumerate(data):
                if row[0] == '@TELESCOPE':
                    telescope = row[1].strip()
                elif row[0] == '@REDSHIFT':
                    add_quantity(name, 'redshift', row[1].strip(), source)
                if r < 14:
                    continue
                specdata.append(list(filter(None, [x.strip(' \t') for x in row])))
            specdata = [list(i) for i in zip(*specdata)]
            wavelengths = specdata[1]
            
            fluxes = [pretty_num(float(x)*1.e-16, sig = get_sig_digits(x)) for x in specdata[2]]
            errors = [pretty_num(float(x)*1.e-16, sig = get_sig_digits(x)) for x in specdata[3]]
    
            add_spectrum(name = name, waveunit = 'Angstrom', fluxunit = 'erg/s/cm^2/Angstrom', wavelengths = wavelengths,
                fluxes = fluxes, u_time = 'MJD' if name in datedict else '', time = datedict[name] if name in datedict else '', telescope = telescope, source = source,
                filename = filename)
            if args.travis and fi >= travislimit:
                break
        journal_events()
    
    if do_task(task, 'cspspectra'): 
        oldname = ''
        for fi, fname in enumerate(sorted(glob.glob('../sne-external-spectra/CSP/*'), key=lambda s: s.lower())):
            filename = os.path.basename(fname)
            sfile = filename.split('.')
            if sfile[1] == 'txt':
                continue
            sfile = sfile[0]
            fileparts = sfile.split('_')
            name = 'SN20' + fileparts[0][2:]
            name = get_preferred_name(name)
            if oldname and name != oldname:
                journal_events()
            oldname = name
            name = add_event(name)
            telescope = fileparts[-2]
            instrument = fileparts[-1]
            source = add_source(name, bibcode = "2013ApJ...773...53F")
            add_quantity(name, 'alias', name, source)
    
            f = open(fname,'r')
            data = csv.reader(f, delimiter=' ', skipinitialspace=True)
            specdata = []
            for r, row in enumerate(data):
                if row[0] == '#JDate_of_observation:':
                    jd = row[1].strip()
                    time = str(jd_to_mjd(Decimal(jd)))
                elif row[0] == '#Redshift:':
                    add_quantity(name, 'redshift', row[1].strip(), source)
                if r < 7:
                    continue
                specdata.append(list(filter(None, [x.strip(' ') for x in row])))
            specdata = [list(i) for i in zip(*specdata)]
            wavelengths = specdata[0]
            fluxes = specdata[1]
    
            add_spectrum(name = name, u_time = 'MJD', time = time, waveunit = 'Angstrom', fluxunit = 'erg/s/cm^2/Angstrom', wavelengths = wavelengths,
                fluxes = fluxes, telescope = telescope, instrument = instrument, source = source, deredshifted = True, filename = filename)
            if args.travis and fi >= travislimit:
                break
        journal_events()
    
    if do_task(task, 'ucbspectra'):
        secondaryreference = "UCB Filippenko Group's Supernova Database (SNDB)"
        secondaryrefurl = "http://heracles.astro.berkeley.edu/sndb/info"
        secondaryrefbib = "2012MNRAS.425.1789S"
        ucbspectracnt = 0
    
        jsontxt = load_cached_url("http://heracles.astro.berkeley.edu/sndb/download?id=allpubspec",
            '../sne-external-spectra/UCB/allpub.json')
        if not jsontxt:
            continue
    
        spectra = json.loads(jsontxt)
        spectra = sorted(spectra, key = lambda k: k['ObjName'])
        oldname = ''
        for spectrum in spectra:
            name = spectrum["ObjName"]
            if oldname and name != oldname:
                journal_events()
            oldname = name
            name = add_event(name)
    
            secondarysource = add_source(name, reference = secondaryreference, url = secondaryrefurl, bibcode = secondaryrefbib, secondary = True)
            add_quantity(name, 'alias', name, secondarysource)
            sources = [secondarysource]
            if spectrum["Reference"]:
                sources += [add_source(name, bibcode = spectrum["Reference"])]
            sources = uniq_cdl(sources)
    
            if spectrum["SNID_Subtype"] and spectrum["SNID_Subtype"].strip() != "NoMatch":
                for ct in spectrum["SNID_Subtype"].strip().split(','):
                    add_quantity(name, 'claimedtype', ct.replace('-norm', '').strip(), sources)
            if spectrum["DiscDate"]:
                add_quantity(name, 'discoverdate', spectrum["DiscDate"].replace('-', '/'), sources)
            if spectrum["HostName"]:
                add_quantity(name, 'host', spectrum["HostName"], sources)
            if spectrum["UT_Date"]:
                epoch = str(spectrum["UT_Date"])
                year = epoch[:4]
                month = epoch[4:6]
                day = epoch[6:]
                sig = get_sig_digits(day) + 5
                mjd = pretty_num(astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day)), sig = sig)
            filename = spectrum["Filename"] if spectrum["Filename"] else ''
            instrument = spectrum["Instrument"] if spectrum["Instrument"] else ''
            reducer = spectrum["Reducer"] if spectrum["Reducer"] else ''
            observer = spectrum["Observer"] if spectrum["Observer"] else ''
            snr = str(spectrum["SNR"]) if spectrum["SNR"] else ''
    
            if not filename:
                raise(ValueError('Filename not found for SNDB spectrum!'))
            if not spectrum["SpecID"]:
                raise(ValueError('ID not found for SNDB spectrum!'))
    
            filepath = '../sne-external-spectra/UCB/' + filename
            if tasks['ucbspectra']['archived'] and os.path.isfile(filepath):
                with open(filepath, 'r') as f:
                    spectxt = f.read()
            else:
                session = requests.Session()
                response = session.get("http://heracles.astro.berkeley.edu/sndb/download?id=ds:" + str(spectrum["SpecID"]))
                spectxt = response.text
                with open(filepath, 'w') as f:
                    f.write(spectxt)
    
            specdata = list(csv.reader(spectxt.splitlines(), delimiter=' ', skipinitialspace=True))
            startrow = 0
            for row in specdata:
                if row[0][0] == '#':
                    startrow += 1
                else:
                    break
            specdata = specdata[startrow:]
    
            haserrors = len(specdata[0]) == 3 and specdata[0][2] and specdata[0][2] != 'NaN'
            specdata = [list(i) for i in zip(*specdata)]
    
            wavelengths = specdata[0]
            fluxes = specdata[1]
            errors = ''
            if haserrors:
                errors = specdata[2]
    
            if not list(filter(None, errors)):
                errors = ''
    
            add_spectrum(name = name, u_time = 'MJD', time = mjd, waveunit = 'Angstrom', fluxunit = 'Uncalibrated',
                wavelengths = wavelengths, filename = filename, fluxes = fluxes, errors = errors, errorunit = 'Uncalibrated',
                instrument = instrument, source = sources, snr = snr, observer = observer, reducer = reducer,
                deredshifted = ('-noz' in filename))
            ucbspectracnt = ucbspectracnt + 1
            if args.travis and ucbspectracnt >= travislimit:
                break
        journal_events()
    
    if do_task(task, 'suspectspectra'): 
        with open('../sne-external-spectra/Suspect/sources.json', 'r') as f:
            sourcedict = json.loads(f.read())
    
        with open('../sne-external-spectra/Suspect/filename-changes.txt', 'r') as f:
            rows = f.readlines()
            changedict = {}
            for row in rows:
                if not row.strip() or row[0] == "#":
                    continue
                items = row.strip().split(' ')
                changedict[items[1]] = items[0]
    
        suspectcnt = 0
        folders = next(os.walk('../sne-external-spectra/Suspect'))[1]
        for folder in folders:
            eventfolders = next(os.walk('../sne-external-spectra/Suspect/'+folder))[1]
            oldname = ''
            for eventfolder in eventfolders:
                name = eventfolder
                if is_number(name[:4]):
                    name = 'SN' + name
                name = get_preferred_name(name)
                if oldname and name != oldname:
                    journal_events()
                oldname = name
                name = add_event(name)
                secondaryreference = "SUSPECT"
                secondaryrefurl = "https://www.nhn.ou.edu/~suspect/"
                secondarybibcode = "2001AAS...199.8408R"
                secondarysource = add_source(name, reference = secondaryreference, url = secondaryrefurl, bibcode = secondarybibcode, secondary = True)
                add_quantity(name, 'alias', name, secondarysource)
                eventspectra = next(os.walk('../sne-external-spectra/Suspect/'+folder+'/'+eventfolder))[2]
                for spectrum in eventspectra:
                    sources = [secondarysource]
                    bibcode = ''
                    if spectrum in changedict:
                        specalias = changedict[spectrum]
                    else:
                        specalias = spectrum
                    if specalias in sourcedict:
                        bibcode = sourcedict[specalias]
                    elif name in sourcedict:
                        bibcode = sourcedict[name]
                    if bibcode:
                        source = add_source(name, bibcode = unescape(bibcode))
                        sources += source
                    sources = uniq_cdl(sources)
    
                    date = spectrum.split('_')[1]
                    year = date[:4]
                    month = date[4:6]
                    day = date[6:]
                    sig = get_sig_digits(day) + 5
                    time = pretty_num(astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day)), sig = sig)
    
                    with open('../sne-external-spectra/Suspect/'+folder+'/'+eventfolder+'/'+spectrum) as f:
                        specdata = list(csv.reader(f, delimiter=' ', skipinitialspace=True))
                        specdata = list(filter(None, specdata))
                        newspec = []
                        oldval = ''
                        for row in specdata:
                            if row[1] == oldval:
                                continue
                            newspec.append(row)
                            oldval = row[1]
                        specdata = newspec
                    haserrors = len(specdata[0]) == 3 and specdata[0][2] and specdata[0][2] != 'NaN'
                    specdata = [list(i) for i in zip(*specdata)]
    
                    wavelengths = specdata[0]
                    fluxes = specdata[1]
                    errors = ''
                    if haserrors:
                        errors = specdata[2]
    
                    add_spectrum(name = name, u_time = 'MJD', time = time, waveunit = 'Angstrom', fluxunit = 'Uncalibrated', wavelengths = wavelengths,
                        fluxes = fluxes, errors = errors, errorunit = 'Uncalibrated', source = sources, filename = spectrum)
                    suspectcnt = suspectcnt + 1
                    if args.travis and suspectcnt % travislimit == 0:
                        break
        journal_events()
    
    if do_task(task, 'snfspectra'): 
        eventfolders = next(os.walk('../sne-external-spectra/SNFactory'))[1]
        bibcodes = {'SN2005gj':'2006ApJ...650..510A', 'SN2006D':'2007ApJ...654L..53T', 'SN2007if':'2010ApJ...713.1073S', 'SN2011fe':'2013A&A...554A..27P'}
        oldname = ''
        snfcnt = 0
        for eventfolder in eventfolders:
            name = eventfolder
            name = get_preferred_name(name)
            if oldname and name != oldname:
                journal_events()
            oldname = name
            name = add_event(name)
            secondaryreference = "Nearby Supernova Factory"
            secondaryrefurl = "http://snfactory.lbl.gov/"
            secondarybibcode = "2002SPIE.4836...61A"
            secondarysource = add_source(name, reference = secondaryreference, url = secondaryrefurl, bibcode = secondarybibcode, secondary = True)
            add_quantity(name, 'alias', name, secondarysource)
            bibcode = bibcodes[name]
            source = add_source(name, bibcode = bibcode)
            sources = uniq_cdl([source,secondarysource])
            eventspectra = glob.glob('../sne-external-spectra/SNFactory/'+eventfolder+'/*.dat')
            for spectrum in eventspectra:
                filename = os.path.basename(spectrum)
                with open(spectrum) as f:
                    specdata = list(csv.reader(f, delimiter=' ', skipinitialspace=True))
                specdata = list(filter(None, specdata))
                newspec = []
                time = ''
                telescope = ''
                instrument = ''
                observer = ''
                observatory = ''
                if 'Keck_20060202_R' in spectrum:
                    time = '53768.23469'
                elif 'Spectrum05_276' in spectrum:
                    time = pretty_num(astrotime('2005-10-03').mjd, sig = 5)
                elif 'Spectrum05_329' in spectrum:
                    time = pretty_num(astrotime('2005-11-25').mjd, sig = 5)
                elif 'Spectrum05_336' in spectrum:
                    time = pretty_num(astrotime('2005-12-02').mjd, sig = 5)
                for row in specdata:
                    if row[0][0] == '#':
                        joinrow = (' '.join(row)).split('=')
                        if len(joinrow) < 2:
                            continue
                        field = joinrow[0].strip('# ')
                        value = joinrow[1].split('/')[0].strip("' ")
                        if not time:
                            if field == 'JD':
                                time = str(jd_to_mjd(Decimal(value)))
                            elif field == 'MJD':
                                time = value
                            elif field == 'MJD-OBS':
                                time = value
                        if field == 'OBSERVER':
                            observer = value.capitalize()
                        if field == 'OBSERVAT':
                            observatory = value.capitalize()
                        if field == 'TELESCOP':
                            telescope = value.capitalize()
                        if field == 'INSTRUME':
                            instrument = value.capitalize()
                    else:
                        newspec.append(row)
                if not time:
                    raise(ValueError('Time missing from spectrum.'))
                specdata = newspec
                haserrors = len(specdata[0]) == 3 and specdata[0][2] and specdata[0][2] != 'NaN'
                specdata = [list(i) for i in zip(*specdata)]
    
                wavelengths = specdata[0]
                fluxes = specdata[1]
                errors = ''
                if haserrors:
                    errors = specdata[2]
    
                add_spectrum(name = name, u_time = 'MJD', time = time, waveunit = 'Angstrom', fluxunit = 'erg/s/cm^2/Angstrom',
                    wavelengths = wavelengths, fluxes = fluxes, errors = errors, observer = observer, observatory = observatory,
                    telescope = telescope, instrument = instrument,
                    errorunit = ('Variance' if name == 'SN2011fe' else 'erg/s/cm^2/Angstrom'), source = sources, filename = filename)
                snfcnt = snfcnt + 1
                if args.travis and snfcnt % travislimit == 0:
                    break
        journal_events()
    
    if do_task(task, 'superfitspectra'):
        sfdirs = glob.glob('../sne-external-spectra/superfit/*')
        for sfdir in sfdirs:
            sffiles = sorted(glob.glob(sfdir + "/*.dat"))
            lastname = ''
            oldname = ''
            for sffile in sffiles:
                basename = os.path.basename(sffile)
                name = basename.split('.')[0]
                if name.startswith('sn'):
                    name = 'SN' + name[2:]
                    if len(name) == 7:
                        name = name[:6] + name[6].upper()
                elif name.startswith('ptf'):
                    name = 'PTF' + name[3:]
    
                if 'theory' in name:
                    continue
                if event_exists(name):
                    prefname = get_preferred_name(name)
                    if 'spectra' in events[prefname] and lastname != prefname:
                        continue
                if oldname and name != oldname:
                    journal_events()
                oldname = name
                name = add_event(name)
                epoch = basename.split('.')[1]
                (mldt, mlmag, mlband, mlsource) = get_max_light(name)
                if mldt:
                    epoff = Decimal(0.0) if epoch == 'max' else (Decimal(epoch[1:]) if epoch[0] == 'p' else -Decimal(epoch[1:]))
                else:
                    epoff = ''
    
                source = add_source(name, reference = 'Superfit', url = 'http://www.dahowell.com/superfit.html', secondary = True)
                add_quantity(name, 'alias', name, source)
    
                with open(sffile) as f:
                    rows = f.read().splitlines()
                specdata = []
                for row in rows:
                    if row.strip():
                        specdata.append(list(filter(None,re.split('\t+|\s+', row, maxsplit=0))))
                specdata = [[x.replace('D','E') for x in list(i)] for i in zip(*specdata)]
                wavelengths = specdata[0]
                fluxes = specdata[1]
    
                mlmjd = str(Decimal(astrotime('-'.join([str(mldt.year), str(mldt.month), str(mldt.day)])).mjd) + epoff) if (epoff != '') else ''
                add_spectrum(name, u_time = 'MJD' if mlmjd else '', time = mlmjd, waveunit = 'Angstrom', fluxunit = 'Uncalibrated',
                    wavelengths = wavelengths, fluxes = fluxes, source = source)
                
                lastname = name
            journal_events()

if args.update and not len(events):
    tprint('No sources changed, event files unchanged in update.')
    sys.exit()

merge_duplicates()
set_preferred_names()

files = []
for rep in repofolders:
    files += glob.glob('../' + rep + "/*.json")

path = '../bibauthors.json'
if os.path.isfile(path):
    with open(path, 'r') as f:
        bibauthordict = json.loads(f.read(), object_pairs_hook=OrderedDict)
else:
    bibauthordict = OrderedDict()
path = '../extinctions.json'
if os.path.isfile(path):
    with open(path, 'r') as f:
        extinctionsdict = json.loads(f.read(), object_pairs_hook=OrderedDict)
else:
    extinctionsdict = OrderedDict()

for fi in tqdm(files, desc = 'Sanitizing and deriving quantities for events'):
    events = OrderedDict()
    name = os.path.basename(os.path.splitext(fi)[0])
    name = add_event(name, loadifempty = False)
    derive_and_sanitize()
    if has_task('writeevents'): 
        write_all_events(empty = True, gz = True, delete = True)

jsonstring = json.dumps(bibauthordict, indent='\t', separators=(',', ':'), ensure_ascii=False)
with codecs.open('../bibauthors.json', 'w', encoding='utf8') as f:
    f.write(jsonstring)
jsonstring = json.dumps(extinctionsdict, indent='\t', separators=(',', ':'), ensure_ascii=False)
with codecs.open('../extinctions.json', 'w', encoding='utf8') as f:
    f.write(jsonstring)

print("Memory used (MBs on Mac, GBs on Linux): " + "{:,}".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024./1024.))
