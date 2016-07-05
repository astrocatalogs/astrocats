#!/usr/local/bin/python3.5

import csv
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
import statistics
import warnings
import subprocess
from datetime import timedelta, datetime
from glob import glob
from hashlib import md5
from html import unescape
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
from math import log10, floor, sqrt, isnan, ceil, hypot, pi
from bs4 import BeautifulSoup, Tag, NavigableString
from string import ascii_letters
from photometry import *
from tq import *
from digits import *
from repos import *
from events import *

parser = argparse.ArgumentParser(description='Generate a catalog JSON file and plot HTML files from SNE data.')
parser.add_argument('--update', '-u',       dest='update',      help='Only update catalog using live sources.',    default=False, action='store_true')
parser.add_argument('--verbose', '-v',      dest='verbose',     help='Print more messages to the screen.',         default=False, action='store_true')
parser.add_argument('--refresh', '-r',      dest='refresh',     help='Ignore most task caches.',                   default=False, action='store_true')
parser.add_argument('--full-refresh', '-f', dest='fullrefresh', help='Ignore all task caches.',                    default=False, action='store_true')
parser.add_argument('--archived', '-a',     dest='archived',    help='Always use task caches.',                    default=False, action='store_true')
parser.add_argument('--travis', '-tr',      dest='travis',      help='Run import script in test mode for Travis.', default=False, action='store_true')
parser.add_argument('--refreshlist', '-rl', dest='refreshlist', help='Comma-delimited list of caches to clear.',   default='')
args = parser.parse_args()

tasks = OrderedDict([
    ("deleteoldevents", {"nicename":"Deleting old events",          "update": False}),
    ("internal",        {"nicename":"%pre metadata and photometry", "update": False}),
    ("radio",           {"nicename":"%pre radio data",              "update": False}),
    ("xray",            {"nicename":"%pre X-ray data",              "update": False}),
    ("simbad",          {"nicename":"%pre SIMBAD",                  "update": False}),
    ("vizier",          {"nicename":"%pre VizieR",                  "update": False}),
    ("donations",       {"nicename":"%pre donations",               "update": False}),
    ("pessto-dr1",      {"nicename":"%pre PESSTO DR1",              "update": False}),
    ("scp",             {"nicename":"%pre SCP",                     "update": False}),
    ("ascii",           {"nicename":"%pre ASCII",                   "update": False}),
    ("cccp",            {"nicename":"%pre CCCP",                    "update": False, "archived": True}),
    ("suspect",         {"nicename":"%pre SUSPECT",                 "update": False}),
    ("cfa",             {"nicename":"%pre CfA archive photometry",  "update": False}),
    ("ucb",             {"nicename":"%pre UCB photometry",          "update": False, "archived": True}),
    ("sdss",            {"nicename":"%pre SDSS photometry",         "update": False}),
    ("csp",             {"nicename":"%pre CSP photometry",          "update": False}),
    ("itep",            {"nicename":"%pre ITEP",                    "update": False}),
    ("asiago",          {"nicename":"%pre Asiago metadata",         "update": False}),
    ("tns",             {"nicename":"%pre TNS metadata",            "update": True,  "archived": True}),
    ("rochester",       {"nicename":"%pre Latest Supernovae",       "update": True,  "archived": False}),
    ("lennarz",         {"nicename":"%pre Lennarz",                 "update": False}),
    ("fermi",           {"nicename":"%pre Fermi",                   "update": False}),
    ("gaia",            {"nicename":"%pre GAIA",                    "update": True,  "archived": False}),
    ("ogle",            {"nicename":"%pre OGLE",                    "update": True,  "archived": False}),
    ("snls",            {"nicename":"%pre SNLS",                    "update": False}),
    ("psthreepi",       {"nicename":"%pre Pan-STARRS 3π",           "update": True,  "archived": False}),
    ("psmds",           {"nicename":"%pre Pan-STARRS MDS",          "update": False}),
    ("psst",            {"nicename":"%pre PSST",                    "update": False}),
    ("grb",             {"nicename":"%pre GRB catalog",             "update": True,  "archived": False}),
    ("crts",            {"nicename":"%pre CRTS",                    "update": True,  "archived": False}),
    ("snhunt",          {"nicename":"%pre SNhunt",                  "update": True,  "archived": False}),
    ("nedd",            {"nicename":"%pre NED-D",                   "update": False}),
    ("cpcs",            {"nicename":"%pre CPCS",                    "update": True,  "archived": False}),
    ("ptf",             {"nicename":"%pre PTF",                     "update": False, "archived": False}),
    ("des",             {"nicename":"%pre DES",                     "update": False, "archived": False}),
    ("asassn",          {"nicename":"%pre ASASSN",                  "update": True }),
    ("snf",             {"nicename":"%pre SNF",                     "update": False}),
    ("asiagospectra",   {"nicename":"%pre Asiago spectra",          "update": True }),
    ("wiserepspectra",  {"nicename":"%pre WISeREP spectra",         "update": False}),
    ("cfaspectra",      {"nicename":"%pre CfA archive spectra",     "update": False}),
    ("snlsspectra",     {"nicename":"%pre SNLS spectra",            "update": False}),
    ("cspspectra",      {"nicename":"%pre CSP spectra",             "update": False}),
    ("ucbspectra",      {"nicename":"%pre UCB spectra",             "update": True,  "archived": True}),
    ("suspectspectra",  {"nicename":"%pre SUSPECT spectra",         "update": False}),
    ("snfspectra",      {"nicename":"%pre SNH spectra",             "update": False}),
    ("superfitspectra", {"nicename":"%pre Superfit spectra",        "update": False}),
    ("mergeduplicates", {"nicename":"Merging duplicates",           "update": False}),
    ("setprefnames",    {"nicename":"Setting preferred names",      "update": False}),
    ("writeevents",     {"nicename":"Writing events",               "update": True })
])

oscbibcode = '2016arXiv160501054G'
oscname = 'The Open Supernova Catalog'
oscurl = 'https://sne.space'

cfaack = ("This research has made use of the CfA Supernova Archive, "
          "which is funded in part by the National Science Foundation "
          "through grant AST 0907903.")

clight = const.c.cgs.value
km = (1.0 * un.km).cgs.value
planckh = const.h.cgs.value
keV = (1.0 * un.keV).cgs.value
travislimit = 10

currenttask = ''

eventnames = []
events = OrderedDict()

warnings.filterwarnings('ignore', r'Warning: converting a masked element to nan.')

with open('type-synonyms.json', 'r') as f:
    typereps = json.loads(f.read(), object_pairs_hook=OrderedDict)
with open('source-synonyms.json', 'r') as f:
    sourcereps = json.loads(f.read(), object_pairs_hook=OrderedDict)
with open('url-redirects.json', 'r') as f:
    urlreps = json.loads(f.read(), object_pairs_hook=OrderedDict)
with open('non-sne-types.json', 'r') as f:
    nonsnetypes = json.loads(f.read(), object_pairs_hook=OrderedDict)
    nonsnetypes = [x.upper() for x in nonsnetypes]

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
    ['V', 'G'],      # if not, V-like bands
    ['R', 'r']       # if not, R-like bands
]

gitrevhash = subprocess.check_output(['git', 'log', '-n', '1', '--format="%H"',
    '--', '../OSC-JSON-format.md']).decode('ascii').strip().strip('"').strip()
def get_schema():
    return 'https://github.com/astrocatalogs/sne/blob/' + gitrevhash + '/OSC-JSON-format.md'

def uniq_cdl(values):
    return ','.join(sorted(list(set(values))))

def rep_chars(string, chars, rep = ''):
    for c in chars:
        if c in string:
            string = string.replace(c, rep)
    return string

def single_spaces(string):
    return ' '.join(list(filter(None, string.split())))

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

def get_source_year(source):
    if 'bibcode' in source:
        if is_number(source['bibcode'][:4]):
            return int(source['bibcode'][:4])
        else:
            return -10000
    raise(ValueError('No bibcode available for source!'))

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
        coordsplit = coords.replace('+','-').split('-')
        if '.' not in coordsplit[0] and len(coordsplit[0]) > 6 and '.' not in coordsplit[1] and len(coordsplit[1]) > 6:
            newname = (prefix + 'J' + coordsplit[0][:6] + '.' + coordsplit[0][6:] +
                decsign + coordsplit[1][:6] + '.' + coordsplit[1][6:])
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
    elif (newname.startswith('SN') and is_number(newname[2:6]) and
        (len(newname) == 8 or len(newname) == 9) and newname[6:].isupper()):
        newname = 'SN' + newname[2:6] + newname[6:].lower()

    newname = (' '.join(newname.split())).strip()
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

def new_event(name, load = True, delete = True, loadifempty = True, refname = '',
    reference = '', url = '', bibcode = '', secondary = '', acknowledgment = ''):
    oldname = name
    name = add_event(name, load = load, delete = delete, loadifempty = loadifempty)
    source = add_source(name, refname = refname, reference = reference, url = url,
        bibcode = bibcode, secondary = secondary, acknowledgment = acknowledgment)
    add_quantity(name, 'alias', oldname, source)
    return (name, source)

def add_event(name, load = True, delete = True, loadifempty = True):
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
        events[newname]['schema'] = get_schema()
        events[newname]['name'] = newname
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

def snname(string):
    newstring = string.replace(' ', '').upper()
    if (newstring[:2] == "SN"):
        head = newstring[:6]
        tail = newstring[6:]
        if len(tail) >= 2 and tail[1] != '?':
            tail = tail.lower()
        newstring = head + tail

    return newstring

def add_source(name, refname = '', reference = '', url = '', bibcode = '', secondary = '', acknowledgment = ''):
    nsources = len(events[name]['sources']) if 'sources' in events[name] else 0
    if not refname:
        if not bibcode:
            raise(ValueError('Bibcode must be specified if name is not.'))

        if bibcode and len(bibcode) != 19:
            raise(ValueError('Bibcode "' + bibcode + '" must be exactly 19 characters long'))

        refname = bibcode

    if refname.upper().startswith('ATEL') and not bibcode:
        refname = refname.replace('ATEL', 'ATel').replace('Atel', 'ATel').replace('ATel #', 'ATel ').replace('ATel#', 'ATel').replace('ATel', 'ATel ')
        refname = ' '.join(refname.split())
        atelnum = refname.split()[-1]
        if is_number(atelnum) and atelnum in atelsdict:
            bibcode = atelsdict[atelnum]

    if refname.upper().startswith('CBET') and not bibcode:
        refname = refname.replace('CBET', 'CBET ')
        refname = ' '.join(refname.split())
        cbetnum = refname.split()[-1]
        if is_number(cbetnum) and cbetnum in cbetsdict:
            bibcode = cbetsdict[cbetnum]

    if refname.upper().startswith('IAUC') and not bibcode:
        refname = refname.replace('IAUC', 'IAUC ')
        refname = ' '.join(refname.split())
        iaucnum = refname.split()[-1]
        if is_number(iaucnum) and iaucnum in iaucsdict:
            bibcode = iaucsdict[iaucnum]

    for rep in sourcereps:
        if refname in sourcereps[rep]:
            refname = rep
            break

    for rep in urlreps:
        if url in urlreps[rep]:
            url = rep
            break

    if 'sources' not in events[name] or (refname not in [x['name'] for x in events[name]['sources']] and
        (not bibcode or bibcode not in [x['bibcode'] if 'bibcode' in x else '' for x in events[name]['sources']])):
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
        newsource['alias'] =  source
        if secondary:
            newsource['secondary'] = True
        events[name].setdefault('sources',[]).append(newsource)
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

def get_source_by_alias(name, alias):
    for source in events[name]['sources']:
        if source['alias'] == alias:
            return source
    print(events[name])
    raise(ValueError('Source alias not found!'))

def same_tag_str(photo, val, tag):
    issame = ((tag not in photo and not val) or (tag in photo and not val) or (tag in photo and photo[tag] == val))
    return issame

def same_tag_num(photo, val, tag, canbelist = False):
    issame = ((tag not in photo and not val) or (tag in photo and not val) or (tag in photo and
        ((not canbelist and Decimal(photo[tag]) == Decimal(val)) or
         (canbelist and
            ((isinstance(photo[tag], str) and isinstance(val, str) and Decimal(photo[tag]) == Decimal(val)) or
             (isinstance(photo[tag], list) and isinstance(val, list) and photo[tag] == val))))))
    return issame

def add_photometry(name, time = "", u_time = "MJD", e_time = "", telescope = "", instrument = "", band = "",
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

    if ((e_magnitude and not is_number(e_magnitude)) or (e_flux and not is_number(e_flux)) or
        (e_fluxdensity and not is_number(e_fluxdensity)) or (e_counts and not is_number(e_counts))):
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
        raise ValueError('Photometry must have source before being added!')

    if is_erroneous(name, 'photometry', source):
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
            if (same_tag_str(photo, sband, 'band') and
                same_tag_str(photo, u_time, 'u_time') and
                same_tag_num(photo, time, 'time', canbelist = True) and
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
                same_tag_str(photo, ssystem, 'system')
                ):
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
    events[name].setdefault('photometry',[]).append(photoentry)

def trim_str_arr(arr, length = 10):
    return [str(round_sig(float(x), length)) if (len(x) > length and len(str(round_sig(float(x), length))) < len(x)) else x for x in arr]

def add_spectrum(name, waveunit, fluxunit, wavelengths = "", fluxes = "", u_time = "", time = "", instrument = "",
    deredshifted = "", dereddened = "", errorunit = "", errors = "", source = "", snr = "", telescope = "",
    observer = "", reducer = "", survey = "", filename = "", observatory = "", data = ""):

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

    if not data and (not wavelengths or not fluxes):
        raise ValueError('Spectrum must have wavelengths and fluxes set, or data set.')

    if not source:
        raise ValueError('Spectrum must have source before being added!')

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
    events[name].setdefault('spectra',[]).append(spectrumentry)

def is_erroneous(name, field, sources):
    if 'errors' in events[name]:
        for alias in sources.split(','):
            source = get_source_by_alias(name, alias)
            if ('bibcode' in source and source['bibcode'] in
                [x['value'] for x in events[name]['errors'] if x['kind'] == 'bibcode' and x['extra'] == field]):
                    return True
            if ('name' in source and source['name'] in
                [x['value'] for x in events[name]['errors'] if x['kind'] == 'name' and x['extra'] == field]):
                    return True
    return False

def add_quantity(name, quantity, value, sources, forcereplacebetter = False, derived = '',
    lowerlimit = '', upperlimit = '', error = '', unit = '', kind = '', probability = '', extra = ''):
    if not quantity:
        raise(ValueError(name + "'s quantity must be specified for add_quantity."))
    if not sources:
        raise(ValueError(name + "'s source must be specified for quantity " + quantity + ' before it is added.'))
    if not isinstance(value, str) and (not isinstance(value, list) or not isinstance(value[0], str)):
        raise(ValueError(name + "'s Quantity " + quantity + " must be a string or an array of strings."))

    if is_erroneous(name, quantity, sources):
        return

    svalue = value.strip()
    serror = error.strip()
    skind = kind.strip()
    sprob = probability.strip()
    sunit = ''

    if not svalue or svalue == '--' or svalue == '-':
        return
    if serror and (not is_number(serror) or float(serror) < 0.):
        raise(ValueError(name + "'s quanta " + quantity + ' error value must be a number and positive.'))

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

        svalue = host_clean(svalue)
        
        if (not skind and ((svalue.lower().startswith('abell') and is_number(svalue[5:].strip())) or
            'cluster' in svalue.lower())):
            skind = 'cluster'
    elif quantity == 'claimedtype':
        isq = False
        svalue = svalue.replace('young', '')
        if svalue.lower() in ['unknown', 'unk', '?', '-']:
            return
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
        (svalue, sunit) = radec_clean(svalue, quantity, unit = unit)
    elif quantity == 'maxdate' or quantity == 'discoverdate':
        # Make sure month and day have leading zeroes
        sparts = svalue.split('/')
        if len(sparts[0]) > 4 and int(sparts[0]) > 0:
            raise ValueError('Date years limited to four digits.')
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
                        if sprob and 'probability' not in events[name][quantity][i]:
                            events[name][quantity][i]['probability'] = sprob
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
    if sprob:
        quantaentry['probability'] = sprob
    if sunit:
        quantaentry['unit'] = sunit
    if lowerlimit:
        quantaentry['lowerlimit'] = lowerlimit
    if upperlimit:
        quantaentry['upperlimit'] = upperlimit
    if derived:
        quantaentry['derived'] = derived
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

def load_cached_url(url, filepath, timeout = 120, write = True, failhard = False):
    filemd5 = ''
    filetxt = ''
    if not args.refresh and os.path.isfile(filepath):
        with codecs.open(filepath, 'r', encoding='utf8') as f:
            filetxt = f.read()
            if args.update:
                filemd5 = md5(filetxt.encode('utf-8')).hexdigest()

    try:
        session = requests.Session()
        response = session.get(url, timeout = timeout)
        response.raise_for_status()
        for x in response.history:
            x.raise_for_status()
            if x.status_code == 500 or x.status_code == 307 or x.status_code == 404:
                raise
        txt = response.text
        newmd5 = md5(txt.encode('utf-8')).hexdigest()
        if args.update and newmd5 == filemd5:
            tprint('Skipping file in "' + currenttask + '," local and remote copies identical [' + newmd5 + '].')
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
            source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
            add_quantity(name, 'maxdate', make_date_string(mldt.year, mldt.month, mldt.day), uniq_cdl([source]+mlsource.split(',')), derived = True)
        if mlmag:
            source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
            add_quantity(name, 'maxappmag', pretty_num(mlmag), uniq_cdl([source]+mlsource.split(',')), derived = True)
        if mlband:
            source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
            add_quantity(name, 'maxband', mlband, uniq_cdl([source]+mlsource.split(',')), derived = True)

    if 'discoverdate' not in events[name] or max([len(x['value'].split('/')) for x in events[name]['discoverdate']]) < 3:
        (fldt, flsource) = get_first_light(name)
        if fldt:
            source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
            add_quantity(name, 'discoverdate', make_date_string(fldt.year, fldt.month, fldt.day), uniq_cdl([source]+flsource.split(',')), derived = True)

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
            source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
            add_quantity(name, 'discoverdate', make_date_string(fldt.year, fldt.month, fldt.day), uniq_cdl([source]+minspecsource.split(',')), derived = True)

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
            bestsrc = z['source']

    return (bestz, bestkind, bestsig, bestsrc)

def jd_to_mjd(jd):
    return jd - Decimal(2400000.5)

def utf8(x):
    return str(x, 'utf-8')

def convert_aq_output(row):
    return OrderedDict([(x, str(row[x]) if is_number(row[x]) else row[x]) for x in row.colnames])

def set_preferred_names():
    currenttask = 'Setting preferred names'
    if not len(events):
        load_stubs()
    for ni, name in enumerate(tq(list(sorted(list(events.keys()))), currenttask = currenttask)):
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
            if load_event_from_file(newname):
                continue
            if load_event_from_file(name, delete = True):
                tprint('Changing event name (' + name + ') to preferred name (' + newname + ').')
                events[newname] = events[name]
                events[newname]['name'] = newname
                del(events[name])
                journal_events()
        if args.travis and ni > travislimit:
            break

# Merge and remove duplicate events
def merge_duplicates():
    if not len(events):
        load_stubs()
    currenttask = 'Merging duplicate events'
    keys = list(sorted(list(events.keys())))
    for n1, name1 in enumerate(tq(keys[:], currenttask)):
        if name1 not in events:
            continue
        allnames1 = set(get_aliases(name1) + (['AT' + name1[2:]] if (name1.startswith('SN') and is_number(name1[2:6])) else []))
        for name2 in keys[n1+1:]:
            if name2 not in events or name1 == name2:
                continue
            allnames2 = set(get_aliases(name2) + (['AT' + name2[2:]] if (name2.startswith('SN') and is_number(name2[2:6])) else []))
            if bool(allnames1 & allnames2):
                tprint('Found single event with multiple entries (' + name1 + ' and ' + name2 + '), merging.')
                load1 = load_event_from_file(name1, delete = True)
                load2 = load_event_from_file(name2, delete = True)
                if load1 and load2:
                    priority1 = 0
                    priority2 = 0
                    for an in list(allnames1):
                        if len(an) >= 2 and an.startswith(('SN', 'AT')):
                            priority1 = priority1 + 1
                    for an in list(allnames2):
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
        if args.travis and n1 > travislimit:
            break

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
                source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
                add_quantity(name, 'alias', name, source)

        if (name.startswith('SN') and is_number(name[2:6]) and 'discoverdate' in events[name] and
            int(events[name]['discoverdate'][0]['value'].split('/')[0]) >= 2016 and not any(['AT' in x for x in aliases])):
            source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
            add_quantity(name, 'alias', 'AT' + name[2:], source)

        events[name]['alias'] = list(sorted(events[name]['alias'], key=lambda key: alias_priority(name, key)))
        aliases = get_aliases(name)

        set_first_max_light(name)

        if 'claimedtype' in events[name]:
            events[name]['claimedtype'] = list(sorted(events[name]['claimedtype'], key=lambda key: ct_priority(name, key)))
        if 'discoverdate' not in events[name]:
            prefixes = ['MLS', 'SSS', 'CSS', 'GRB ']
            for alias in aliases:
                for prefix in prefixes:
                    if alias.startswith(prefix) and is_number(alias.replace(prefix, '')[:2]):
                        discoverdate = '/'.join(['20' + alias.replace(prefix, '')[:2],
                            alias.replace(prefix, '')[2:4], alias.replace(prefix, '')[4:6]])
                        if args.verbose:
                            tprint ('Added discoverdate from name [' + alias + ']: ' + discoverdate)
                        source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
                        add_quantity(name, 'discoverdate', discoverdate, source, derived = True)
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
                        source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
                        add_quantity(name, 'discoverdate', discoverdate, source, derived = True)
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
                        if args.verbose:
                            tprint ('Added discoverdate from name [' + alias + ']: ' + discoverdate)
                        source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
                        add_quantity(name, 'discoverdate', discoverdate, source, derived = True)
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
                        source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
                        add_quantity(name, 'discoverdate', discoverdate, source, derived = True)
                        break
                if 'discoverdate' in events[name]:
                    break
        if 'discoverdate' not in events[name]:
            prefixes = ['AT', 'SN', 'OGLE-', 'SM ', 'KSN-']
            for alias in aliases:
                for prefix in prefixes:
                    if alias.startswith(prefix) and is_number(alias.replace(prefix, '')[:4]) and '.' not in alias.replace(prefix, '')[:4]:
                        discoverdate = alias.replace(prefix, '')[:4]
                        if args.verbose:
                            tprint ('Added discoverdate from name [' + alias + ']: ' + discoverdate)
                        source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
                        add_quantity(name, 'discoverdate', discoverdate, source, derived = True)
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
                        noprefix = noprefix.replace('+','|').replace('-','|')
                        nops = noprefix.split('|')
                        if len(nops) < 2:
                            continue
                        rastr = nops[0]
                        decstr = nops[1]
                        ra = ':'.join([rastr[:2], rastr[2:4], rastr[4:6]]) + ('.' + rastr[6:] if len(rastr) > 6 else '') 
                        dec = decsign + ':'.join([decstr[:2], decstr[2:4], decstr[4:6]]) + ('.' + decstr[6:] if len(decstr) > 6 else '')
                        if args.verbose:
                            tprint ('Added ra/dec from name: ' + ra + ' ' + dec)
                        source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
                        add_quantity(name, 'ra', ra, source, derived = True)
                        add_quantity(name, 'dec', dec, source, derived = True)
                        break
                if 'ra' in events[name]:
                    break
        if ('ra' in events[name] and 'dec' in events[name] and 
            (not 'host' in events[name] or not any([x['value'] == 'Milky Way' for x in events[name]['host']]))):
            if name not in extinctionsdict:
                try:
                    result = IrsaDust.get_query_table(events[name]['ra'][0]['value'] + " " + events[name]['dec'][0]['value'], section = 'ebv')
                except (KeyboardInterrupt, SystemExit):
                    raise
                except:
                    warnings.warn("Coordinate lookup for " + name + " failed in IRSA.")
                else:
                    ebv = result['ext SandF mean'][0]
                    ebverr = result['ext SandF std'][0]
                    extinctionsdict[name] = [ebv, ebverr]
            if name in extinctionsdict:
                sources = uniq_cdl([add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True), 
                    add_source(name, bibcode = '2011ApJ...737..103S')])
                add_quantity(name, 'ebv', str(extinctionsdict[name][0]), sources, error = str(extinctionsdict[name][1]), derived = True)
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
                    source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
                    add_quantity(name, 'hostra', hostra, source, derived = True)
                    add_quantity(name, 'hostdec', hostdec, source, derived = True)
                    break
                if 'hostra' in events[name]:
                    break
        if 'claimedtype' in events[name]:
            events[name]['claimedtype'][:] = [ct for ct in events[name]['claimedtype'] if (ct['value'] != '?' and ct['value'] != '-')]
            if not len(events[name]['claimedtype']):
                del(events[name]['claimedtype'])
        if 'claimedtype' not in events[name] and name.startswith('AT'):
            source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
            add_quantity(name, 'claimedtype', 'Candidate', source)
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
                voc = float(besthv)*1.e5/clight
                source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
                sources = uniq_cdl([source] + bestsrc.split(','))
                add_quantity(name, 'redshift', pretty_num(sqrt((1. + voc)/(1. - voc)) - 1., sig = bestsig), sources,
                    kind = 'heliocentric', derived = True)
        if 'redshift' not in events[name] and has_task('nedd') and 'host' in events[name]:
            reference = "NED-D"
            refurl = "http://ned.ipac.caltech.edu/Library/Distances/"
            for host in events[name]['host']:
                if host['value'] in nedddict:
                    source = add_source(name, bibcode = '2015arXiv150201589P')
                    secondarysource = add_source(name, refname = reference, url = refurl, secondary = True)
                    meddist = statistics.median(nedddict[host['value']])
                    redshift = pretty_num(z_at_value(cosmo.comoving_distance, float(meddist) * un.Mpc), sig = get_sig_digits(str(meddist)))
                    add_quantity(name, 'redshift', redshift, uniq_cdl([source,secondarysource]), kind = 'host', derived = True)
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
                source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
                sources = uniq_cdl([source] + bestsrc.split(','))
                add_quantity(name, 'maxabsmag', pretty_num(float(events[name]['maxappmag'][0]['value']) -
                    5.0*(log10(float(bestld)*1.0e6) - 1.0), sig = bestsig), sources, derived = True)
        if 'redshift' in events[name]:
            # Find the "best" redshift to use for this
            (bestz, bestkind, bestsig, bestsrc) = get_best_redshift(name)
            if bestsig > 0:
                bestz = float(bestz)
                if 'velocity' not in events[name]:
                    source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
                    sources = uniq_cdl([source] + bestsrc.split(','))
                    add_quantity(name, 'velocity', pretty_num(clight/km*((bestz + 1.)**2. - 1.)/
                        ((bestz + 1.)**2. + 1.), sig = bestsig), sources, kind = prefkinds[bestkind], derived = True)
                if bestz > 0.:
                    if 'lumdist' not in events[name]:
                        dl = cosmo.luminosity_distance(bestz)
                        sources = [add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True),
                            add_source(name, bibcode = '2015arXiv150201589P')]
                        sources = uniq_cdl(sources + bestsrc.split(','))
                        add_quantity(name, 'lumdist', pretty_num(dl.value, sig = bestsig), sources,
                            kind = prefkinds[bestkind], derived = True)
                        if 'maxabsmag' not in events[name] and 'maxappmag' in events[name]:
                            source = add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
                            add_quantity(name, 'maxabsmag', pretty_num(float(events[name]['maxappmag'][0]['value']) -
                                5.0*(log10(dl.to('pc').value) - 1.0), sig = bestsig), sources, derived = True)
                    if 'comovingdist' not in events[name]:
                        cd = cosmo.comoving_distance(bestz)
                        sources = [add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True),
                            add_source(name, bibcode = '2015arXiv150201589P')]
                        sources = uniq_cdl(sources + bestsrc.split(','))
                        add_quantity(name, 'comovingdist', pretty_num(cd.value, sig = bestsig), sources, derived = True)
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
                sources = uniq_cdl([add_source(name, bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)] +
                    events[name]['ra'][0]['source'].split(',') + events[name]['dec'][0]['source'].split(',') +
                    events[name]['hostra'][0]['source'].split(',') + events[name]['hostdec'][0]['source'].split(','))
                if 'hostoffsetang' not in events[name]:
                    add_quantity(name, 'hostoffsetang', pretty_num(Decimal(hypot(c1.ra.degree - c2.ra.degree,
                        c1.dec.degree - c2.dec.degree))*Decimal(3600.)), sources, derived = True, unit = 'arcseconds')
                if 'comovingdist' in events[name] and 'redshift' in events[name] and 'hostoffsetdist' not in events[name]:
                    offsetsig = get_sig_digits(events[name]['hostoffsetang'][0]['value'])
                    sources = uniq_cdl(sources.split(',') +
                        events[name]['comovingdist'][0]['source'].split(',') + events[name]['redshift'][0]['source'].split(','))
                    add_quantity(name, 'hostoffsetdist',
                        pretty_num(float(events[name]['hostoffsetang'][0]['value']) / 3600. * (pi / 180.) *
                        float(events[name]['comovingdist'][0]['value']) * 1000. / (1.0 + float(events[name]['redshift'][0]['value'])),
                        sig = offsetsig), sources)
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
                    source['reference'] = bibauthordict[source['bibcode']]
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

def write_all_events(empty = False, gz = False, bury = False):
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

        outdir = "../" + str(repofolders[0])
        if 'discoverdate' in events[name]:
            for r, year in enumerate(repoyears):
                eyr = events[name]['discoverdate'][0]['value'].split('/')[0]
                if is_integer(eyr) and int(eyr) <= year:
                    outdir = "../" + str(repofolders[r])
                    break

        # Bury non-SN events here if only claimed type is non-SN type, or if primary
        # name starts with a non-SN prefix.
        if bury:
            buryevent = False
            nonsneprefixes = ('PNVJ', 'PNV J', 'OGLE-2013-NOVA', 'EV*', 'V*', "Nova")
            if name.startswith(nonsneprefixes):
                tprint('Burying ' + name + ', non-SNe prefix.')
                continue
            if 'claimedtype' in events[name]:
                for ct in events[name]['claimedtype']:
                    if ct['value'].upper() not in nonsnetypes and ct['value'].upper() != 'CANDIDATE':
                        buryevent = False
                        break
                    if ct['value'].upper() in nonsnetypes:
                        buryevent = True
                if buryevent:
                    tprint('Burying ' + name + ' (' + ct['value'] + ').')
                    outdir = '../sne-boneyard'

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
            newsourcealiases[source['alias']] = (add_source(destname,
                bibcode = source['bibcode'] if 'bibcode' in source else '',
                refname = source['name'] if 'name' in source else '',
                reference = source['reference'] if 'reference' in source else '',
                url = source['url'] if 'url' in source else ''))

    if 'errors' in events[fromname]:
        for err in events[fromname]['errors']:
            events[destname].setdefault('errors',[]).append(err)

    for key in keys:
        if key not in ['schema', 'name', 'sources', 'errors']:
            for item in events[fromname][key]:
                isd = False
                sources = []
                if 'source' not in item:
                    raise ValueError("Item has no source!")
                for sid in item['source'].split(','):
                    if sid == 'D':
                        sources.append('D')
                    elif sid in newsourcealiases:
                        sources.append(newsourcealiases[sid])
                    else:
                        raise ValueError("Couldn't find source alias!")
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
                        unit = null_field(item, "unit"), probability = null_field(item, "probability"), kind = null_field(item, "kind"))

def load_event_from_file(name = '', location = '', clean = False, delete = True):
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
        if path or namepath:
            if name in events:
                del events[name]

        if path:
            with open(path, 'r') as f:
                newevent = json.loads(f.read(), object_pairs_hook=OrderedDict)
        elif namepath:
            with open(namepath, 'r') as f:
                newevent = json.loads(f.read(), object_pairs_hook=OrderedDict)

        if newevent:
            if clean:
                newevent = clean_event(newevent)
            name = next(reversed(newevent))

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

    if 'schema' not in events['temp']:
        events['temp']['schema'] = get_schema()
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
                add_source('temp', refname = source['name'], url = source['url'])

    # Clean some legacy fields
    if 'aliases' in events['temp'] and isinstance(events['temp']['aliases'], list):
        source = add_source('temp', bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
        for alias in events['temp']['aliases']:
            add_quantity('temp', 'alias', alias, source)
        del(events['temp']['aliases'])

    if ('distinctfrom' in events['temp'] and isinstance(events['temp']['distinctfrom'], list) and
        isinstance(events['temp']['distinctfrom'][0], str)):
        distinctfroms = [x for x in events['temp']['distinctfrom']]
        del(events['temp']['distinctfrom'])
        source = add_source('temp', bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
        for df in distinctfroms:
            add_quantity('temp', 'distinctfrom', df, source)

    if ('errors' in events['temp'] and isinstance(events['temp']['errors'], list) and
        'sourcekind' in events['temp']['errors'][0]):
        source = add_source('temp', bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
        for err in events['temp']['errors']:
            add_quantity('temp', 'error', err['quantity'], source, kind = err['sourcekind'], extra = err['id'])
        del(events['temp']['errors'])

    if not bibcodes:
        add_source('temp', bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
        bibcodes = [oscbibcode]

    for key in list(events['temp'].keys()):
        if key in ['name', 'schema', 'sources', 'errors']:
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

def archived_task(task):
    if 'archived' in tasks[task] and args.archived and task not in args.refreshlist.split(','):
        return True
    if ('archived' in tasks[task] and tasks[task]['archived'] and
        task not in args.refreshlist.split(',') and not args.fullrefresh):
        return True
    return False

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
    global currenttask
    currenttask = 'Loading event stubs'

    files = repo_file_list()

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
    for fi in tq(files, currenttask):
        fname = fi
        if '.gz' in fi:
            fname = fi.replace('.gz', '')
            with gzip.open(fi, 'rb') as f_in, open(fname, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(fi)
        name = os.path.basename(os.path.splitext(fname)[0]).replace('.json', '')
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
        currenttask = 'Deleting old events'
        delete_old_event_files()

    # Import data provided directly to OSC
    if do_task(task, 'internal'):
        for datafile in tq(sorted(glob("../sne-internal/*.json"), key=lambda s: s.lower()), currenttask):
            if not load_event_from_file(location = datafile, clean = True, delete = False):
                raise IOError('Failed to find specified file.')
        journal_events()
    
    if do_task(task, 'radio'):
        for datafile in tq(sorted(glob("../sne-external-radio/*.txt"), key=lambda s: s.lower()), currenttask):
            oldname = os.path.basename(datafile).split('.')[0]
            name = add_event(oldname)
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
                        add_quantity(name, 'alias', oldname, source)
        journal_events()
    
    if do_task(task, 'xray'):
        for datafile in tq(sorted(glob("../sne-external-xray/*.txt"), key=lambda s: s.lower()), currenttask):
            oldname = os.path.basename(datafile).split('.')[0]
            name = add_event(oldname)
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
                        add_quantity(name, 'alias', oldname, source)
        journal_events()
    
    if do_task(task, 'simbad'):
        #Simbad.list_votable_fields()
        # Some coordinates that SIMBAD claims belong to the SNe actually belong to the host.
        simbadmirrors = ['http://simbad.harvard.edu/simbad/sim-script', 'http://simbad.u-strasbg.fr/simbad/sim-script']
        simbadbadcoordbib = ['2013ApJ...770..107C']
        simbadbadnamebib = ['2004AJ....127.2809W', '2005MNRAS.364.1419Z', '2015A&A...574A.112D', '2011MNRAS.417..916G', '2002ApJ...566..880G']
        simbadbannedcats = ['[TBV2008]', 'OGLE-MBR']
        customSimbad = Simbad()
        customSimbad.ROW_LIMIT = -1
        customSimbad.TIMEOUT = 120
        customSimbad.add_votable_fields('otype', 'sptype', 'sp_bibcode', 'id')
        for mirror in simbadmirrors:
            customSimbad.SIMBAD_URL = mirror
            try:
                table = customSimbad.query_criteria('maintype=SN | maintype="SN?"')
            except:
                continue
            else:
                break

        # 2000A&AS..143....9W
        for brow in tq(table, currenttask = currenttask):
            row = {x:re.sub(r'b\'(.*)\'', r'\1', str(brow[x])) for x in brow.colnames}
            # Skip items with no bibliographic info aside from SIMBAD, too error-prone
            if row['OTYPE'] == 'Candidate_SN*' and not row['SP_TYPE']:
                continue
            if not row['COO_BIBCODE'] and not row['SP_BIBCODE'] and not row['SP_BIBCODE_2']:
                continue
            if any([x in row['MAIN_ID'] for x in simbadbannedcats]):
                continue
            if row['COO_BIBCODE'] and row['COO_BIBCODE'] in simbadbadnamebib:
                continue
            name = single_spaces(re.sub(r'\[[^)]*\]', '', row['MAIN_ID']).strip())
            if name == 'SN':
                continue
            if is_number(name):
                continue
            name = add_event(name)
            source = add_source(name, refname = 'SIMBAD astronomical database', bibcode = "2000A&AS..143....9W",
                url = "http://simbad.u-strasbg.fr/", secondary = True)
            aliases = row['ID'].split(',')
            for alias in aliases:
                if any([x in alias for x in simbadbannedcats]):
                    continue
                ali = single_spaces(re.sub(r'\[[^)]*\]', '', alias).strip())
                if is_number(ali):
                    continue
                ali = name_clean(ali)
                add_quantity(name, 'alias', ali, source)
            if row['COO_BIBCODE'] and row['COO_BIBCODE'] not in simbadbadcoordbib:
                csources = ','.join([source, add_source(name, bibcode = row['COO_BIBCODE'])])
                add_quantity(name, 'ra', row['RA'], csources)
                add_quantity(name, 'dec', row['DEC'], csources)
            if row['SP_BIBCODE']:
                ssources = uniq_cdl([source, add_source(name, bibcode = row['SP_BIBCODE'])] +
                    ([add_source(name, bibcode = row['SP_BIBCODE_2'])] if row['SP_BIBCODE_2'] else []))
                add_quantity(name, 'claimedtype', row['SP_TYPE'].replace('SN.', '').replace('SN', '').replace('(~)', '').strip(': '), ssources)
        journal_events()
    
    # Import primary data sources from Vizier
    if do_task(task, 'vizier'):
        Vizier.ROW_LIMIT = -1
    
        # 2012ApJS..200...12H
        result = Vizier.get_catalogs("J/ApJS/200/12/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        oldname = ''
        for row in tq(table, currenttask):
            name = row['SN']
            if is_number(name[:4]):
                name = 'SN' + name
            (name, source) = new_event(name, bibcode = "2012ApJS..200...12H")
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
        for row in tq(table, currenttask):
            name = row['Name'].replace('SCP', 'SCP-')
            (name, source) = new_event(name, bibcode = "2012ApJ...746...85S")
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
        for row in tq(table, currenttask):
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
            (name, source) = new_event(name, bibcode = "2012ApJ...746...85S")
            add_photometry(name, time = str(row['MJD']), band = row['Filter'], instrument = row['Inst'],
                magnitude = magnitude, e_magnitude = e_magnitude, source = source)
    
        # 2004ApJ...602..571B
        result = Vizier.get_catalogs("J/ApJ/602/571/table8")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        oldname = ''
        for row in tq(table, currenttask):
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
            (name, source) = new_event(name, bibcode = "2004ApJ...602..571B")
            band = row['Filt']
            system = ''
            telescope = ''
            if band in ['R', 'I']:
                system = 'Cousins'
            if band == 'Z':
                telescope = 'Subaru'
            add_photometry(name, time = str(row['MJD']), band = band, system = system, telescope = telescope,
                magnitude = magnitude, e_magnitude = e_magnitude, source = source)
        
        # 2014MNRAS.444.3258M
        result = Vizier.get_catalogs("J/MNRAS/444/3258/SNe")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        oldname = ''
        for row in tq(table, currenttask):
            name = row['SN']
            if name == oldname:
                continue
            oldname = name
            (name, source) = new_event(name, bibcode = '2014MNRAS.444.3258M')
            add_quantity(name, 'redshift', str(row['z']), source, kind = 'heliocentric', error = str(row['e_z']))
            add_quantity(name, 'ra', str(row['_RA']), source, unit = 'floatdegrees')
            add_quantity(name, 'dec', str(row['_DE']), source, unit = 'floatdegrees')
        journal_events()
    
        # 2014MNRAS.438.1391P
        result = Vizier.get_catalogs("J/MNRAS/438/1391/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            name = row['SN']
            (name, source) = new_event(name, bibcode = '2014MNRAS.438.1391P')
            add_quantity(name, 'redshift', str(row['zh']), source, kind = 'heliocentric')
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
        journal_events()
    
        # 2012ApJ...749...18B
        result = Vizier.get_catalogs("J/ApJ/749/18/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            name = row['Name'].replace(' ','')
            (name, source) = new_event(name, bibcode = '2012ApJ...749...18B')
            mjd = str(astrotime(2450000.+row['JD'], format='jd').mjd)
            band = row['Filt']
            magnitude = str(row['mag'])
            e_magnitude = str(row['e_mag'])
            e_magnitude = '' if e_magnitude == '--' else e_magnitude
            upperlimit = True if row['l_mag'] == '>' else False
            add_photometry(name, time = mjd, band = band, magnitude = magnitude, e_magnitude = e_magnitude, instrument = 'UVOT',
                source = source, upperlimit = upperlimit, telescope = 'Swift', system = 'Swift')
        journal_events()
    
        # 2010A&A...523A...7G
        result = Vizier.get_catalogs("J/A+A/523/A7/table9")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            name = 'SNLS-' + row['SNLS']
            (name, source) = new_event(name, bibcode = '2010A&A...523A...7G')
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
        for row in tq(table, currenttask):
            name = 'SN' + row['SN']
            (name, source) = new_event(name, bibcode = '2004A&A...415..863G')
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
        for row in tq(table, currenttask):
            name = 'SDSS-II SN ' + str(row['SNID'])
            (name, source) = new_event(name, bibcode = '2008AJ....136.2306H')
            add_quantity(name, 'claimedtype', row['SpType'].replace('SN.', '').strip(':'), source)
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
    
        # 2010ApJ...708..661D
        result = Vizier.get_catalogs("J/ApJ/708/661/sn")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            name = row['SN']
            if not name:
                name = 'SDSS-II SN ' + str(row['SDSS-II'])
            else:
                name = 'SN' + name
            (name, source) = new_event(name, bibcode = '2010ApJ...708..661D')
            add_quantity(name, 'alias', 'SDSS-II SN ' + str(row['SDSS-II']), source)
            add_quantity(name, 'claimedtype', 'II P', source)
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
    
        result = Vizier.get_catalogs("J/ApJ/708/661/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            if row['f_SN'] == 'a':
                name = 'SDSS-II SN ' + str(row['SN'])
            else:
                name = 'SN' + row['SN']
            (name, source) = new_event(name, bibcode = '2010ApJ...708..661D')
            add_quantity(name, 'redshift', str(row['z']), source, error = str(row['e_z']))
        journal_events()
    
        # 2014ApJ...795...44R
        result = Vizier.get_catalogs("J/ApJ/795/44/ps1_snIa")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            name = row['SN']
            (name, source) = new_event(name, bibcode = '2014ApJ...795...44R')
            astrot = astrotime(row['tdisc'], format='mjd').datetime
            add_quantity(name, 'discoverdate',  make_date_string(astrot.year, astrot.month, astrot.day), source)
            add_quantity(name, 'redshift', str(row['z']), source, error = str(row['e_z']), kind = 'heliocentric')
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
            add_quantity(name, 'claimedtype', 'Ia', source)
    
        result = Vizier.get_catalogs("J/ApJ/795/44/table6")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            name = row['SN']
            (name, source) = new_event(name, bibcode = '2014ApJ...795...44R')
            if row['mag'] != '--':
                add_photometry(name, time = str(row['MJD']), band = row['Filt'], magnitude = str(row['mag']),
                    e_magnitude = str(row['e_mag']), source = source, system = 'AB', telescope = 'PS1', instrument = 'PS1')
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
    
        for row in tq(table, currenttask):
            if row['band'][0] == '(':
                continue
            oldname = 'SN' + row['SN']
            name = add_event(oldname)
            source = ''
            secsource = add_source(name, bibcode = '1990A&AS...82..145C', secondary = True)
            mjd = str(jd_to_mjd(Decimal(row['JD'])))
            mag = str(row['m'])
            band = row['band'].strip("'")
            if row['r_m'] in ii189bibdict:
                source = add_source(name, bibcode = ii189bibdict[row['r_m']])
            else:
                source = add_source(name, refname = ii189refdict[row['r_m']])
            add_quantity(name, 'alias', oldname, source)
    
            add_photometry(name, time = mjd, band = band, magnitude = mag, source = uniq_cdl([source,secsource]))
        journal_events()
    
        # 2014yCat.7272....0G
        result = Vizier.get_catalogs("VII/272/snrs")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
    
        for row in tq(table, currenttask):
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
    
            oldname = name
            name = add_event(name)
            source = (add_source(name, bibcode = '2014BASI...42...47G') + ',' +
                      add_source(name, refname = 'Galactic SNRs', url = 'https://www.mrao.cam.ac.uk/surveys/snrs/snrs.data.html'))
            add_quantity(name, 'alias', oldname, source)
    
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
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = 'SN' + row['SN']
            (name, source) = new_event(name, bibcode = '2014MNRAS.442..844F')
            add_quantity(name, 'redshift', str(row['zhost']), source, kind = 'host')
            add_quantity(name, 'ebv', str(row['E_B-V_']), source)
        journal_events()
    
        result = Vizier.get_catalogs("J/MNRAS/442/844/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = 'SN' + str(row['SN'])
            (name, source) = new_event(name, bibcode = "2014MNRAS.442..844F")
            for band in ['B', 'V', 'R', 'I']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = row['MJD'], band = band, magnitude = row[bandtag],
                        e_magnitude = row['e_' + bandtag], source = source, telescope = 'KAIT', instrument = 'KAIT')
        journal_events()
    
        # 2012MNRAS.425.1789S
        result = Vizier.get_catalogs("J/MNRAS/425/1789/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = ''.join(row['SimbadName'].split(' '))
            (name, source) = new_event(name, bibcode = '2012MNRAS.425.1789S')
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
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = u'LSQ' + str(row['LSQ'])
            (name, source) = new_event(name, bibcode = "2015ApJS..219...13W")
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
            add_quantity(name, 'redshift', row['z'], source, error = row['e_z'], kind = 'heliocentric')
            add_quantity(name, 'ebv', row['E_B-V_'], source)
            add_quantity(name, 'claimedtype', 'Ia', source)
        result = Vizier.get_catalogs("J/ApJS/219/13/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = 'LSQ' + row['LSQ']
            (name, source) = new_event(name, bibcode = "2015ApJS..219...13W")
            add_photometry(name, time = str(jd_to_mjd(Decimal(row['JD']))), instrument = 'QUEST', telescope = 'ESO Schmidt',
                observatory = 'La Silla', band = row['Filt'],
                magnitude = row['mag'], e_magnitude = row['e_mag'], system = "Swope", source = source)
        journal_events()
    
        # 2012Natur.491..228C
        result = Vizier.get_catalogs("J/other/Nat/491.228/tablef1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        name = 'SN2213-1745'
        (name, source) = new_event(name, bibcode = "2012Natur.491..228C")
        add_quantity(name, 'claimedtype', 'SLSN-R', source)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            for band in ['g', 'r', 'i']:
                bandtag = band + '_mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = row["MJD" + band + "_"], band = band + "'", magnitude = row[bandtag],
                        e_magnitude = row["e_" + bandtag], source = source)
    
        result = Vizier.get_catalogs("J/other/Nat/491.228/tablef2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        name = 'SN1000+0216'
        (name, source) = new_event(name, bibcode = "2012Natur.491..228C")
        add_quantity(name, 'claimedtype', 'SLSN-II?', source)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            for band in ['g', 'r', 'i']:
                bandtag = band + '_mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = row["MJD" + band + "_"], band = band + "'", magnitude = row[bandtag],
                        e_magnitude = row["e_" + bandtag], source = source)
        journal_events()
    
        # 2011Natur.474..484Q
        result = Vizier.get_catalogs("J/other/Nat/474.484/tables1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = str(row['Name'])
            (name, source) = new_event(name, bibcode = "2011Natur.474..484Q")
            add_photometry(name, time = row['MJD'], band = row['Filt'], telescope = row['Tel'],
                magnitude = row['mag'], e_magnitude = row['e_mag'], source = source)
        journal_events()
    
        # 2011ApJ...736..159G
        result = Vizier.get_catalogs("J/ApJ/736/159/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        name = 'PTF10vdl'
        (name, source) = new_event(name, bibcode = "2011ApJ...736..159G")
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            add_photometry(name, time = str(jd_to_mjd(Decimal(row['JD']))), band = row['Filt'], telescope = row['Tel'], magnitude = row['mag'],
                           e_magnitude = row['e_mag'] if is_number(row['e_mag']) else '', upperlimit = (not is_number(row['e_mag'])), source = source)
        journal_events()
    
        # 2012ApJ...760L..33B
        result = Vizier.get_catalogs("J/ApJ/760/L33/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        name = 'PTF12gzk'
        (name, source) = new_event(name, bibcode = "2012ApJ...760L..33B")
        for row in tq(table, currenttask):
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
        (name, source) = new_event(name, bibcode = "2013ApJ...769...39S")
        for row in tq(table, currenttask):
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
        (name, source) = new_event(name, bibcode = "2009MNRAS.394.2266P")
        result = Vizier.get_catalogs("J/MNRAS/394/2266/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            for band in ['U', 'B', 'V', 'R', 'I']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), band = band, magnitude = row[bandtag],
                        e_magnitude = (row["e_" + bandtag] if row['l_' + bandtag] != '>' else ''),
                        source = source, upperlimit = (row['l_' + bandtag] == '>'))
            if "zmag" in row and is_number(row["zmag"]) and not isnan(float(row["zmag"])):
                add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), band = "z", magnitude = row["zmag"],
                               e_magnitude = row["e_zmag"], source = source)
    
        result = Vizier.get_catalogs("J/MNRAS/394/2266/table3")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            for band in ['B', 'V', 'R']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), band = band, magnitude = row[bandtag],
                        e_magnitude = (row["e_" + bandtag] if row['l_' + bandtag] != '>' else ''),
                        source = source, upperlimit = (row['l_' + bandtag] == '>'))
    
        result = Vizier.get_catalogs("J/MNRAS/394/2266/table4")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            for band in ['J', 'H', 'K']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), band = band, magnitude = row[bandtag],
                                   e_magnitude = row["e_" + bandtag], source = source)
        journal_events()
    
        # 2013AJ....145...99A
        result = Vizier.get_catalogs("J/AJ/145/99/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        name = 'SN2003ie'
        (name, source) = new_event(name, bibcode = "2013AJ....145...99A")
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            for band in ['B', 'R']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = row["MJD"], band = band, magnitude = row[bandtag],
                                   e_magnitude = row["e_" + bandtag] if not row["l_" + bandtag] else '',
                                   upperlimit = (row['l_' + bandtag] == '>'), source = source)
            for band in ['V', 'I']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = row["MJD"], band = band, magnitude = row[bandtag],
                                   e_magnitude = row["e_" + bandtag] if is_number(row["e_" + bandtag]) else '',
                                   upperlimit = (not is_number(row["e_" + bandtag])), source = source)
        journal_events()
    
        # 2011ApJ...729..143C
        name = 'SN2008am'
        (name, source) = new_event(name, bibcode = "2011ApJ...729..143C")
    
        result = Vizier.get_catalogs("J/ApJ/729/143/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            add_photometry(name, time = row['MJD'], band = 'ROTSE', telescope = 'ROTSE', magnitude = row['mag'],
                           e_magnitude = row['e_mag'] if not row['l_mag'] else '', upperlimit = (row['l_mag'] == '<'), source = source)
    
        result = Vizier.get_catalogs("J/ApJ/729/143/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            for band in ['J', 'H', 'Ks']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = row["MJD"], telescope = "PAIRITEL", band = band, magnitude = row[bandtag],
                                   e_magnitude = row["e_" + bandtag], source = source)
    
        result = Vizier.get_catalogs("J/ApJ/729/143/table4")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            add_photometry(name, time = row['MJD'], band = row['Filt'], telescope = 'P60', magnitude = row['mag'],
                           e_magnitude = row['e_mag'], source = source)
    
        result = Vizier.get_catalogs("J/ApJ/729/143/table5")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            add_photometry(name, time = row['MJD'], band = row['Filt'], instrument = 'UVOT', telescope = 'Swift', magnitude = row['mag'],
                           e_magnitude = row['e_mag'], source = source)
        journal_events()
    
        # 2011ApJ...728...14P
        name = 'SN2009bb'
        (name, source) = new_event(name, bibcode = "2011ApJ...728...14P")
    
        result = Vizier.get_catalogs("J/ApJ/728/14/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            for band in ['B', 'V', 'R', 'I']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), telescope = row["Tel"], band = band, magnitude = row[bandtag],
                                   e_magnitude = row["e_" + bandtag], source = source)
    
        result = Vizier.get_catalogs("J/ApJ/728/14/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            for band in ['u', 'g', 'r', 'i', 'z']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), telescope = row["Tel"], band = band + "'", magnitude = row[bandtag],
                                   e_magnitude = row["e_" + bandtag], source = source)
    
        result = Vizier.get_catalogs("J/ApJ/728/14/table3")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            for band in ['Y', 'J', 'H']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = str(jd_to_mjd(Decimal(row["JD"]))), instrument = row['Inst'], band = band, magnitude = row[bandtag],
                                   e_magnitude = row["e_" + bandtag], source = source)
        journal_events()
    
        # 2011PAZh...37..837T
        name = 'SN2009nr'
        (name, source) = new_event(name, bibcode = "2011PAZh...37..837T")
    
        result = Vizier.get_catalogs("J/PAZh/37/837/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            mjd = str(jd_to_mjd(Decimal(row["JD"]) + 2455000))
            for band in ['U', 'B', 'V', 'R', 'I']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = mjd, telescope = row["Tel"], band = band, magnitude = row[bandtag],
                                   e_magnitude = row["e_" + bandtag], source = source)
        journal_events()
    
        # 2013MNRAS.433.1871B
        name = 'SN2012aw'
        (name, source) = new_event(name, bibcode = "2013MNRAS.433.1871B")
    
        result = Vizier.get_catalogs("J/MNRAS/433/1871/table3a")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            mjd = str(jd_to_mjd(Decimal(row["JD"]) + 2456000))
            for band in ['U', 'B', 'V', 'Rc', 'Ic']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = mjd, telescope = row["Tel"], band = band, magnitude = row[bandtag],
                                   e_magnitude = row["e_" + bandtag], source = source)
    
        result = Vizier.get_catalogs("J/MNRAS/433/1871/table3b")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            mjd = str(jd_to_mjd(Decimal(row["JD"]) + 2456000))
            for band in ['g', 'r', 'i', 'z']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = mjd, telescope = row["Tel"], band = band, magnitude = row[bandtag],
                                   e_magnitude = row["e_" + bandtag], source = source)
        journal_events()
    
        # 2014AJ....148....1Z
        name = 'SN2012fr'
        (name, source) = new_event(name, bibcode = "2014AJ....148....1Z")
    
        result = Vizier.get_catalogs("J/AJ/148/1/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            mjd = row['MJD']
            for band in ['B', 'V', 'R', 'I']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = mjd, telescope = "LJT", instrument = "YFOSC", band = band, magnitude = row[bandtag],
                                   e_magnitude = row["e_" + bandtag], source = source)
    
        result = Vizier.get_catalogs("J/AJ/148/1/table3")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            mjd = row['MJD']
            for band in ['U', 'B', 'V', 'UVW1', 'UVW2', 'UVM2']:
                bandtag = band + 'mag' if len(band) == 1 else band
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = mjd, telescope = "Swift", instrument = "UVOT", band = band, magnitude = row[bandtag],
                                   e_magnitude = row["e_" + bandtag], source = source)
    
        result = Vizier.get_catalogs("J/AJ/148/1/table5")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            mjd = row['MJD']
            for band in ['B', 'V', 'R', 'I']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = mjd, telescope = "LJT", band = band, magnitude = row[bandtag],
                                   e_magnitude = row["e_" + bandtag], source = source)
        journal_events()
    
        # 2015ApJ...805...74B
        name = 'SN2014J'
        (name, source) = new_event(name, bibcode = "2014AJ....148....1Z")
    
        result = Vizier.get_catalogs("J/ApJ/805/74/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
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
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = str(row['SN'])
            (name, source) = new_event(name, bibcode = "2011ApJ...741...97D")
            add_photometry(name, time = str(jd_to_mjd(Decimal(row['JD']))), band = row['Filt'], magnitude = row['mag'],
                           e_magnitude = row['e_mag'] if is_number(row['e_mag']) else '', upperlimit = (not is_number(row['e_mag'])), source = source)
        journal_events()
    
        # 2015MNRAS.448.1206M
        # Note: Photometry from two SN can also be added from this source.
        result = Vizier.get_catalogs("J/MNRAS/448/1206/table3")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = str(row['Name'])
            (name, source) = new_event(name, bibcode = "2015MNRAS.448.1206M")
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
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = str(row['Name'])
            (name, source) = new_event(name, bibcode = "2015MNRAS.448.1206M")
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
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = str(row['Name'])
            (name, source) = new_event(name, bibcode = "2015MNRAS.448.1206M")
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
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = str(row['Name'])
            (name, source) = new_event(name, bibcode = "2015MNRAS.448.1206M")
            add_quantity(name, 'discoverdate', '20' + name[4:6], source)
            add_quantity(name, 'ra', row['RAJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'dec', row['DEJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'maxappmag', row['rP1mag'], source, error = row['e_rP1mag'])
            add_quantity(name, 'maxband', 'r', source)
            add_quantity(name, 'claimedtype', row['Type'], source)
        result = Vizier.get_catalogs("J/MNRAS/448/1206/tablea2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = str(row['Name'])
            (name, source) = new_event(name, bibcode = "2015MNRAS.448.1206M")
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
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = str(row['Name'])
            (name, source) = new_event(name, bibcode = "2015MNRAS.448.1206M")
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
        for row in tq(table, currenttask):
            if not row['Wcl'] or row['Wcl'] == 'N':
                continue
            row = convert_aq_output(row)
            name = str(row['SN']).replace(' ', '')
            (name, source) = new_event(name, bibcode = "2012AJ....143..126B")
            add_quantity(name, 'claimedtype', 'Ia-' + row['Wcl'], source)
        journal_events()

        # 2015ApJS..220....9F
        for viztab in ['1', '2']:
            result = Vizier.get_catalogs("J/ApJS/220/9/table" + viztab)
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for row in tq(table, currenttask):
                row = convert_aq_output(row)
                (name, source) = new_event(row['SN'], bibcode = "2015ApJS..220....9F")
                add_quantity(name, 'claimedtype', row['Type'], source)
                add_quantity(name, 'ra', row['RAJ2000'], source, unit = 'floatdegrees')
                add_quantity(name, 'dec', row['DEJ2000'], source, unit = 'floatdegrees')
                if '?' not in row['Host']:
                    add_quantity(name, 'host', row['Host'].replace('_', ' '), source)
                kind = ''
                if 'Host' in row['n_z']:
                    kind = 'host'
                elif 'Spectrum' in row['n_z']:
                    kind = 'spectroscopic'
                add_quantity(name, 'redshift', row['z'], source, error = row['e_z'], kind = kind)
        result = Vizier.get_catalogs("J/ApJS/220/9/table8")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            (name, source) = new_event(row['SN'], bibcode = "2015ApJS..220....9F")
            add_quantity(name, 'claimedtype', row['Type'], source)
            add_photometry(name, time = row['MJD'], band = row['Band'], magnitude = row['mag'],
                e_magnitude = row["e_mag"], telescope = row["Tel"], source = source)
        journal_events()

        # 2008ApJ...673..999P
        result = Vizier.get_catalogs("J/ApJ/673/999/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            (name, source) = new_event('SN'+row['SN'], bibcode = "2008ApJ...673..999P")
            add_quantity(name, 'ra', row['RAJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'dec', row['DEJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'redshift', row['z'], source, kind = 'host')
            add_quantity(name, 'hostra', row['RAGdeg'], source, unit = 'floatdegrees')
            add_quantity(name, 'hostdec', row['DEGdeg'], source, unit = 'floatdegrees')
            add_quantity(name, 'claimedtype', row['Type'].strip(':'), source)
        journal_events()

        # 2011MNRAS.417..916G
        result = Vizier.get_catalogs("J/MNRAS/417/916/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            (name, source) = new_event('SNSDF'+row['SNSDF'], bibcode = "2011MNRAS.417..916G")
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
            add_quantity(name, 'redshift', row['zsp'] if row['zsp'] else row['zph'], source, kind = 'host')
            add_quantity(name, 'discoverdate', '20' + row['SNSDF'][:2] + '/' + row['SNSDF'][2:4], source, kind = 'host')
            add_quantity(name, 'hostoffsetang', row['Offset'], source, unit = 'arcseconds')
            add_quantity(name, 'claimedtype', row['Type'], source)
        journal_events()

        # 2013MNRAS.430.1746G
        result = Vizier.get_catalogs("J/MNRAS/430/1746/table4")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            (name, source) = new_event('SDSS'+row['SDSS'], bibcode = "2013MNRAS.430.1746G")
            add_quantity(name, 'ra', row['RAJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'dec', row['DEJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'discoverdate', row['Date'].replace('-', '/'), source)
            add_quantity(name, 'redshift', row['z'], source)
            add_quantity(name, 'claimedtype', row['Type'], source)
        journal_events()

        # 2014AJ....148...13R
        result = Vizier.get_catalogs("J/AJ/148/13/high_z")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            (name, source) = new_event(row['Name'], bibcode = "2014AJ....148...13R")
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
            add_quantity(name, 'discoverdate', '20' + row['Name'][3:5], source)
            add_quantity(name, 'redshift', row['zSN'], source, kind = 'heliocentric', error = row['e_zSN'])
            add_quantity(name, 'hostra', row['RAG'], source)
            add_quantity(name, 'hostdec', row['DEG'], source)
            add_quantity(name, 'hostoffsetang', row['ASep'], source, unit = 'arcseconds')
            add_quantity(name, 'redshift', row['zhost'], source, kind = 'host', error = row['e_zhost'])
        result = Vizier.get_catalogs("J/AJ/148/13/low_z")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            (name, source) = new_event(row['Name'], bibcode = "2014AJ....148...13R")
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
            add_quantity(name, 'discoverdate', '20' + row['Name'][3:5], source)
            add_quantity(name, 'redshift', row['zSN'], source, kind = 'heliocentric', error = row['e_zSN'])
            add_quantity(name, 'hostra', row['RAG'], source)
            add_quantity(name, 'hostdec', row['DEG'], source)
            add_quantity(name, 'hostoffsetang', row['ASep'], source, unit = 'arcseconds')
            add_quantity(name, 'redshift', row['zhost'], source, kind = 'host', error = row['e_zhost'])
        journal_events()

        # 2007ApJ...666..674M
        result = Vizier.get_catalogs("J/ApJ/666/674/table3")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            essname = 'ESSENCE '+row['ESSENCE']
            if row['SN']:
                name = 'SN'+row['SN']
            else:
                name = essname
            (name, source) = new_event(name, bibcode = "2007ApJ...666..674M")
            add_quantity(name, 'alias', essname, source)
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
            add_quantity(name, 'redshift', row['zSN'], source, error = row['e_zSN'], kind = 'heliocentric')
            add_quantity(name, 'redshift', row['zGal'], source, kind = 'host')
            add_quantity(name, 'claimedtype', row['SType'] if row['SType'] else row['Type'], source)
        journal_events()

        # 2013AcA....63....1K
        result = Vizier.get_catalogs("J/AcA/63/1/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            if 'OGLE' not in row['Name']:
                continue
            (name, source) = new_event(row['Name'], bibcode = "2013AcA....63....1K")
            add_quantity(name, 'alias', row['OGLEIV'], source)
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
            astrot = astrotime(float(row['Tmax']), format = 'jd').datetime
            add_quantity(name, 'maxdate', make_date_string(astrot.year, astrot.month, astrot.day), source)
        journal_events()
    
        # 2011MNRAS.410.1262W
        result = Vizier.get_catalogs("J/MNRAS/410/1262/tablea2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            (name, source) = new_event('SNLS-' + row['SN'], bibcode = "2011MNRAS.410.1262W")
            add_quantity(name, 'ra', row['_RA'], source, unit = 'floatdegrees')
            add_quantity(name, 'dec', row['_DE'], source, unit = 'floatdegrees')
            add_quantity(name, 'redshift', row['z'], source, error = row['e_z'], kind = 'heliocentric')
        journal_events()

        # 2012ApJ...755...61S
        result = Vizier.get_catalogs("J/ApJ/755/61/table3")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            sdssname = 'SDSS-II SN ' + row['SNID']
            if row['SN']:
                name = 'SN' + row['SN']
            else:
                name = sdssname
            (name, source) = new_event(name, bibcode = "2012ApJ...755...61S")
            add_quantity(name, 'alias', sdssname, source)
            add_quantity(name, 'hostra', row['RAJ2000'], source)
            add_quantity(name, 'hostdec', row['DEJ2000'], source)
            add_quantity(name, 'redshift', row['z'], source, error = row['e_z'] if is_number(row['e_z']) else '', kind = 'host')
        journal_events()

        # 2008AJ....135..348S
        result = Vizier.get_catalogs("J/AJ/135/348/SNe")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            sdssname = 'SDSS-II SN ' + row['SNID']
            if row['SN']:
                name = 'SN' + row['SN']
            else:
                name = sdssname
            (name, source) = new_event(name, bibcode = "2008AJ....135..348S")
            add_quantity(name, 'alias', sdssname, source)
            fra = Decimal(row['RAJ2000'])
            if fra < Decimal(0.0):
                fra = Decimal(360.0) + fra
            add_quantity(name, 'ra', str(fra), source, unit = 'floatdegrees')
            add_quantity(name, 'dec', row['DEJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'redshift', row['zsp'], source, kind = 'spectroscopic')
            add_quantity(name, 'claimedtype', row['Type'].replace('SN', '').strip(), source)
        journal_events()

        # 2010ApJ...713.1026D
        result = Vizier.get_catalogs("J/ApJ/713/1026/SNe")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            sdssname = 'SDSS-II SN ' + row['ID']
            if row['IAU']:
                name = 'SN' + row['IAU']
            else:
                name = sdssname
            (name, source) = new_event(name, bibcode = "2010ApJ...713.1026D")
            add_quantity(name, 'alias', sdssname, source)
            add_quantity(name, 'ra', row['RAJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'dec', row['DEJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'redshift', row['z'], source, kind = 'heliocentric')
        journal_events()

        # 2013ApJ...770..107C
        result = Vizier.get_catalogs("J/ApJ/770/107/galaxies")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            (name, source) = new_event(row['SN'], bibcode = "2013ApJ...770..107C")
            add_quantity(name, 'hostra', row['RAJ2000'], source)
            add_quantity(name, 'hostdec', row['DEJ2000'], source)
            add_quantity(name, 'redshift', row['z'], source, error = row['e_z'] if is_number(row['e_z']) else '', kind = 'host')
        journal_events()

        # 2011ApJ...738..162S
        result = Vizier.get_catalogs("J/ApJ/738/162/table3")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = 'SDSS-II SN ' + row['CID']
            (name, source) = new_event(name, bibcode = "2011ApJ...738..162S")
            fra = Decimal(row['RAJ2000'])
            if fra < Decimal(0.0):
                fra = Decimal(360.0) + fra
            add_quantity(name, 'ra', str(fra), source, unit = 'floatdegrees')
            add_quantity(name, 'dec', row['DEJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'redshift', row['z'], source, kind = 'spectroscopic', error = row['e_z'])
            add_quantity(name, 'claimedtype', 'Ia', source, probability = row['PzIa'])
        result = Vizier.get_catalogs("J/ApJ/738/162/table4")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = 'SDSS-II SN ' + row['CID']
            (name, source) = new_event(name, bibcode = "2011ApJ...738..162S")
            fra = Decimal(row['RAJ2000'])
            if fra < Decimal(0.0):
                fra = Decimal(360.0) + fra
            add_quantity(name, 'ra', str(fra), source, unit = 'floatdegrees')
            add_quantity(name, 'dec', row['DEJ2000'], source, unit = 'floatdegrees')
            add_quantity(name, 'redshift', row['zph'], source, kind = 'photometric')
            add_quantity(name, 'claimedtype', 'Ia', source, probability = row['PIa'])
        journal_events()

        # 2015MNRAS.446..943V
        snrtabs = ["ngc2403","ngc2903","ngc300","ngc3077","ngc4214","ngc4395","ngc4449","ngc5204",
            "ngc5585","ngc6946","ngc7793","m33","m74","m81","m82","m83","m101","m31"]
        for tab in tq(snrtabs, currenttask):
            result = Vizier.get_catalogs("J/MNRAS/446/943/" + tab)
            table = result[list(result.keys())[0]]
            table.convert_bytestring_to_unicode(python3_only=True)
            for ri, row in enumerate(tq(table, currenttask)):
                ra = row['RAJ2000'] if isinstance(row['RAJ2000'], str) else radec_clean(str(row['RAJ2000']), 'ra', unit = 'floatdegrees')[0]
                dec = row['DEJ2000'] if isinstance(row['DEJ2000'], str) else radec_clean(str(row['DEJ2000']), 'dec', unit = 'floatdegrees')[0]
                name = tab.upper() + 'SNR J' + rep_chars(ra, ' :.') + rep_chars(dec, ' :.')
                (name, source) = new_event(name, bibcode = "2015MNRAS.446..943V")
                add_quantity(name, 'ra', ra, source)
                add_quantity(name, 'dec', dec, source)
                add_quantity(name, 'host', tab.upper(), source)
        journal_events()

        # 2009ApJ...703..370C
        result = Vizier.get_catalogs("J/ApJ/703/370/tables")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            ra = row['RAJ2000']
            dec = row['DEJ2000']
            name = row['Gal'].replace(' ', '') + 'SNR J' + rep_chars(ra, ' .') + rep_chars(dec, ' .')
            (name, source) = new_event(name, bibcode = "2009ApJ...703..370C")
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
            add_quantity(name, 'host', row['Gal'], source)
        journal_events()

        # 2016ApJ...821...57D
        (name, source) = new_event('SN2013ge', bibcode = "2016ApJ...821...57D")
        result = Vizier.get_catalogs("J/ApJ/821/57/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            for band in ['UVW2', 'UVM2', 'UVW1', 'U', 'B', 'V']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = str(row["MJD"]), band = band, magnitude = row[bandtag],
                                   e_magnitude = row["e_" + bandtag], telescope = 'Swift', instrument = 'UVOT', source = source)
        result = Vizier.get_catalogs("J/ApJ/821/57/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            for band in ['B', 'V', 'R', 'I']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = str(row["MJD"]), band = band, magnitude = row[bandtag],
                                   e_magnitude = row["e_" + bandtag], instrument = 'CAO', source = source)
        result = Vizier.get_catalogs("J/ApJ/821/57/table3")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            for band in ['B', 'V', "r'", "i'"]:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    add_photometry(name, time = str(row["MJD"]), band = band, magnitude = row[bandtag],
                                   e_magnitude = row["e_" + bandtag], instrument = 'FLWO', source = source)
        result = Vizier.get_catalogs("J/ApJ/821/57/table4")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            for band in ['r', 'i', 'z']:
                bandtag = band + 'mag'
                if bandtag in row and is_number(row[bandtag]) and not isnan(float(row[bandtag])):
                    upp = False
                    if "l_" + bandtag in row and row["l_" + bandtag] == ">":
                        upp = True
                    add_photometry(name, time = str(row["MJD"]), band = band, magnitude = row[bandtag], upperlimit = upp,
                                   e_magnitude = row["e_" + bandtag] if is_number(row["e_" + bandtag]) else '',
                                   instrument = row["Inst"], source = source)
        journal_events()

        # 2004ApJ...607..665R
        result = Vizier.get_catalogs("J/ApJ/607/665/table1")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = row['Name'].replace('SN ', 'SN')
            (name, source) = new_event(name, bibcode = "2004ApJ...607..665R")
            add_quantity(name, 'alias', row['OName'], source)
            add_quantity(name, 'ra', row['RAJ2000'], source)
            add_quantity(name, 'dec', row['DEJ2000'], source)
        result = Vizier.get_catalogs("J/ApJ/607/665/table2")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = row['Name'].replace('SN ', 'SN')
            (name, source) = new_event(name, bibcode = "2004ApJ...607..665R")
            mjd = str(jd_to_mjd(Decimal(row['HJD'])))
            add_photometry(name, time = mjd, band = row['Filt'], magnitude = row['Vega'], system = 'Vega',
                           e_magnitude = row['e_Vega'], source = source)
        result = Vizier.get_catalogs("J/ApJ/607/665/table5")
        table = result[list(result.keys())[0]]
        table.convert_bytestring_to_unicode(python3_only=True)
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = row['Name'].replace('SN ', 'SN')
            (name, source) = new_event(name, bibcode = "2004ApJ...607..665R")
            add_quantity(name, 'redshift', row['z'], source, kind = 'spectroscopic')
        journal_events()

    if do_task(task, 'donations'):
        # Nicholl 04-01-16 donation
        with open("../sne-external/Nicholl-04-01-16/bibcodes.json", 'r') as f:
            bcs = json.loads(f.read())
    
        for datafile in sorted(glob("../sne-external/Nicholl-04-01-16/*.txt"), key=lambda s: s.lower()):
            inpname = os.path.basename(datafile).split('_')[0]
            name = add_event(inpname)
            bibcode = ''
            for bc in bcs:
                if inpname in bcs[bc]:
                    bibcode = bc
            if not bibcode:
                raise(ValueError('Bibcode not found!'))
            source = add_source(name, bibcode = bibcode)
            add_quantity(name, 'alias', inpname, source)
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
            for row in tq(tsvin, currenttask):
                name = 'MCSNR ' + row[0]
                (name, source) = new_event(name, bibcode = '2016A&A...585A.162M')
                ra = row[2]
                dec = row[3]
                add_quantity(name, 'alias', 'LMCSNR J' + rep_chars(ra, ' :.') + rep_chars(dec, ' :.'), source)
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
            for row in tq(tsvin, currenttask):
                name = 'MCSNR ' + row[0]
                (name, source) = new_event(name, refname = 'Pierre Maggi')
                ra = row[3]
                dec = row[4]
                add_quantity(name, 'alias', 'SMCSNR J' + ra.replace(':', '')[:6] + dec.replace(':', '')[:7], source)
                add_quantity(name, 'alias', row[1], source)
                add_quantity(name, 'alias', row[2], source)
                add_quantity(name, 'ra', row[3], source)
                add_quantity(name, 'dec', row[4], source)
                add_quantity(name, 'host', 'SMC', source)
        journal_events()
    
        # Galbany 04-18-16 donation
        folders = next(os.walk('../sne-external/galbany-04-18-16/'))[1]
        bibcode = '2016AJ....151...33G'
        for folder in folders:
            infofiles = glob("../sne-external/galbany-04-18-16/" + folder + "/*.info")
            photfiles = glob("../sne-external/galbany-04-18-16/" + folder + "/*.out*")
    
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
                            (name, source) = new_event(name, bibcode = bibcode)
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

        # Brown 05-14-16
        files = glob("../sne-external/brown-05-14-16/*.dat")
        for fi in tq(files, currenttask):
            name = os.path.basename(fi).split('_')[0]
            (name, source) = new_event(name, refname = 'Swift Supernovae', bibcode = '2014Ap&SS.354...89B',
                url = 'http://people.physics.tamu.edu/pbrown/SwiftSN/swift_sn.html')
            with open(fi, 'r') as f:
                lines = f.read().splitlines()
                for line in lines:
                    if not line or line[0] == '#':
                        continue
                    cols = list(filter(None, line.split()))
                    band = cols[0]
                    mjd = cols[1]
                    # Skip lower limit entries for now
                    if cols[2] == 'NULL' and cols[6] == 'NULL':
                        continue
                    isupp = cols[2] == 'NULL' and cols[6] != 'NULL'
                    mag = cols[2] if not isupp else cols[4]
                    e_mag = cols[3] if not isupp else ''
                    upp = '' if not isupp else True
                    add_photometry(name, time = mjd, magnitude = mag, e_magnitude = e_mag,
                                   upperlimit = upp, band = band, source = source,
                                   telescope = 'Swift', instrument = 'UVOT', system = 'Vega')
        journal_events()

        # Nicholl 05-03-16
        files = glob("../sne-external/nicholl-05-03-16/*.txt")
        (name, source) = new_event('SN2015bn', bibcode = '2016arXiv160304748N')
        add_quantity(name, 'alias', 'PS15ae', source)
        for fi in tq(files, currenttask):
            telescope = os.path.basename(fi).split('_')[1]
            with open(fi, 'r') as f:
                lines = f.read().splitlines()
                for li, line in enumerate(lines):
                    if not line or (line[0] == '#' and li != 0):
                        continue
                    cols = list(filter(None, line.split()))
                    if not cols:
                        continue
                    if li == 0:
                        bands = cols[1:]
                        continue

                    mjd = cols[0]
                    for ci, col in enumerate(cols[1::2]):
                        if not is_number(col):
                            continue

                        emag = cols[2*ci+2]
                        upp = ''
                        if not is_number(emag):
                            emag = ''
                            upp = True
                        add_photometry(name, time = mjd, magnitude = col, e_magnitude = emag,
                                       upperlimit = upp, band = bands[ci], source = source,
                                       telescope = telescope, instrument = 'UVOT' if telescope == 'Swift' else '',
                                       system = 'Vega' if telescope == 'Swift' else 'AB')
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
                (name, source) = new_event(name, bibcode = "2015A&A...579A..40S")
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
            for ri, row in enumerate(tq(tsvin, currenttask)):
                if ri == 0:
                    continue
                name = row[0].replace('SCP', 'SCP-')
                (name, source) = new_event(name, refname = 'Supernova Cosmology Project', url = 'http://supernova.lbl.gov/2009ClusterSurvey/')
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
            for ri, row in enumerate(tq(tsvin, currenttask)):
                name = 'SNLS-' + row[0]
                (name, source) = new_event(name, bibcode = '2006ApJ...645..841N')
                add_quantity(name, 'redshift', row[1], source, kind = 'spectroscopic')
                astrot = astrotime(float(row[4]) + 2450000., format = 'jd').datetime
                add_quantity(name, 'discoverdate', make_date_string(astrot.year, astrot.month, astrot.day), source)
        journal_events()
    
        # Anderson 2014
        for datafile in tq(sorted(glob("../sne-external/SNII_anderson2014/*.dat"), key=lambda s: s.lower()), currenttask):
            basename = os.path.basename(datafile)
            if not is_number(basename[:2]):
                continue
            if basename == '0210_V.dat':
                name = 'SN0210'
            else:
                name = ('SN20' if int(basename[:2]) < 50 else 'SN19') + basename.split('_')[0]
            (name, source) = new_event(name, bibcode = '2014ApJ...786...67A')
    
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
            for row in tq(tsvin, currenttask):
                name = row[0]
                (name, source) = new_event(name, bibcode = "2004A&A...415..863G")
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
            for r, row in enumerate(tq(data, currenttask)):
                if r == 0:
                    continue
                namesplit = row[0].split('/')
                name = namesplit[-1]
                if name.startswith('SN'):
                    name = name.replace(' ', '')
                (name, source) = new_event(name, bibcode = '2015MNRAS.449..451W')
                if len(namesplit) > 1:
                    add_quantity(name, 'alias', namesplit[0], source)
                add_quantity(name, 'claimedtype', row[1], source)
                add_photometry(name, time = row[2], band = row[4], magnitude = row[3], source = source)
        journal_events()

        # 2016MNRAS.459.1039T
        with open("../sne-external/2016MNRAS.459.1039T.tsv", 'r') as f:
            data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace = True)
            (name, source) = new_event('LSQ13zm', bibcode = '2016MNRAS.459.1039T')
            for r, row in enumerate(tq(data, currenttask)):
                if row[0][0] == '#':
                    bands = [x.replace('(err)', '') for x in row[3:-1]]
                    continue
                mjd = row[1]
                mags = [re.sub(r'\([^)]*\)', '', x) for x in row[3:-1]]
                upps = [True if '>' in x else '' for x in mags]
                mags = [x.replace('>', '') for x in mags]
                errs = [x[x.find("(")+1:x.find(")")] if "(" in x else '' for x in row[3:-1]]
                for mi, mag in enumerate(mags):
                    if not is_number(mag):
                        continue
                    add_photometry(name, time = mjd, band = bands[mi], magnitude = mag, e_magnitude = errs[mi],
                        instrument = row[-1], upperlimit = upps[mi], source = source)
        journal_events()

        # 2015ApJ...804...28G
        with open("../sne-external/2015ApJ...804...28G.tsv", 'r') as f:
            data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace = True)
            (name, source) = new_event('PS1-13arp', bibcode = '2015ApJ...804...28G')
            for r, row in enumerate(tq(data, currenttask)):
                if r == 0:
                    continue    
                mjd = row[1]
                mag = row[3]
                upp = True if '<' in mag else ''
                mag = mag.replace('<', '')
                err = row[4] if is_number(row[4]) else ''
                ins = row[5]
                add_photometry(name, time = mjd, band = row[0], magnitude = mag, e_magnitude = err,
                    instrument = ins, upperlimit = upp, source = source)
        journal_events()

        # 2016ApJ...819...35A
        with open("../sne-external/2016ApJ...819...35A.tsv", 'r') as f:
            data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace = True)
            for r, row in enumerate(tq(data, currenttask)):
                if row[0][0] == '#':
                    continue    
                (name, source) = new_event(row[0], bibcode = '2016ApJ...819...35A')
                add_quantity(name, 'ra', row[1], source)
                add_quantity(name, 'dec', row[2], source)
                add_quantity(name, 'redshift', row[3], source)
                add_quantity(name, 'discoverdate',
                    datetime.strptime(row[4], '%Y %b %d').isoformat().split('T')[0].replace('-', '/'), source)
        journal_events()

        # 2014ApJ...784..105W
        with open("../sne-external/2014ApJ...784..105W.tsv", 'r') as f:
            data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace = True)
            for r, row in enumerate(tq(data, currenttask)):
                if row[0][0] == '#':
                    continue    
                (name, source) = new_event(row[0], bibcode = '2014ApJ...784..105W')
                mjd = row[1]
                band = row[2]
                mag = row[3]
                err = row[4]
                add_photometry(name, time = mjd, band = row[2], magnitude = mag, e_magnitude = err,
                    instrument = 'WHIRC', telescope = 'WIYN 3.5 m', observatory = 'NOAO',
                    system = 'WHIRC', source = source)
        journal_events()

        # 2012MNRAS.425.1007B
        with open("../sne-external/2012MNRAS.425.1007B.tsv", 'r') as f:
            data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace = True)
            for r, row in enumerate(tq(data, currenttask)):
                if row[0][0] == '#':
                    bands = row[2:]
                    continue    
                (name, source) = new_event(row[0], bibcode = '2012MNRAS.425.1007B')
                mjd = row[1]
                mags = [x.split('±')[0].strip() for x in row[2:]] 
                errs = [x.split('±')[1].strip() if '±' in x else '' for x in row[2:]] 
                if row[0] == 'PTF09dlc':
                    ins = 'HAWK-I'
                    tel = 'VLT 8.1m'
                    obs = 'ESO'
                else:
                    ins = 'NIRI'
                    tel = 'Gemini North 8.2m'
                    obs = 'Gemini'
                    
                for mi, mag in enumerate(mags):
                    if not is_number(mag):
                        continue
                    add_photometry(name, time = mjd, band = bands[mi], magnitude = mag, e_magnitude = errs[mi],
                        instrument = ins, telescope = tel, observatory = obs,
                        system = 'Natural', source = source)
        journal_events()

        # 2014ApJ...783...28G
        with open("../sne-external/apj490105t2_ascii.txt", 'r') as f:
            data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace = True)
            for r, row in enumerate(tq(data, currenttask)):
                if row[0][0] == '#':
                    continue    
                (name, source) = new_event(row[0], bibcode = '2014ApJ...783...28G')
                add_quantity(name, 'alias', row[1], source)
                add_quantity(name, 'discoverdate', '20' + row[0][3:5], source)
                add_quantity(name, 'ra', row[2], source)
                add_quantity(name, 'dec', row[3], source)
                add_quantity(name, 'redshift', row[13] if is_number(row[13]) else row[10], source)
        journal_events()

        # 2005ApJ...634.1190H
        with open("../sne-external/2005ApJ...634.1190H.tsv", 'r') as f:
            data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace = True)
            for r, row in enumerate(tq(data, currenttask)):
                (name, source) = new_event('SNLS-' + row[0], bibcode = '2005ApJ...634.1190H')
                add_quantity(name, 'discoverdate', '20' + row[0][:2], source)
                add_quantity(name, 'ra', row[1], source)
                add_quantity(name, 'dec', row[2], source)
                add_quantity(name, 'redshift', row[5].replace('?', ''), source, error = row[6], kind = 'host')
                add_quantity(name, 'claimedtype', row[7].replace('SN', '').strip(':* '), source)
        journal_events()

        # 2014MNRAS.444.2133S
        with open("../sne-external/2014MNRAS.444.2133S.tsv", 'r') as f:
            data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace = True)
            for r, row in enumerate(tq(data, currenttask)):
                if row[0][0] == '#':
                    continue    
                name = row[0]
                if is_number(name[:4]):
                    name = 'SN' + name
                (name, source) = new_event(name, bibcode = '2014MNRAS.444.2133S')
                add_quantity(name, 'ra', row[1], source)
                add_quantity(name, 'dec', row[2], source)
                add_quantity(name, 'redshift', row[3], source, kind = 'host')
        journal_events()

        # 2009MNRAS.398.1041B
        with open("../sne-external/2009MNRAS.398.1041B.tsv", 'r') as f:
            data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace = True)
            for r, row in enumerate(tq(data, currenttask)):
                if row[0][0] == '#':
                    bands = row[2:-1]
                    continue    
                (name, source) = new_event('SN2008S', bibcode = '2009MNRAS.398.1041B')
                mjd = str(jd_to_mjd(Decimal(row[0])))
                mags = [x.split('±')[0].strip() for x in row[2:]] 
                upps = [('<' in x.split('±')[0]) for x in row[2:]] 
                errs = [x.split('±')[1].strip() if '±' in x else '' for x in row[2:]] 

                instrument = row[-1]
                    
                for mi, mag in enumerate(mags):
                    if not is_number(mag):
                        continue
                    add_photometry(name, time = mjd, band = bands[mi], magnitude = mag, e_magnitude = errs[mi],
                        instrument = ins, source = source)
        journal_events()

        # 2010arXiv1007.0011P
        with open("../sne-external/2010arXiv1007.0011P.tsv", 'r') as f:
            data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace = True)
            for r, row in enumerate(tq(data, currenttask)):
                if row[0][0] == '#':
                    bands = row[1:]
                    continue    
                (name, source) = new_event('SN2008S', bibcode = '2010arXiv1007.0011P')
                mjd = row[0]
                mags = [x.split('±')[0].strip() for x in row[1:]] 
                errs = [x.split('±')[1].strip() if '±' in x else '' for x in row[1:]] 

                for mi, mag in enumerate(mags):
                    if not is_number(mag):
                        continue
                    add_photometry(name, time = mjd, band = bands[mi], magnitude = mag, e_magnitude = errs[mi],
                        instrument = 'LBT', source = source)
        journal_events()

        # 2000ApJ...533..320G
        with open("../sne-external/2000ApJ...533..320G.tsv", 'r') as f:
            data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace = True)
            (name, source) = new_event('SN1997cy', bibcode = '2000ApJ...533..320G')
            for r, row in enumerate(tq(data, currenttask)):
                if row[0][0] == '#':
                    bands = row[1:-1]
                    continue
                mjd = str(jd_to_mjd(Decimal(row[0])))
                mags = row[1:len(bands)]
                for mi, mag in enumerate(mags):
                    if not is_number(mag):
                        continue
                    add_photometry(name, time = mjd, band = bands[mi], magnitude = mag,
                        observatory = 'Mount Stromlo', telescope = 'MSSSO', source = source, kcorrected = True)
        journal_events()

    # CCCP
    if do_task(task, 'cccp'):
        cccpbands = ['B', 'V', 'R', 'I']
        for datafile in sorted(glob("../sne-external/CCCP/apj407397*.txt"), key=lambda s: s.lower()):
            with open(datafile,'r') as f:
                tsvin = csv.reader(f, delimiter='\t', skipinitialspace=True)
                for r, row in enumerate(tsvin):
                    if r == 0:
                        continue
                    elif r == 1:
                        name = 'SN' + row[0].split('SN ')[-1]
                        (name, source) = new_event(name, bibcode = '2012ApJ...744...10K')
                    elif r >= 5:
                        mjd = str(Decimal(row[0]) + 53000)
                        for b, band in enumerate(cccpbands):
                            if row[2*b + 1]:
                                if not row[2*b + 2]:
                                    upplim = True
                                add_photometry(name, time = mjd, band = band, magnitude = row[2*b + 1].strip('>'),
                                    e_magnitude = row[2*b + 2], upperlimit = (not row[2*b + 2]), source = source)
    
        if archived_task('cccp'):
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
        for link in tq(links, currenttask):
            if 'sc_sn' in link['href']:
                (name, source) = new_event(link.text.replace(' ', ''), refname = 'CCCP', url = 'https://webhome.weizmann.ac.il/home/iair/sc_cccp.html')
    
                if archived_task('cccp'):
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
                        if archived_task('cccp'):
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
    
        for datafile in tq(sorted(glob("../sne-external/SUSPECT/*.html"), key=lambda s: s.lower()), currenttask):
            basename = os.path.basename(datafile)
            basesplit = basename.split('-')
            oldname = basesplit[1]
            name = add_event(oldname)
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
            secondarysource = add_source(name, refname = secondaryreference, url = secondaryrefurl, secondary = True)
            add_quantity(name, 'alias', oldname, secondarysource)
    
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
        for fname in tq(sorted(glob("../sne-external/cfa-input/*.dat"), key=lambda s: s.lower()), currenttask):
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
    
            oldname = snname(eventparts[0])
            name = add_event(oldname)
            secondaryname = 'CfA Supernova Archive'
            secondaryurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
            secondarysource = add_source(name, refname = secondaryname, url = secondaryurl, secondary = True, acknowledgment = cfaack)
            add_quantity(name, 'alias', oldname, secondarysource)
    
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
        for r, row in enumerate(tq(tsvin, currenttask)):
            if r <= 47:
                continue
    
            if row[0][:2] != 'sn':
                name = 'SN' + row[0].strip()
            else:
                name = row[0].strip()
    
            (name, source) = new_event(name, bibcode = '2012ApJS..200...12H')
            add_quantity(name, 'claimedtype', 'Ia', source)
            add_photometry(name, u_time = 'MJD', time = row[2].strip(), band = row[1].strip(),
                magnitude = row[6].strip(), e_magnitude = row[7].strip(), source = source)
        
        # Bianco 2014
        tsvin = open("../sne-external/bianco-2014-standard.dat", 'r')
        tsvin = csv.reader(tsvin, delimiter=' ', skipinitialspace=True)
        for row in tq(tsvin, currenttask):
            name = 'SN' + row[0]
    
            (name, source) = new_event(name, bibcode = '2014ApJS..213...19B')
            add_photometry(name, u_time = 'MJD', time = row[2], band = row[1], magnitude = row[3],
                e_magnitude = row[4], telescope = row[5], system = "Standard", source = source)
        f.close()
        journal_events()
    
    # New UCB import
    if do_task(task, 'ucb'):
        secondaryreference = "UCB Filippenko Group's Supernova Database (SNDB)"
        secondaryrefurl = "http://heracles.astro.berkeley.edu/sndb/info"
        secondaryrefbib = "2012MNRAS.425.1789S"
    
        jsontxt = load_cached_url("http://heracles.astro.berkeley.edu/sndb/download?id=allpubphot",
            '../sne-external/SNDB/allpub.json')
        if not jsontxt:
            continue
    
        photom = json.loads(jsontxt)
        photom = sorted(photom, key = lambda k: k['ObjName'])
        for phot in tq(photom, currenttask = currenttask):
            oldname = phot["ObjName"]
            name = add_event(oldname)
    
            secondarysource = add_source(name, refname = secondaryreference, url = secondaryrefurl, bibcode = secondaryrefbib, secondary = True)
            add_quantity(name, 'alias', oldname, secondarysource)
            sources = [secondarysource]
            if phot["Reference"]:
                sources += [add_source(name, bibcode = phot["Reference"])]
            sources = uniq_cdl(sources)
    
            if phot["Type"] and phot["Type"].strip() != "NoMatch":
                for ct in phot["Type"].strip().split(','):
                    add_quantity(name, 'claimedtype', ct.replace('-norm', '').strip(), sources)
            if phot["DiscDate"]:
                add_quantity(name, 'discoverdate', phot["DiscDate"].replace('-', '/'), sources)
            if phot["HostName"]:
                add_quantity(name, 'host', urllib.parse.unquote(phot["HostName"]).replace('*', ''), sources)
            filename = phot["Filename"] if phot["Filename"] else ''
    
            if not filename:
                raise(ValueError('Filename not found for SNDB phot!'))
            if not phot["PhotID"]:
                raise(ValueError('ID not found for SNDB phot!'))
    
            filepath = '../sne-external/SNDB/' + filename
            if archived_task('ucb') and os.path.isfile(filepath):
                with open(filepath, 'r') as f:
                    phottxt = f.read()
            else:
                session = requests.Session()
                response = session.get("http://heracles.astro.berkeley.edu/sndb/download?id=dp:" + str(phot["PhotID"]))
                phottxt = response.text
                with open(filepath, 'w') as f:
                    f.write(phottxt)
    
            tsvin = csv.reader(phottxt.splitlines(), delimiter=' ', skipinitialspace=True)
    
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
                    e_magnitude = e_magnitude, source = sources)

        journal_events()
        
    # Import SDSS
    if do_task(task, 'sdss'): 
        with open('../sne-external/SDSS/2010ApJ...708..661D.txt', 'r') as f:
            bibcodes2010 = f.read().split("\n")
        sdssbands = ['u', 'g', 'r', 'i', 'z']
        for fname in tq(sorted(glob("../sne-external/SDSS/*.sum"), key=lambda s: s.lower()), currenttask):
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
                        name = "SDSS-II SN " + row[3]
                    else:
                        name = "SN" + row[5]
                    (name, source) = new_event(name, bibcode = bibcode)
                    add_quantity(name, 'alias', "SDSS-II SN " + row[3], source)
    
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
        for ri, row in enumerate(tq(tsvin, currenttask)):
            if ri == 0 or not row:
                continue
            (name, source) = new_event(row[0], refname = reference, url = refurl)
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
            if not args.fullrefresh and archived_task('gaia') and os.path.isfile(fname):
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
    # VizieR catalogs exist for this: J/AJ/139/519, J/AJ/142/156. Should replace eventually.
    if do_task(task, 'csp'): 
        cspbands = ['u', 'B', 'V', 'g', 'r', 'i', 'Y', 'J', 'H', 'K']
        for fname in tq(sorted(glob("../sne-external/CSP/*.dat"), key=lambda s: s.lower()), currenttask):
            f = open(fname,'r')
            tsvin = csv.reader(f, delimiter='\t', skipinitialspace=True)
    
            eventname = os.path.basename(os.path.splitext(fname)[0])
    
            eventparts = eventname.split('opt+')
    
            name = snname(eventparts[0])
    
            reference = "Carnegie Supernova Project"
            refbib = "2010AJ....139..519C"
            refurl = "http://csp.obs.carnegiescience.edu/data"
            (name, source) = new_event(name, bibcode = refbib, refname = reference, url = refurl)
    
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
        itepbadsources = ['2004ApJ...602..571B']
    
        needsbib = []
        with open("../sne-external/itep-refs.txt",'r') as f:
            refrep = f.read().splitlines()
        refrepf = dict(list(zip(refrep[1::2], refrep[::2])))
        f = open("../sne-external/itep-lc-cat-28dec2015.txt",'r')
        tsvin = csv.reader(f, delimiter='|', skipinitialspace=True)
        curname = ''
        for r, row in enumerate(tq(tsvin, currenttask)):
            if r <= 1 or len(row) < 7:
                continue
            oldname = 'SN' + row[0].strip()
            mjd = str(jd_to_mjd(Decimal(row[1].strip())))
            band = row[2].strip()
            magnitude = row[3].strip()
            e_magnitude = row[4].strip()
            reference = row[6].strip().strip(',')
    
            if curname != oldname:
                curname = oldname
                name = add_event(oldname)
    
                secondaryreference = "Sternberg Astronomical Institute Supernova Light Curve Catalogue"
                secondaryrefurl = "http://dau.itep.ru/sn/node/72"
                secondarysource = add_source(name, refname = secondaryreference, url = secondaryrefurl, secondary = True)
                add_quantity(name, 'alias', oldname, secondarysource)
    
                year = re.findall(r'\d+', name)[0]
                add_quantity(name, 'discoverdate', year, secondarysource)
            if reference in refrepf:
                bibcode = unescape(refrepf[reference])
                source = add_source(name, bibcode = bibcode)
            else:
                needsbib.append(reference)
                source = add_source(name, refname = reference) if reference else ''
    
            if bibcode not in itepbadsources:
                add_photometry(name, time = mjd, band = band, magnitude = magnitude, e_magnitude = e_magnitude, source = secondarysource + ',' + source)
        f.close()
        
        # Write out references that could use a bibcode
        needsbib = list(OrderedDict.fromkeys(needsbib))
        with open('../itep-needsbib.txt', 'w') as f:
            f.writelines(["%s\n" % i for i in needsbib])
        journal_events()
    
    # Now import the Asiago catalog
    if do_task(task, 'asiago'): 
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
    
        for record in tq(records, currenttask):
            if len(record) > 1 and record[1] != '':
                oldname = snname("SN" + record[1]).strip('?')
    
                reference = 'Asiago Supernova Catalogue'
                refurl = 'http://graspa.oapd.inaf.it/cgi-bin/sncat.php'
                refbib = '1989A&AS...81..421B'
                (name, source) = new_event(oldname, refname = reference, url = refurl, bibcode = refbib, secondary = True)
    
                year = re.findall(r'\d+', oldname)[0]
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
                if (ra != ''):
                    add_quantity(name, 'ra', ra, source, unit = 'nospace')
                if (dec != ''):
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
        for row in tq(table, currenttask):
            row = convert_aq_output(row)
            name = 'SN' + row['SN']
    
            (name, source) = new_event(name, bibcode = bibcode)
    
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

    if do_task(task, 'fermi'):
        with open("../sne-external/1SC_catalog_v01.asc", 'r') as f:
            tsvin = csv.reader(f, delimiter=',')
            for ri, row in enumerate(tq(tsvin, currenttask)):
                if row[0].startswith('#'):
                    if len(row) > 1 and 'UPPER_LIMITS' in row[1]:
                        break
                    continue
                if 'Classified' not in row[1]:
                    continue
                name = row[0].replace('SNR', 'G')
                (name, source) = new_event(name, bibcode = '2016ApJS..224....8A')
                add_quantity(name, 'alias', row[0].replace('SNR', 'MWSNR'), source)
                add_quantity(name, 'ra', row[2], source, unit = 'floatdegrees')
                add_quantity(name, 'dec', row[3], source, unit = 'floatdegrees')
        journal_events()
    
    if do_task(task, 'tns'):
        session = requests.Session()
        csvtxt = load_cached_url("https://wis-tns.weizmann.ac.il/search?&num_page=1&format=html&sort=desc&order=id&format=csv&page=0",
            "../sne-external/TNS/index.csv")
        if not csvtxt:
            continue
        maxid = csvtxt.splitlines()[1].split(",")[0].strip('"')
        maxpages = ceil(int(maxid)/1000.)
    
        for page in tq(range(maxpages), currenttask):
            fname = '../sne-external/TNS/page-' + str(page).zfill(2) + '.csv'
            if archived_task('tns') and os.path.isfile(fname) and page < 7:
                with open(fname, 'r') as f:
                    csvtxt = f.read()
            else:
                with open(fname, 'w') as f:
                    session = requests.Session()
                    response = session.get("https://wis-tns.weizmann.ac.il/search?&num_page=1000&format=html&edit[type]=&edit[objname]=&edit[id]=&sort=asc&order=id&display[redshift]=1&display[hostname]=1&display[host_redshift]=1&display[source_group_name]=1&display[programs_name]=1&display[internal_name]=1&display[isTNS_AT]=1&display[public]=1&display[end_pop_period]=0&display[spectra_count]=1&display[discoverymag]=1&display[discmagfilter]=1&display[discoverydate]=1&display[discoverer]=1&display[sources]=1&display[bibcode]=1&format=csv&page=" + str(page))
                    csvtxt = response.text
                    f.write(csvtxt)
    
            tsvin = csv.reader(csvtxt.splitlines(), delimiter=',')
            for ri, row in enumerate(tq(tsvin, currenttask, leave = False)):
                if ri == 0:
                    continue
                if row[4] and 'SN' not in row[4]:
                    continue
                name = row[1].replace(' ', '')
                (name, source) = new_event(name, refname = 'Transient Name Server', url = 'https://wis-tns.weizmann.ac.il')
                if row[2] and row[2] != '00:00:00.00':
                    add_quantity(name, 'ra', row[2], source)
                if row[3] and row[3] != '+00:00:00.00':
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
                # Currently, all events listing all possible observers. TNS bug?
                #if row[9]:
                #    observers = row[9].split(',')
                #    for observer in observers:
                #        add_quantity(name, 'observer', observer.strip(), source)
                if row[10]:
                    add_quantity(name, 'alias', row[10], source)
                if row[8] and row[14] and row[15] and row[16]:
                    survey = row[8]
                    magnitude = row[14]
                    band = row[15].split('-')[0]
                    mjd = astrotime(row[16]).mjd
                    add_photometry(name, time = mjd, magnitude = magnitude, band = band, survey = survey, source = source)
                if row[16]:
                    date = row[16].split()[0].replace('-', '/')
                    if date != '0000/00/00':
                        date = date.replace('/00', '')
                        time = row[16].split()[1]
                        if time != '00:00:00':
                            ts = time.split(':')
                            date += pretty_num(timedelta(hours = int(ts[0]), minutes = int(ts[1]), seconds = int(ts[2])).total_seconds()/(24*60*60), sig=6).lstrip('0')
                        add_quantity(name, 'discoverdate', date, source)
                if args.update:
                    journal_events()
        journal_events()
    
    if do_task(task, 'rochester'): 
        rochestermirrors = ['http://www.rochesterastronomy.org/', 'http://www.supernova.thistlethwaites.com/']
        rochesterpaths = ['snimages/snredshiftall.html', 'sn2016/snredshift.html', 'snimages/snredboneyard.html']
        rochesterupdate = [False, True, True]
    
        for p, path in enumerate(tq(rochesterpaths, currenttask)):
            if args.update and not rochesterupdate[p]:
                continue
    
            for mirror in rochestermirrors:
                filepath = '../sne-external/rochester/' + os.path.basename(path)
                html = load_cached_url(mirror + path, filepath, failhard = (mirror != rochestermirrors[-1]))
                if html:
                    break
            if not html:
                continue
    
            soup = BeautifulSoup(html, "html5lib")
            rows = soup.findAll('tr')
            secondaryreference = "Latest Supernovae"
            secondaryrefurl = "http://www.rochesterastronomy.org/snimages/snredshiftall.html"
            for r, row in enumerate(tq(rows, currenttask)):
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
                        oldname = aka
                        name = add_event(aka)
                    elif len(aka) == 4 and is_number(aka[:4]):
                        aka = 'SN' + aka
                        oldname = aka
                        name = add_event(aka)
    
                ra = str(cols[3].contents[0]).strip()
                dec = str(cols[4].contents[0]).strip()

                sn = re.sub('<[^<]+?>', '', str(cols[0].contents[0])).strip()
                if is_number(sn.strip('?')):
                    sn = 'SN' + sn.strip('?') + 'A'
                elif len(sn) == 4 and is_number(sn[:4]):
                    sn = 'SN' + sn
                if not name:
                    if not sn:
                        continue
                    if sn[:8] == 'MASTER J':
                        sn = sn.replace('MASTER J', 'MASTER OT J').replace('SNHunt', 'SNhunt')
                    if 'POSSIBLE' in sn.upper() and ra and dec:
                        sn = 'PSN J' + ra.replace(':', '').replace('.', '') + dec.replace(':', '').replace('.', '')
                    oldname = sn
                    name = add_event(sn)
    
                reference = cols[12].findAll('a')[0].contents[0].strip()
                refurl = cols[12].findAll('a')[0]['href'].strip()
                source = add_source(name, refname = reference, url = refurl)
                secondarysource = add_source(name, refname = secondaryreference, url = secondaryrefurl, secondary = True)
                sources = uniq_cdl(list(filter(None, [source, secondarysource])))
                add_quantity(name, 'alias', oldname, sources)
                add_quantity(name, 'alias', sn, sources)
    
                if cols[14].contents:
                    if aka == 'SNR G1.9+0.3':
                        aka = 'G001.9+00.3'
                    if aka[:4] == 'PS1 ':
                        aka = 'PS1-' + aka[4:]
                    if aka[:8] == 'MASTER J':
                        aka = aka.replace('MASTER J', 'MASTER OT J').replace('SNHunt', 'SNhunt')
                    if 'POSSIBLE' in aka.upper() and ra and dec:
                        aka = 'PSN J' + ra.replace(':', '').replace('.', '') + dec.replace(':', '').replace('.', '')
                    add_quantity(name, 'alias', aka, sources)
    
                if str(cols[1].contents[0]).strip() != 'unk':
                    add_quantity(name, 'claimedtype', str(cols[1].contents[0]).strip(' :,'), sources)
                if str(cols[2].contents[0]).strip() != 'anonymous':
                    add_quantity(name, 'host', str(cols[2].contents[0]).strip(), sources)
                add_quantity(name, 'ra', ra, sources)
                add_quantity(name, 'dec', dec, sources)
                if str(cols[6].contents[0]).strip() not in ['2440587', '2440587.292']:
                    astrot = astrotime(float(str(cols[6].contents[0]).strip()), format='jd').datetime
                    add_quantity(name, 'discoverdate', make_date_string(astrot.year, astrot.month, astrot.day), sources)
                if str(cols[7].contents[0]).strip() not in ['2440587', '2440587.292']:
                    astrot = astrotime(float(str(cols[7].contents[0]).strip()), format='jd')
                    if (float(str(cols[8].contents[0]).strip()) <= 90.0 and
                        not any('GRB' in x for x in get_aliases(name))):
                        add_photometry(name, time = str(astrot.mjd), magnitude = str(cols[8].contents[0]).strip(), source = sources)
                if cols[11].contents[0] != 'n/a':
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
                    (name, secondarysource) = new_event(name, refname = secondaryreference, url = secondaryrefurl, secondary = True)
    
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
                            source = add_source(name, refname = reference)
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
        for b, bn in enumerate(tq(basenames, currenttask)):
            if args.update and not ogleupdate[b]:
                continue
    
            filepath = '../sne-external/OGLE-' + bn.replace('/', '-') + '-transients.html'
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
                        datafnames.append(bn.replace('/', '-') + '-' + a['href'].replace('/', '-'))
    
            ec = -1
            reference = 'OGLE-IV Transient Detection System'
            refurl = 'http://ogle.astrouw.edu.pl/ogle4/transients/transients.html'
            for br in tq(breaks, currenttask):
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
                    if not args.fullrefresh and archived_task('ogle') and os.path.isfile(fname):
                        with open(fname, 'r') as f:
                            csvtxt = f.read()
                    else:
                        response = urllib.request.urlopen(datalinks[ec])
                        with open(fname, 'w') as f:
                            csvtxt = response.read().decode('utf-8')
                            f.write(csvtxt)
    
                    lcdat = csvtxt.splitlines()
                    sources = [add_source(name, refname = reference, url = refurl)]
                    add_quantity(name, 'alias', name, sources[0])
                    if atelref and atelref != 'ATel#----':
                        sources.append(add_source(name, refname = atelref, url = atelurl))
                    sources = uniq_cdl(sources)
    
                    if name.startswith('OGLE'):
                        if name[4] == '-':
                            if is_number(name[5:9]):
                                add_quantity(name, 'discoverdate', name[5:9], sources)
                        else:
                            if is_number(name[4:6]):
                                add_quantity(name, 'discoverdate', '20' + name[4:6], sources)
    
                    # RA and Dec from OGLE pages currently not reliable
                    #add_quantity(name, 'ra', ra, sources)
                    #add_quantity(name, 'dec', dec, sources)
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
                        add_photometry(name, time = mjd, band = 'I', magnitude = magnitude, e_magnitude = e_magnitude,
                            system = 'Vega', source = sources, upperlimit = upperlimit)
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
                (name, source) = new_event(name, bibcode = '2010A&A...523A...7G')
                band = row[1]
                mjd = row[2]
                sig = get_sig_digits(flux.split('E')[0])+1
                # Conversion comes from SNLS-Readme
                # NOTE: Datafiles available for download suggest different zeropoints than 30, need to inquire.
                magnitude = pretty_num(30.0-2.5*log10(float(flux)), sig = sig)
                e_magnitude = pretty_num(2.5*log10(1.0 + float(err)/float(flux)), sig = sig)
                #e_magnitude = pretty_num(2.5*(log10(float(flux) + float(err)) - log10(float(flux))), sig = sig)
                add_photometry(name, time = mjd, band = band, magnitude = magnitude, e_magnitude = e_magnitude, counts = flux,
                    e_counts = err, source = source)
        journal_events()
    
    if do_task(task, 'psthreepi'):
        fname = '../sne-external/3pi/page00.html'
        html = load_cached_url("http://psweb.mp.qub.ac.uk/ps1threepi/psdb/public/?page=1&sort=followup_flag_date", fname, write = False)
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
            if args.update:
                continue
            warnings.warn("Pan-STARRS 3pi offline, using local files only.")
            with open(fname, 'r') as f:
                html = f.read()
            bs = BeautifulSoup(html, "html5lib")
            div = bs.find('div', {"class":"pagination"})
            links = div.findAll('a')
        else:
            with open(fname, 'w') as f:
                f.write(html)
    
        numpages = int(links[-2].contents[0])
        oldnumpages = len(glob('../sne-external/3pi/page*'))
        for page in tq(range(1,numpages), currenttask):
            fname = '../sne-external/3pi/page' + str(page).zfill(2) + '.html'
            if offline:
                if not os.path.isfile(fname):
                    continue
                with open(fname, 'r') as f:
                    html = f.read()
            else:
                if not args.fullrefresh and archived_task('psthreepi') and page < oldnumpages and os.path.isfile(fname) :
                    with open(fname, 'r') as f:
                        html = f.read()
                else:
                    response = urllib.request.urlopen("http://psweb.mp.qub.ac.uk/ps1threepi/psdb/public/?page=" + str(page) + "&sort=followup_flag_date")
                    with open(fname, 'w') as f:
                        html = response.read().decode('utf-8')
                        f.write(html)
    
            bs = BeautifulSoup(html, "html5lib")
            trs = bs.findAll('tr')
            for tr in tq(trs, currenttask):
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
                    continue
    
                name = ''
                for alias in aliases:
                    if alias[:2] == 'SN':
                        name = alias
                if not name:
                    name = psname
                name = add_event(name)
                sources = [add_source(name, refname = 'Pan-STARRS 3Pi', url = 'http://psweb.mp.qub.ac.uk/ps1threepi/psdb/')]
                add_quantity(name, 'alias', name, sources[0])
                for ref in refs:
                    sources.append(add_source(name, refname = ref[0], url = ref[1]))
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
                if offline:
                    if not os.path.isfile(fname2):
                        continue
                    with open(fname2, 'r') as f:
                        html2 = f.read()
                else:
                    if archived_task('psthreepi') and os.path.isfile(fname2):
                        with open(fname2, 'r') as f:
                            html2 = f.read()
                    else:
                        pslink = 'http://psweb.mp.qub.ac.uk/ps1threepi/psdb/public/' + pslink
                        try:
                            session2 = requests.Session()
                            response2 = session2.get(pslink)
                        except:
                            offline = True
                            if not os.path.isfile(fname2):
                                continue
                            with open(fname2, 'r') as f:
                                html2 = f.read()
                        else:
                            html2 = response2.text
                            with open(fname2, 'w') as f:
                                f.write(html2)
    
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
                        add_photometry(name, time = str(obs[0]), band = nslabels[li], magnitude = str(obs[1]), e_magnitude = str(obs[2]), source = source,
                            telescope = 'Pan-STARRS1')
                for li, line in enumerate(nslines[2*len(nslabels):]):
                    if not line:
                        continue
                    for obs in line:
                        add_photometry(name, time = str(obs[0]), band = nslabels[li], magnitude = str(obs[1]), upperlimit = True, source = source,
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

            # Only run first page for Travis
            if args.travis:
                break
    
    if do_task(task, 'psmds'):
        with open('../sne-external/MDS/apj506838t1_mrt.txt') as f:
            for ri, row in enumerate(tq(f.read().splitlines(), currenttask)):
                if ri < 35:
                    continue
                cols = [x.strip() for x in row.split(',')]
                (name, source) = new_event(cols[0], bibcode = '2015ApJ...799..208S')
                add_quantity(name, 'ra', cols[2], source)
                add_quantity(name, 'dec', cols[3], source)
                astrot = astrotime(float(cols[4]), format='mjd').datetime
                add_quantity(name, 'discoverdate', make_date_string(astrot.year, astrot.month, astrot.day), source)
                add_quantity(name, 'redshift', cols[5], source, kind = 'spectroscopic')
                add_quantity(name, 'claimedtype', 'II P', source)
        journal_events()

    if do_task(task, 'psst'):
        # 2016arXiv160204156S
        with open("../sne-external/2016arXiv160204156S-tab1.tsv", 'r') as f:
            data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace = True)
            for r, row in enumerate(tq(data, currenttask)):
                if row[0][0] == '#':
                    continue
                (name, source) = new_event(row[0], bibcode = '2016arXiv160204156S')
                add_quantity(name, 'claimedtype', row[3].replace('SN', '').strip('() '), source)
                add_quantity(name, 'redshift', row[5].strip('() '), source, kind = 'spectroscopic')
        with open("../sne-external/2016arXiv160204156S-tab2.tsv", 'r') as f:
            data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace = True)
            for r, row in enumerate(tq(data, currenttask)):
                if row[0][0] == '#':
                    continue
                (name, source) = new_event(row[0], bibcode = '2016arXiv160204156S')
                add_quantity(name, 'ra', row[1], source)
                add_quantity(name, 'dec', row[2], source)
                mldt = astrotime(float(row[4]), format = 'mjd').datetime
                discoverdate = make_date_string(mldt.year, mldt.month, mldt.day)
                add_quantity(name, 'discoverdate', discoverdate, source)
        journal_events()

        # 1606.04795
        with open("../sne-external/1606.04795.tsv", 'r') as f:
            data = csv.reader(f, delimiter='\t', quotechar='"', skipinitialspace = True)
            for r, row in enumerate(tq(data, currenttask)):
                if row[0][0] == '#':
                    continue
                (name, source) = new_event(row[0], refname = 'Smartt et al. 2016', url = 'http://arxiv.org/abs/1606.04795')
                add_quantity(name, 'ra', row[1], source)
                add_quantity(name, 'dec', row[2], source)
                mldt = astrotime(float(row[3]), format = 'mjd').datetime
                discoverdate = make_date_string(mldt.year, mldt.month, mldt.day)
                add_quantity(name, 'discoverdate', discoverdate, source)
                add_quantity(name, 'claimedtype', row[6], source)
                add_quantity(name, 'redshift', row[7], source, kind = 'spectroscopic')
                for alias in [x.strip() for x in row[8].split(',')]:
                    add_quantity(name, 'alias', alias, source)
        journal_events()

    if do_task(task, 'grb'):
        csvtxt = load_cached_url('http://grb.pa.msu.edu/grbcatalog/download_data?cut_0_min=10&cut_0=BAT%20T90&cut_0_max=100000&num_cuts=1&no_date_cut=True',
            '../sne-external/GRB-catalog/catalog.csv')
        if not csvtxt:
            continue
        data = csv.reader(csvtxt.splitlines(), delimiter=',', quotechar='"', skipinitialspace = True)
        for r, row in enumerate(tq(data, currenttask)):
            if r == 0:
                continue
            (name, source) = new_event('GRB ' + row[0], refname = 'Gamma-ray Bursts Catalog', url = 'http://grbcatalog.org')
            add_quantity(name, 'ra', row[2], source, unit = 'floatdegrees')
            add_quantity(name, 'dec', row[3], source, unit = 'floatdegrees')
            add_quantity(name, 'redshift', row[8], source)
        journal_events()
    
    if do_task(task, 'crts'):
        crtsnameerrors = ['2011ax']
    
        folders = ["catalina", "MLS", "SSS"]
        for fold in tq(folders, currenttask):
            html = load_cached_url("http://nesssi.cacr.caltech.edu/" + fold + "/AllSN.html", '../sne-external/CRTS/' + fold + '.html')
            if not html:
                continue
            bs = BeautifulSoup(html, "html5lib")
            trs = bs.findAll('tr')
            for tri, tr in enumerate(tq(trs, currenttask)):
                tds = tr.findAll('td')
                if not tds:
                    continue
                refs = []
                aliases = []
                crtsname = ''
                ra = ''
                dec = ''
                lclink = ''
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
                            ind = ai+1
                            if aliases[ai+1] in ['SDSS']:
                                ind = ai+2
                            elif aliases[ai+1] in ['gal', 'obj', 'object', 'source']:
                                ind = ai-1
                            if '>' in aliases[ind]:
                                hostupper = True
                            hostmag = aliases[ind].strip('>~').replace(',', '.').replace('m', '.')
                        continue
                    if is_number(alias[:4]) and alias[:2] == '20' and len(alias) > 4:
                        name = 'SN' + alias
                    if (('asassn' in alias and len(alias) > 6) or ('ptf' in alias and len(alias) > 3) or
                        ('ps1' in alias and len(alias) > 3) or 'snhunt' in alias or
                        ('mls' in alias and len(alias) > 3) or 'gaia' in alias or ('lsq' in alias and len(alias) > 3)):
                        alias = alias.replace('SNHunt', 'SNhunt')
                        validaliases.append(alias)
                if not name:
                    name = crtsname
                (name, source) = new_event(name, refname = 'Catalina Sky Survey', bibcode = '2009ApJ...696..870D',
                    url = 'http://nesssi.cacr.caltech.edu/catalina/AllSN.html')
                for alias in validaliases:
                    add_quantity(name, 'alias', alias, source)
                add_quantity(name, 'ra', ra, source, unit = 'floatdegrees')
                add_quantity(name, 'dec', dec, source, unit = 'floatdegrees')
    
                if hostmag:
                    # 1.0 magnitude error based on Drake 2009 assertion that SN are only considered real if they are 2 mags brighter than host.
                    add_photometry(name, band = 'C', magnitude = hostmag, e_magnitude = 1.0, source = source, host = True,
                        telescope = 'Catalina Schmidt', upperlimit = hostupper)
    
                fname2 = '../sne-external/' + fold + '/' + lclink.split('.')[-2].rstrip('p').split('/')[-1] + '.html'
                if not args.fullrefresh and archived_task('crts') and os.path.isfile(fname2):
                    with open(fname2, 'r') as f:
                        html2 = f.read()
                else:
                    try:
                        with open(fname2, 'w') as f:
                            response2 = urllib.request.urlopen(lclink)
                            html2 = response2.read().decode('utf-8')
                            f.write(html2)
                    except:
                        continue
    
                lines = html2.splitlines()
                for line in lines:
                    if 'javascript:showx' in line:
                        search = re.search("showx\('(.*?)'\)", line)
                        if not search:
                            continue
                        mjdstr = search.group(1).split('(')[0].strip()
                        if not is_number(mjdstr):
                            continue
                        mjd = str(Decimal(mjdstr) + Decimal(53249.0))
                    else:
                        continue
                    if 'javascript:showy' in line:
                        mag = re.search("showy\('(.*?)'\)", line).group(1)
                    if 'javascript:showz' in line:
                        err = re.search("showz\('(.*?)'\)", line).group(1)
                    if not is_number(mag) or (err and not is_number(err)):
                        continue
                    add_photometry(name, time = mjd, band = 'C', magnitude = mag, source = source, includeshost = True,
                        telescope = 'Catalina Schmidt', e_magnitude = err if float(err) > 0.0 else '', upperlimit = (float(err) == 0.0))
                if args.update:
                    journal_events()
            if args.travis and tri > travislimit:
                break
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
        for tr in tq(trs, currenttask):
            cols = [str(x.text) for x in tr.findAll('td')]
            if not cols:
                continue
            name = re.sub('<[^<]+?>', '', cols[4]).strip().replace(' ', '').replace('SNHunt', 'SNhunt')
            (name, source) = new_event(name, refname = 'Supernova Hunt', url = 'http://nesssi.cacr.caltech.edu/catalina/current.html')
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
        with open("../sne-external/NED26.05.1-D-12.1.0-20160501.csv", 'r') as f:
            data = sorted(list(csv.reader(f, delimiter=',', quotechar='"'))[13:], key=lambda x: (x[9], x[3]))
            reference = "NED-D"
            refurl = "http://ned.ipac.caltech.edu/Library/Distances/"
            nedddict = OrderedDict()
            olddistname = ''
            for r, row in enumerate(tq(data, currenttask = currenttask)):
                if r <= 12:
                    continue
                distname = row[3]
                name = name_clean(distname) 
                distmod = row[4]
                moderr = row[5]
                dist = row[6]
                bibcode = unescape(row[8])
                snname = name_clean(row[9])
                redshift = row[10]
                cleanhost = ''

                if name != snname and (name + ' HOST' != snname):
                    cleanhost = host_clean(distname)
                    if cleanhost.endswith(' HOST'):
                        cleanhost = ''
                    if not is_number(dist):
                        print(dist)
                    if dist:
                        nedddict.setdefault(cleanhost,[]).append(Decimal(dist))
    
                if snname and 'HOST' not in snname:
                    (snname, secondarysource) = new_event(snname, refname = reference, url = refurl, secondary = True)
                    if bibcode:
                        source = add_source(snname, bibcode = bibcode)
                        sources = uniq_cdl([source, secondarysource])
                    else:
                        sources = secondarysource

                    if name == snname:
                        if redshift:
                            add_quantity(snname, 'redshift', redshift, sources)
                        if dist:
                            add_quantity(snname, 'comovingdist', dist, sources)
                            if not redshift:
                                try:
                                    redshift = pretty_num(z_at_value(cosmo.comoving_distance, float(dist) * un.Mpc, zmax = 5.0), sig = get_sig_digits(str(dist)))
                                except (KeyboardInterrupt, SystemExit):
                                    raise
                                except:
                                    pass
                                else:
                                    cosmosource = add_source(name, bibcode = '2015arXiv150201589P')
                                    add_quantity(snname, 'redshift', redshift, uniq_cdl(sources.split(',') + [cosmosource])) 

                    if cleanhost:
                        add_quantity(snname, 'host', cleanhost, sources)

                    if args.update and olddistname != distname:
                        journal_events()
                olddistname = distname
            journal_events()
    
    # Import CPCS
    if do_task(task, 'cpcs'):
        jsontxt = load_cached_url("http://gsaweb.ast.cam.ac.uk/followup/list_of_alerts?format=json&num=100000&published=1&observed_only=1&hashtag=JG_530ad9462a0b8785bfb385614bf178c6",
            "../sne-external/CPCS/index.json")
        if not jsontxt:
            continue
        alertindex = json.loads(jsontxt, object_pairs_hook=OrderedDict)
        ids = [x["id"] for x in alertindex]
    
        for i, ai in enumerate(tq(ids, currenttask)):
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
                oldname = name
                name = add_event(name)
            else:
                continue
    
            secondarysource = add_source(name, refname = 'Cambridge Photometric Calibration Server', url = 'http://gsaweb.ast.cam.ac.uk/followup/', secondary = True)
            add_quantity(name, 'alias', oldname, secondarysource)
            add_quantity(name, 'ra', str(alertindex[i]['ra']), secondarysource, unit = 'floatdegrees')
            add_quantity(name, 'dec', str(alertindex[i]['dec']), secondarysource, unit = 'floatdegrees')
    
            alerturl = "http://gsaweb.ast.cam.ac.uk/followup/get_alert_lc_data?alert_id=" + str(ai)
            source = add_source(name, refname = 'CPCS Alert ' + str(ai), url = alerturl)
            fname = '../sne-external/CPCS/alert-' + str(ai).zfill(2) + '.json'
            if archived_task('cpcs') and os.path.isfile(fname):
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
    
        if archived_task('ptf'):
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
                    (name, source) = new_event(name, bibcode = '2012PASP..124..668Y')
                    add_quantity(name, 'alias', alias, source)
                else:
                    (name, source) = new_event(name, bibcode = '2012PASP..124..668Y')
        
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
                (name, source) = new_event(name, bibcode = '2016arXiv160408207P')
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

    if do_task(task, 'des'):
        html = load_cached_url("https://portal.nersc.gov/des-sn/transients/", "../sne-external/DES/transients.html")
        if not html:
            continue
        bs = BeautifulSoup(html, "html5lib")
        trs = bs.find('tbody').findAll('tr')
        for tri, tr in enumerate(tq(trs, currenttask)):
            name = ''
            source = ''
            if tri == 0:
                continue
            tds = tr.findAll('td')
            for tdi, td in enumerate(tds):
                if tdi == 0:
                    name = add_event(td.text.strip())
                if tdi == 1:
                    (ra, dec) = [x.strip() for x in td.text.split('\xa0')]
                if tdi == 6:
                    atellink = td.find('a')
                    if atellink:
                        atellink = atellink['href']
                    else:
                        atellink = ''

            sources = [add_source(name, url = 'https://portal.nersc.gov/des-sn/', refname = 'DES Bright Transients',
                acknowledgment = 'http://www.noao.edu/noao/library/NOAO_Publications_Acknowledgments.html#DESdatause')]
            if atellink:
                sources.append(add_source(name, refname = 'ATel ' + atellink.split('=')[-1], url = atellink))
            sources += [add_source(name, bibcode = '2012ApJ...753..152B'),
                add_source(name, bibcode = '2015AJ....150..150F'),
                add_source(name, bibcode = '2015AJ....150...82G'),
                add_source(name, bibcode = '2015AJ....150..172K')]
            sources = ','.join(sources)
            add_quantity(name, 'alias', name, sources)
            add_quantity(name, 'ra', ra, sources)
            add_quantity(name, 'dec', dec, sources)

            html2 = load_cached_url("https://portal.nersc.gov/des-sn/transients/" + name, "../sne-external/DES/" + name + ".html")
            if not html2:
                continue
            lines = html2.splitlines()
            for line in lines:
                if 'var data = ' in line:
                    jsontxt = json.loads(line.split('=')[-1].rstrip(';'))
                    for i, band in enumerate(jsontxt['band']):
                        add_photometry(name, time = jsontxt['mjd'][i], magnitude = jsontxt['mag'][i], e_magnitude = jsontxt['mag_error'][i],
                            band = band, observatory = 'CTIO', telescope = 'Blanco 4m', instrument = 'DECam',
                            upperlimit = True if float(jsontxt['snr'][i]) <= 3.0 else '', source = sources)
        journal_events()

    if do_task(task, 'asassn'):
        html = load_cached_url("http://www.astronomy.ohio-state.edu/~assassin/sn_list.html", "../sne-external/ASASSN/sn_list.html")
        if not html:
            continue
        bs = BeautifulSoup(html, "html5lib")
        trs = bs.find('table').findAll('tr')
        for tri, tr in enumerate(tq(trs, currenttask)):
            name = ''
            source = ''
            ra = ''
            dec = ''
            redshift = ''
            hostoff = ''
            claimedtype = ''
            host = ''
            atellink = ''
            typelink = ''
            if tri == 0:
                continue
            tds = tr.findAll('td')
            for tdi, td in enumerate(tds):
                if tdi == 1:
                    name = add_event(td.text.strip())
                    atellink = td.find('a')
                    if atellink:
                        atellink = atellink['href']
                    else:
                        atellink = ''
                if tdi == 2:
                    discdate = td.text.replace('-', '/')
                if tdi == 3:
                    ra = td.text
                if tdi == 4:
                    dec = td.text
                if tdi == 5:
                    redshift = td.text
                if tdi == 8:
                    hostoff = td.text
                if tdi == 9:
                    claimedtype = td.text
                    typelink = td.find('a')
                    if typelink:
                        typelink = typelink['href']
                    else:
                        typelink = ''
                if tdi == 12:
                    host = td.text

            sources = [add_source(name, url = 'http://www.astronomy.ohio-state.edu/~assassin/sn_list.html', refname = 'ASAS-SN Supernovae')]
            typesources = sources[:]
            if atellink:
                sources.append(add_source(name, refname = 'ATel ' + atellink.split('=')[-1], url = atellink))
            if typelink:
                typesources.append(add_source(name, refname = 'ATel ' + typelink.split('=')[-1], url = typelink))
            sources = ','.join(sources)
            typesources = ','.join(typesources)
            add_quantity(name, 'alias', name, sources)
            add_quantity(name, 'discoverdate', discdate, sources)
            add_quantity(name, 'ra', ra, sources, unit = 'floatdegrees')
            add_quantity(name, 'dec', dec, sources, unit = 'floatdegrees')
            add_quantity(name, 'redshift', redshift, sources)
            add_quantity(name, 'hostoffsetang', hostoff, sources, unit = 'arcseconds')
            for ct in claimedtype.split('/'):
                if ct != 'Unk':
                    add_quantity(name, 'claimedtype', ct, typesources)
            if host != 'Uncatalogued':
                add_quantity(name, 'host', host, sources)
        journal_events()

    if do_task(task, 'snf'):
        with open('../sne-external/SNF/snf-aliases.csv') as f:
            for row in [x.split(',') for x in f.read().splitlines()]:
                (name, source) = new_event(row[0], bibcode = oscbibcode, refname = oscname, url = oscurl, secondary = True)
                add_quantity(name, 'alias', row[1], source)
        journal_events()
    
    if do_task(task, 'asiagospectra'):
        html = load_cached_url("http://sngroup.oapd.inaf.it./cgi-bin/output_class.cgi?sn=1990", "../sne-external-spectra/Asiago/spectra.html")
        if not html:
            continue
        bs = BeautifulSoup(html, "html5lib")
        trs = bs.findAll('tr')
        for tr in tq(trs, currenttask):
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
                    oldname = name
                    name = add_event(name)
                    reference = 'Asiago Supernova Catalogue'
                    refurl = 'http://graspa.oapd.inaf.it/cgi-bin/sncat.php'
                    secondarysource = add_source(name, refname = reference, url = refurl, secondary = True)
                    add_quantity(name, 'alias', oldname, secondarysource)
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
                        source = add_source(name, refname = reference, url = refurl)
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
        for folder in tq(sorted(next(os.walk("../sne-external-WISEREP"))[1], key=lambda s: s.lower()), currenttask):
            files = glob("../sne-external-WISEREP/" + folder + '/*')
            for fname in tq(files, currenttask):
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
                                            except (KeyboardInterrupt, SystemExit):
                                                raise
                                            except StopIteration:
                                                pass
                                            #if not specpath:
                                            #    warnings.warn('Spectrum file not found, "' + specfile + '"')
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
    
                                #print(name + " " + claimedtype + " " + epoch + " " + observer + " " + reducer + " " + specfile + " " + bibcode + " " + redshift)
    
                                (name, secondarysource) = new_event(name, refname = secondaryreference, url = secondaryrefurl, bibcode = secondarybibcode, secondary = True)
                                if bibcode:
                                    newbibcode = bibcode
                                    if bibcode in wiserepbibcorrectdict:
                                        newbibcode = wiserepbibcorrectdict[bibcode]
                                    if newbibcode:
                                        source = add_source(name, bibcode = unescape(newbibcode))
                                    else:
                                        source = add_source(name, refname = unescape(bibcode))
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
        for name in tq(sorted(next(os.walk("../sne-external-spectra/CfA_SNIa"))[1], key=lambda s: s.lower()), currenttask):
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
            reference = 'CfA Supernova Archive'
            refurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
            (name, source) = new_event(name, refname = reference, url = refurl, secondary = True, acknowledgment = cfaack)
            for fi, fname in enumerate(sorted(glob(fullpath + '/*'), key=lambda s: s.lower())):
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
                time = str(astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day)))
                f = open(fname,'r')
                data = csv.reader(f, delimiter=' ', skipinitialspace=True)
                data = [list(i) for i in zip(*data)]
                wavelengths = data[0]
                fluxes = data[1]
                errors = data[2]
                sources = uniq_cdl([source, add_source(name, bibcode = '2012AJ....143..126B'), add_source(name, bibcode = '2008AJ....135.1598M')])
                add_spectrum(name = name, waveunit = 'Angstrom', fluxunit = 'erg/s/cm^2/Angstrom', filename = filename,
                    wavelengths = wavelengths, fluxes = fluxes, u_time = 'MJD' if time else '', time = time, instrument = instrument,
                    errorunit = "ergs/s/cm^2/Angstrom", errors = errors, source = sources, dereddened = False, deredshifted = False)
                if args.travis and fi >= travislimit:
                    break
        journal_events()
    
        # Ibc spectra
        oldname = ''
        for name in tq(sorted(next(os.walk("../sne-external-spectra/CfA_SNIbc"))[1], key=lambda s: s.lower()), currenttask):
            fullpath = "../sne-external-spectra/CfA_SNIbc/" + name
            if name.startswith('sn') and is_number(name[2:6]):
                name = 'SN' + name[2:]
            name = get_preferred_name(name)
            if oldname and name != oldname:
                journal_events()
            oldname = name
            reference = 'CfA Supernova Archive'
            refurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
            (name, source) = new_event(name, refname = reference, url = refurl, secondary = True, acknowledgment = cfaack)
            for fi, fname in enumerate(sorted(glob(fullpath + '/*'), key=lambda s: s.lower())):
                filename = os.path.basename(fname)
                fileparts = filename.split('-')
                instrument = ''
                year = fileparts[1][:4]
                month = fileparts[1][4:6]
                day = fileparts[1][6:].split('.')[0]
                if len(fileparts) > 2:
                    instrument = fileparts[-1].split('.')[0]
                time = str(astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day)))
                f = open(fname,'r')
                data = csv.reader(f, delimiter=' ', skipinitialspace=True)
                data = [list(i) for i in zip(*data)]
                wavelengths = data[0]
                fluxes = data[1]
                sources = uniq_cdl([source, add_source(name, bibcode = '2014AJ....147...99M')])
                add_spectrum(name = name, waveunit = 'Angstrom', fluxunit = 'erg/s/cm^2/Angstrom', wavelengths = wavelengths, filename = filename,
                    fluxes = fluxes, u_time = 'MJD' if time else '', time = time, instrument = instrument, source = sources,
                    dereddened = False, deredshifted = False)
                if args.travis and fi >= travislimit:
                    break
        journal_events()

        # Other spectra
        oldname = ''
        for name in tq(sorted(next(os.walk("../sne-external-spectra/CfA_Extra"))[1], key=lambda s: s.lower()), currenttask):
            fullpath = "../sne-external-spectra/CfA_Extra/" + name
            if name.startswith('sn') and is_number(name[2:6]):
                name = 'SN' + name[2:]
            name = get_preferred_name(name)
            if oldname and name != oldname:
                journal_events()
            oldname = name
            reference = 'CfA Supernova Archive'
            refurl = 'https://www.cfa.harvard.edu/supernova/SNarchive.html'
            (name, source) = new_event(name, refname = reference, url = refurl, secondary = True, acknowledgment = cfaack)
            for fi, fname in enumerate(sorted(glob(fullpath + '/*'), key=lambda s: s.lower())):
                if not os.path.isfile(fname):
                    continue
                filename = os.path.basename(fname)
                if (not filename.startswith('sn') or not filename.endswith('flm') or
                    any(x in filename for x in ['-interp', '-z', '-dered', '-obj', '-gal'])):
                    continue
                fileparts = filename.split('.')[0].split('-')
                instrument = ''
                time = ''
                if len(fileparts) > 1:
                    year = fileparts[1][:4]
                    month = fileparts[1][4:6]
                    day = fileparts[1][6:]
                    if is_number(year) and is_number(month) and is_number(day):
                        if len(fileparts) > 2:
                            instrument = fileparts[-1]
                        time = str(astrotime(year + '-' + month + '-' + str(floor(float(day))).zfill(2)).mjd + float(day) - floor(float(day)))
                f = open(fname,'r')
                data = csv.reader(f, delimiter=' ', skipinitialspace=True)
                data = [list(i) for i in zip(*data)]
                wavelengths = data[0]
                fluxes = [str(Decimal(x)*Decimal(1.0e-15)) for x in data[1]]
                add_spectrum(name = name, waveunit = 'Angstrom', fluxunit = 'erg/s/cm^2/Angstrom', wavelengths = wavelengths, filename = filename,
                    fluxes = fluxes, u_time = 'MJD' if time else '', time = time, instrument = instrument, source = source,
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
        for fi, fname in enumerate(tq(sorted(glob('../sne-external-spectra/SNLS/*'), key=lambda s: s.lower()), currenttask = currenttask)):
            filename = os.path.basename(fname)
            fileparts = filename.split('_')
            name = 'SNLS-' + fileparts[1]
            name = get_preferred_name(name)
            if oldname and name != oldname:
                journal_events()
            oldname = name
            (name, source) = new_event(name, bibcode = "2009A&A...507...85B")
    
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
        for fi, fname in enumerate(tq(sorted(glob('../sne-external-spectra/CSP/*'), key=lambda s: s.lower()), currenttask = currenttask)):
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
            telescope = fileparts[-2]
            instrument = fileparts[-1]
            (name, source) = new_event(name, bibcode = "2013ApJ...773...53F")
    
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
        for spectrum in tq(spectra, currenttask = currenttask):
            name = spectrum["ObjName"]
            if oldname and name != oldname:
                journal_events()
            oldname = name
    
            (name, secondarysource) = new_event(name, refname = secondaryreference, url = secondaryrefurl, bibcode = secondaryrefbib, secondary = True)
            sources = [secondarysource]
            if spectrum["Reference"]:
                sources += [add_source(name, bibcode = spectrum["Reference"])]
            sources = uniq_cdl(sources)
    
            if spectrum["Type"] and spectrum["Type"].strip() != "NoMatch":
                for ct in spectrum["Type"].strip().split(','):
                    add_quantity(name, 'claimedtype', ct.replace('-norm', '').strip(), sources)
            if spectrum["DiscDate"]:
                add_quantity(name, 'discoverdate', spectrum["DiscDate"].replace('-', '/'), sources)
            if spectrum["HostName"]:
                add_quantity(name, 'host', urllib.parse.unquote(spectrum["HostName"]).replace('*', ''), sources)
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
            if archived_task('ucbspectra') and os.path.isfile(filepath):
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
        for folder in tq(folders, currenttask):
            eventfolders = next(os.walk('../sne-external-spectra/Suspect/'+folder))[1]
            oldname = ''
            for eventfolder in tq(eventfolders, currenttask):
                name = eventfolder
                if is_number(name[:4]):
                    name = 'SN' + name
                name = get_preferred_name(name)
                if oldname and name != oldname:
                    journal_events()
                oldname = name
                secondaryreference = "SUSPECT"
                secondaryrefurl = "https://www.nhn.ou.edu/~suspect/"
                secondarybibcode = "2001AAS...199.8408R"
                (name, secondarysource) = new_event(name, refname = secondaryreference, url = secondaryrefurl, bibcode = secondarybibcode, secondary = True)
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
                        sources += [source]
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
            secondaryreference = "Nearby Supernova Factory"
            secondaryrefurl = "http://snfactory.lbl.gov/"
            secondarybibcode = "2002SPIE.4836...61A"
            (name, secondarysource) = new_event(name, refname = secondaryreference, url = secondaryrefurl, bibcode = secondarybibcode, secondary = True)
            bibcode = bibcodes[name]
            source = add_source(name, bibcode = bibcode)
            sources = uniq_cdl([source,secondarysource])
            eventspectra = glob('../sne-external-spectra/SNFactory/'+eventfolder+'/*.dat')
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
        sfdirs = glob('../sne-external-spectra/superfit/*')
        for sfdir in tq(sfdirs, currenttask = currenttask):
            sffiles = sorted(glob(sfdir + "/*.dat"))
            lastname = ''
            oldname = ''
            for sffile in tq(sffiles, currenttask = currenttask):
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
    
                source = add_source(name, refname = 'Superfit', url = 'http://www.dahowell.com/superfit.html', secondary = True)
                add_quantity(name, 'alias', oldname, source)
    
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

    if do_task(task, 'mergeduplicates'):
        if args.update and not len(events):
            tprint('No sources changed, event files unchanged in update.')
            sys.exit(1)
        merge_duplicates()

    if do_task(task, 'setprefnames'):
        set_preferred_names()

files = repo_file_list()

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

for i, fi in enumerate(tq(files, 'Sanitizing and deriving quantities for events')):
    events = OrderedDict()
    name = os.path.basename(os.path.splitext(fi)[0]).replace('.json', '')
    name = add_event(name, loadifempty = False)
    derive_and_sanitize()
    if has_task('writeevents'): 
        write_all_events(empty = True, gz = True, bury = True)
    if args.travis and i > travislimit:
        break

jsonstring = json.dumps(bibauthordict, indent='\t', separators=(',', ':'), ensure_ascii=False)
with codecs.open('../bibauthors.json', 'w', encoding='utf8') as f:
    f.write(jsonstring)
jsonstring = json.dumps(extinctionsdict, indent='\t', separators=(',', ':'), ensure_ascii=False)
with codecs.open('../extinctions.json', 'w', encoding='utf8') as f:
    f.write(jsonstring)

print("Memory used (MBs on Mac, GBs on Linux): " + "{:,}".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024./1024.))

sys.exit(0)
