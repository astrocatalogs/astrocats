import csv
import os
from glob import glob

from astropy.time import Time as astrotime
from astroquery.vizier import Vizier

from scripts import PATH

from .. import Events
from ...utils import get_sig_digits, pbar_strings, pretty_num
from ..constants import TRAVIS_QUERY_LIMIT
from ..funcs import add_spectrum, get_preferred_name


def do_snls_spectra(events, stubs, args, tasks, task_obj, log):
    """
    """

    current_task = task_obj.current_task(args)
    result = Vizier.get_catalogs('J/A+A/507/85/table1')
    table = result[list(result.keys())[0]]
    table.convert_bytestring_to_unicode(python3_only=True)
    datedict = {}
    for row in table:
        datedict['SNLS-' + row['SN']] = str(astrotime(row['Date']).mjd)

    oldname = ''
    file_names = glob(os.path.join(PATH.REPO_EXTERNAL_SPECTRA, 'SNLS/*'))
    for fi, fname in enumerate(pbar_strings(file_names,
                                            current_task=current_task)):
        filename = os.path.basename(fname)
        fileparts = filename.split('_')
        name = 'SNLS-' + fileparts[1]
        name = get_preferred_name(events, name)
        if oldname and name != oldname:
            events, stubs = Events.journal_events(
                tasks, args, events, stubs, log)
        oldname = name
        events, name = Events.add_event(tasks, args, events, name, log)
        source = events[name].add_source(bibcode='2009A&A...507...85B')
        events[name].add_quantity('alias', name, source)

        events[name].add_quantity(
            'discoverdate', '20' + fileparts[1][:2], source)

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

        fluxes = [pretty_num(float(x) * 1.e-16, sig=get_sig_digits(x))
                  for x in specdata[2]]
        # FIX: this isnt being used
        # errors = [pretty_num(float(x)*1.e-16, sig=get_sig_digits(x)) for x in specdata[3]]

        add_spectrum(
            events, name, 'Angstrom', 'erg/s/cm^2/Angstrom', wavelengths=wavelengths,
            fluxes=fluxes, u_time='MJD' if name in datedict else '',
            time=datedict[name] if name in datedict else '', telescope=telescope, source=source,
            filename=filename)
        if args.travis and fi >= TRAVIS_QUERY_LIMIT:
            break
    events, stubs = Events.journal_events(tasks, args, events, stubs, log)
    return events
