"""Import tasks for the PESSTO spectroscopic program.
"""
import csv
import os

from scripts import PATH


def do_pessto(catalog):
    pessto_path = os.path.join(PATH.REPO_EXTERNAL, 'PESSTO_MPHOT.csv')
    tsvin = list(csv.reader(open(pessto_path, 'r'), delimiter=','))
    for ri, row in enumerate(tsvin):
        if ri == 0:
            bands = [xx.split('_')[0] for xx in row[3::2]]
            systems = [xx.split('_')[1].capitalize().replace(
                'Ab', 'AB') for xx in row[3::2]]
            continue
        name = row[1]
        name = catalog.add_entry(name)
        source = catalog.events[name].add_source(bibcode='2015A&A...579A..40S')
        catalog.events[name].add_quantity('alias', name, source)
        for hi, ci in enumerate(range(3, len(row) - 1, 2)):
            if not row[ci]:
                continue
            teles = 'Swift' if systems[hi] == 'Swift' else ''
            (catalog.events[name]
             .add_photometry(time=row[2], magnitude=row[ci],
                             e_magnitude=row[ci + 1],
                             band=bands[hi], system=systems[hi],
                             telescope=teles,
                             source=source))

    catalog.journal_events()
    return
