"""
"""
from collections import OrderedDict
import json
import os

from scripts import FILENAME
from .funcs import copy_to_event, get_aliases, journal_events, load_event_from_file, load_stubs
from ..utils import is_number, pbar, tprint


def merge_duplicates(tasks, args, events):
    """Merge and remove duplicate events
    """
    if len(events) == 0:
        if args.update:
            import sys
            tprint('No sources changed, event files unchanged in update.')
            sys.exit(1)

        load_stubs(tasks, args, events)

    currenttask = 'Merging duplicate events'
    keys = list(sorted(list(events.keys())))
    for n1, name1 in enumerate(pbar(keys[:], currenttask)):
        if name1 not in events:
            continue
        # allnames1 = get_aliases(events, name1) + (['AT' + name1[2:]] if
        #     (name1.startswith('SN') and is_number(name1[2:6])) else [])
        allnames1 = get_aliases(events, name1)
        if name1.startswith('SN') and is_number(name1[2:6]):
            allnames1 += ['AT' + name1[2:]]

        for name2 in keys[n1+1:]:
            if name2 not in events or name1 == name2:
                continue
            # allnames2 = get_aliases(events, name2) + (['AT' + name2[2:]]
            #    if (name2.startswith('SN') and is_number(name2[2:6])) else [])
            allnames2 = get_aliases(events, name2)
            if name2.startswith('SN') and is_number(name2[2:6]):
                allnames2 += ['AT' + name2[2:]]
            if any(ii in allnames1 for ii in allnames2):
                tprint("Found single event with multiple entries ('{}' and '{}'), merging.".format(
                    name1, name2))

                load1 = load_event_from_file(events, args, tasks, name1, delete=True)
                load2 = load_event_from_file(events, args, tasks, name2, delete=True)
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
                        copy_to_event(events, name2, name1)
                        keys.append(name1)
                        del(events[name2])
                    else:
                        copy_to_event(events, name1, name2)
                        keys.append(name2)
                        del(events[name1])
                else:
                    print('Duplicate already deleted')
                journal_events(tasks, args, events)

    return events


def set_preferred_names(tasks, args, events):
    if not len(events):
        load_stubs(tasks, args, events)

    for name in list(sorted(list(events.keys()))):
        if name not in events:
            continue
        newname = ''
        aliases = get_aliases(events, name)
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
            if load_event_from_file(events, args, tasks, newname):
                continue
            if load_event_from_file(events, args, tasks, name, delete=True):
                tprint('Changing event name (' + name + ') to preferred name (' + newname + ').')
                events[newname] = events[name]
                events[newname]['name'] = newname
                del(events[name])
                journal_events(tasks, args, events)

    return events


def write_all_events(events, args, empty=False, gz=False, bury=False):
    """Save all `events` to files.
    """
    import codecs
    from scripts import PATH
    from .. utils import get_event_filename, get_repo_folders, get_repo_years, is_number, tprint
    repo_folders = get_repo_folders()
    non_sne_types = None
    if bury:
        with open(FILENAME.NON_SNE_TYPES, 'r') as f:
            non_sne_types = json.loads(f.read(), object_pairs_hook=OrderedDict)
            non_sne_types = [x.upper() for x in non_sne_types]

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

        # outdir = '../'
        outdir = str(PATH.ROOT)
        if 'discoverdate' in events[name]:
            repo_years = get_repo_years(repo_folders)
            for r, year in enumerate(repo_years):
                if int(events[name]['discoverdate'][0]['value'].split('/')[0]) <= year:
                    # outdir += repo_folders[r]
                    outdir = os.path.join(outdir, repo_folders[r])
                    break
        else:
            # outdir += str(repo_folders[0])
            outdir = os.path.join(outdir, repo_folders[0])

        # Delete non-SN events here without IAU designations (those with only banned types)
        if bury:
            buryevent = False
            nonsneprefixes = ('PNVJ', 'PNV J', 'OGLE-2013-NOVA')
            if name.startswith(nonsneprefixes):
                tprint('Burying ' + name + ', non-SNe prefix.')
                continue
            if 'claimedtype' in events[name] and not (name.startswith('SN') and is_number(name[2:6])):
                for ct in events[name]['claimedtype']:
                    if ct['value'].upper() not in non_sne_types and ct['value'].upper() != 'CANDIDATE':
                        buryevent = False
                        break
                    if ct['value'].upper() in non_sne_types:
                        buryevent = True
                if buryevent:
                    tprint('Burying ' + name + ' (' + ct['value'] + ').')
                    # outdir = '../sne-boneyard'
                    outdir = str(PATH.REPO_BONEYARD)  # os.path.join(PATH.ROOT, 'sne-boneyard')

        jsonstring = json.dumps({name: events[name]}, indent='\t', separators=(',', ':'), ensure_ascii=False)

        # path = outdir + '/' + filename + '.json'
        path = os.path.join(outdir, filename + '.json')
        with codecs.open(path, 'w', encoding='utf8') as f:
            f.write(jsonstring)

        if gz:
            if os.path.getsize(path) > 90000000:
                import shutil
                import gzip
                if not args.travis:
                    tprint('Compressing ' + name)
                with open(path, 'rb') as f_in, gzip.open(path + '.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                os.remove(path)
                os.system('cd ' + outdir + '; git rm ' + filename + '.json; git add -f ' + filename + '.json.gz; cd ' + '../scripts')

    return
