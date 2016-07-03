"""Try to compare between the old and new OSC scripts.
"""
import copy
import os
from datetime import datetime
import json
from glob import glob
import numpy as np
from itertools import zip_longest

OLD_REPO_DIR = '../sne/'
NEW_REPO_DIR = './astrocats/supernovae/output/'
LIMIT_TRIES = 20
FAILURE_LIMIT = 2


def main():
    beg = datetime.now()
    start_str = "Starting {} at {}".format(__file__, beg.ctime())
    print("{}\n{}".format(start_str, '='*len(start_str)))
    good = 0
    bad = 0

    OLD_FILE = "../sne/sne-pre-1990/SDSS-II 4577.json"

    # Keep going until we reach the maximum number of failes
    while bad < FAILURE_LIMIT:
        # Choose a random file and compare it between old and new
        retval = check_random()
        if retval:
            good += 1
            return
        else:
            bad += 1
        tot = good + bad
        print("Good: {} ({:.4f}), Bad: {} ({:.4f})".format(
            good, good/tot, bad, bad/tot))

    end = datetime.now()
    print("Done at {} after {}".format(end.ctime(), end-beg))
    return


def check_random():
    # Get a list of the 'old' repositories
    old_repo_dirs = glob(os.path.join(OLD_REPO_DIR, 'sne-*'))
    # Choose a random output directory
    check_file = None
    count = 0
    # Try to find a json file to compare (give up after `LIMIT_TRIES` failures)
    while check_file is None and count < LIMIT_TRIES:
        # Choose one of the repo dirs at random
        check_dir = np.random.choice(old_repo_dirs)
        # Get a list of json files for this directory
        file_pattern = os.path.join(check_dir, '*.json')
        check_file = glob(file_pattern)
        # If no json files found, repeat
        if not check_file:
            check_file = None
        count += 1

    if count >= LIMIT_TRIES:
        raise RuntimeError("Couldn't find a `check_file`.  '{}'".format(
            file_pattern))

    # Choose a random json file from the list
    check_file = os.path.abspath(np.random.choice(check_file))
    return compare_files(check_file)


def compare_files(old_file):
    # Find the same file in the new script directories
    dir_name, file_name = os.path.split(old_file)
    new_dir = os.path.join(NEW_REPO_DIR, os.path.split(dir_name)[-1])
    new_file = os.path.join(new_dir, file_name)

    print("{} ====> {}".format(old_file, new_file))
    if not os.path.exists(new_file):
        print("NEW FILE DOES NOT EXIST")
        return False

    # Load json data into 'dicts' for each file version
    old_data = json.load(open(old_file, 'r'))
    new_data = json.load(open(new_file, 'r'))

    # Compare data recursively
    if not compare_dicts(old_data, new_data, copy.deepcopy(old_data), copy.deepcopy(new_data)):
        print("NEW FILE DOES NOT MATCH")
        print("\nOLD:")
        print(pprint(old_data))
        print("\nNEW:")
        print(pprint(new_data))
        print("\n")
        return False

    return True


def compare_dicts(old_full, new_full, old_data, new_data, depth=0):
    """Function compares dictionaries by key-value recursively.

    Old and new input data are both dictionaries
    """
    depth = depth + 1
    indent = "  "*depth

    # Print with an indentation matching the nested-dictionary depth
    def my_print(str):
        print("{}{}".format(indent, str))

    old_keys = list(old_data.keys())
    # Compare data key by key, in *this* dictionary level
    # Note: since we're comparing by keys explicity, order doesnt matter
    for key in old_keys:
        # Remove elements as we go
        old_vals = old_data.pop(key)
        # Current key
        my_print("{}".format(key))
        # If `new_data` doesnt also have this key, return False
        if key not in new_data:
            my_print("Key '{}' not in new_data.".format(key))
            my_print("Old:")
            my_print(pprint(new_data))
            my_print("New:")
            my_print(pprint(new_data))
            return False

        # If it does have the key, extract the values (remove as we go)
        new_vals = new_data.pop(key)
        # If these values are a sub-dictionary, compare those
        if isinstance(old_vals, dict) and isinstance(new_vals, dict):
            # If the sub-dictionary are not the same, return False
            if not compare_dicts(old_full, new_full, old_vals, new_vals, depth=depth):
                return False
        # If these values are a list of sub-dictionaries, compare each of those
        elif (isinstance(old_vals, list) and isinstance(old_vals[0], dict) and
              isinstance(old_vals, list) and isinstance(old_vals[0], dict)):
            for old_elem, new_elem in zip_longest(old_vals, new_vals):
                # If one or the other has extra elements, print message, but
                # continue on
                if old_elem is None or new_elem is None:
                    my_print("Missing element!")
                    my_print("\tOld: '{}'".format(old_elem))
                    my_print("\tNew: '{}'".format(new_elem))
                else:
                    if not compare_dicts(old_full, new_full, old_elem, new_elem, depth=depth):
                        return False

        # At the lowest-dictionary level, compare the values themselves
        else:
            # Turn everything into a list for convenience (most things should be
            # already)
            if  (not isinstance(old_vals, list) and
                 not isinstance(new_vals, list)):
                old_vals = [old_vals]
                new_vals = [new_vals]

            # Sort both lists
            old_vals = sorted(old_vals)
            new_vals = sorted(new_vals)

            for oldv, newv in zip_longest(old_vals, new_vals):
                # If one or the other has extra elements, print message, but
                # continue on
                if oldv is None or newv is None:
                    my_print("Missing element!")
                    my_print("\tOld: '{}'".format(oldv))
                    my_print("\tNew: '{}'".format(newv))
                # If values match, continue
                elif oldv == newv:
                    my_print("Good Match: '{}'".format(key))
                # If values dont match, return False
                else:
                    my_print("Bad  Match: '{}'".format(key))
                    my_print("\tOld: '{}'".format(oldv))
                    my_print("\tNew: '{}'".format(newv))
                    return False

    return True


def pprint(dat):
    return json.dumps(dat, sort_keys=True, indent=4, separators=(',', ': '))

if __name__ == "__main__":
    main()
