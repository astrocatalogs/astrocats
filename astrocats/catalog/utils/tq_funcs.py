from datetime import datetime

from tqdm import tqdm

__all__ = ['tq', 'tprint', 'pbar', 'pbar_strings']


def tq(li, currenttask='', leave=True):
    return tqdm(
        list(li),
        desc=('<' + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + '> ' +
              currenttask),
        leave=leave,
        dynamic_ncols=True)


def pbar(iter, desc='', **kwargs):
    """Wrapper for `tqdm` progress bar.
    """
    return tqdm(
        iter,
        desc=('<' + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + '> ' +
              desc),
        dynamic_ncols=True,
        **kwargs)


def pbar_strings(files, desc='', **kwargs):
    """Wrapper for `tqdm` progress bar which also sorts list of strings
    """
    return tqdm(
        sorted(files, key=lambda s: s.lower()),
        desc=('<' + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + '> ' +
              desc),
        dynamic_ncols=True,
        **kwargs)


def tprint(string):
    """Print string via `tqdm` so that it doesnt interfere with a progressbar.
    """
    try:
        tqdm.write(string)
    except:
        print(string)
