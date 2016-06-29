from tqdm import tqdm

__all__ = ['tq', 'tprint', 'pbar', 'pbar_strings']


def tq(li, currenttask='', leave=True):
    return tqdm(list(li), desc=currenttask, leave=leave)


def pbar(iter, desc='', **kwargs):
    """Wrapper for `tqdm` progress bar.
    """
    return tqdm(iter, desc=desc, **kwargs)


def pbar_strings(files, desc='', **kwargs):
    """Wrapper for `tqdm` progress bar which also sorts list of strings
    """
    return tqdm(sorted(files, key=lambda s: s.lower()), desc=desc, **kwargs)


def tprint(string):
    """Print string via `tqdm` so that it doesnt interfere with a progressbar.
    """
    try:
        tqdm.write(string)
    except:
        print(string)
