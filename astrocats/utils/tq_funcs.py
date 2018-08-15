from datetime import datetime

from tqdm import tqdm

__all__ = ['tprint', 'pbar']


'''
def pbar(li, currenttask='', leave=True):
    desc = '<' + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + '> ' + currenttask
    return tqdm(list(li), desc=desc, leave=leave, dynamic_ncols=True)
'''


def pbar(iters, desc='', sort=False, **kwargs):
    """Wrapper for `tqdm` progress bar.
    """
    kwargs.setdefault('leave', True)
    kwargs.setdefault('dynamic_ncols', True)
    if sort:
        iters = sorted(iters, key=lambda s: s.lower())
    desc = '<' + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + '> ' + desc
    return tqdm(iters, desc=desc, **kwargs)


'''
def pbar_strings(files, desc='', **kwargs):
    """Wrapper for `tqdm` progress bar which also sorts list of strings
    """
    return tqdm(
        sorted(files, key=lambda s: s.lower()),
        desc=('<' + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + '> ' +
              desc),
        dynamic_ncols=True,
        **kwargs)
'''


def tprint(string):
    """Print string via `tqdm` so that it doesnt interfere with a progressbar.
    """
    try:
        tqdm.write(string)
    except:
        print(string)
