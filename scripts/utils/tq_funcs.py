from tqdm import tqdm, trange

__all__ = ['tq', 'tprint']


def tq(li, currenttask='', leave=True):
    return tqdm(list(li), desc=currenttask, leave=leave)


def tprint(string):
    try:
        tqdm.write(string)
    except:
        print(string)
