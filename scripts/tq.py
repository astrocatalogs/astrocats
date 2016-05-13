from tqdm import tqdm, trange

currenttask = ''

def tq(li, leave = True):
    global currenttask
    return tqdm(list(li), desc = currenttask, leave = leave)

def tprint(string):
    try:
        tqdm.write(string)
    except:
        print(string)
