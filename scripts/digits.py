from math import log10, floor

def get_sig_digits(x):
    return len((''.join(x.split('.'))).strip('0'))

def round_sig(x, sig=4):
    if x == 0.0:
        return 0.0
    return round(x, sig-int(floor(log10(abs(x))))-1)

def pretty_num(x, sig=4):
    return str('%g'%(round_sig(x, sig)))

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
