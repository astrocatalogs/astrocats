"""
"""

__all__ = ['rep_chars']

def rep_chars(string, chars, rep = ''):
    for c in chars:
        if c in string:
            string = string.replace(c, rep)
    return string

def single_spaces(string):
    return ' '.join(list(filter(None, string.split())))
