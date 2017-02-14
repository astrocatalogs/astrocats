"""
"""

__all__ = [
    'listify'
]


def listify(x):
    if not isinstance(x, list):
        return [x]
    return x
