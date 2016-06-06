import gzip

__all__ = ['get_event_filename', 'get_event_text']


def get_event_filename(name):
    return(name.replace('/', '_'))


def get_event_text(eventfile):
    if eventfile.split('.')[-1] == 'gz':
        with gzip.open(eventfile, 'rt') as f:
            filetext = f.read()
    else:
        with open(eventfile, 'r') as f:
            filetext = f.read()
    return filetext
