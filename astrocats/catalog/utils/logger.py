"""Logging submodule and related functions.
"""

import inspect
import logging
from logging import DEBUG, INFO, WARNING

import numpy as np

_FILE_LEVEL_DEF = DEBUG
_STREAM_LEVEL_DEF = WARNING
_LOADED_LEVEL = INFO


__all__ = ["get_logger", "log_raise", "DEBUG", "WARNING", "INFO", "log_memory"]


class IndentFormatter(logging.Formatter):
    """Logging formatter where the depth of the stack sets the message
    indentation level.
    """

    def __init__(self, fmt=None, datefmt=None):
        logging.Formatter.__init__(self, fmt, datefmt)
        self.baseline = None

    def format(self, rec):
        stack = inspect.stack()
        if self.baseline is None:
            self.baseline = len(stack)
        indent = (len(stack) - self.baseline)
        addSpace = ((indent > 0) & (not rec.msg.startswith(" -")))
        rec.indent = ' -' * indent + ' ' * addSpace
        out = logging.Formatter.format(self, rec)
        del rec.indent
        return out


def get_logger(name=None, stream_fmt=None, file_fmt=None, date_fmt=None,
               stream_level=None, file_level=None,
               tofile=None, tostr=True):
    """Create a standard logger object which logs to file and or stdout stream.

    If a logger has already been created in this session, it is returned
    (unless `name` is given).

    Arguments
    ---------
    name : str,
        Handle for this logger, must be distinct for a distinct logger.
    stream_fmt : str or `None`,
        Format of log messages to stream (stdout).  If `None`, default settings
        are used.
    file_fmt : str or `None`,
        Format of log messages to file.  If `None`, default settings are used.
    date_fmt : str or `None`
        Format of time stamps to stream and/or file.  If `None`, default
        settings are used.
    stream_level : int,
        Logging level for stream.
    file_level : int,
        Logging level for file.
    tofile : str or `None`,
        Filename to log to (turned off if `None`).
    tostr : bool,
        Log to stdout stream.

    Returns
    -------
    logger : ``logging.Logger`` object,
        Logger object to use for logging.

    """
    if tofile is None and not tostr:
        raise ValueError(
            "Must log to something: `tofile` or `tostr` must be `True`.")

    logger = logging.getLogger(name)
    # Add a custom attribute to this `logger` so that we know when an existing
    # one is being returned
    if hasattr(logger, '_OSC_LOGGER'):
        return logger
    else:
        logger._OSC_LOGGER = True

    # Set other custom parameters
    logger._LOADED = _LOADED_LEVEL

    # Make sure handlers don't get duplicated (ipython issue)
    while len(logger.handlers) > 0:
        logger.handlers.pop()
    # Prevents duplication or something something...
    logger.propagate = 0

    # Determine and Set Logging Levels
    if file_level is None:
        file_level = _FILE_LEVEL_DEF
    if stream_level is None:
        stream_level = _STREAM_LEVEL_DEF
    # Logger object must be at minimum level
    logger.setLevel(int(np.min([file_level, stream_level])))

    if date_fmt is None:
        date_fmt = '%Y/%m/%d %H:%M:%S'

    # Log to file
    # -----------
    if tofile is not None:
        if file_fmt is None:
            file_fmt = "%(asctime)s %(levelname)8.8s [%(filename)20.20s:"
            file_fmt += "%(funcName)-20.20s]%(indent)s%(message)s"

        fileFormatter = IndentFormatter(file_fmt, datefmt=date_fmt)
        fileHandler = logging.FileHandler(tofile, 'w')
        fileHandler.setFormatter(fileFormatter)
        fileHandler.setLevel(file_level)
        logger.addHandler(fileHandler)
        #     Store output filename to `logger` object
        logger.filename = tofile

    # Log To stdout
    # -------------
    if tostr:
        if stream_fmt is None:
            stream_fmt = "%(indent)s%(message)s"

        strFormatter = IndentFormatter(stream_fmt, datefmt=date_fmt)
        strHandler = logging.StreamHandler()
        strHandler.setFormatter(strFormatter)
        strHandler.setLevel(stream_level)
        logger.addHandler(strHandler)

    return logger


def log_raise(log, err_str, err_type=RuntimeError):
    """Log an error message and raise an error.

    Arguments
    ---------
    log : `logging.Logger` object
    err_str : str
        Error message to be logged and raised.
    err_type : `Exception` object
        Type of error to raise.

    """
    log.error(err_str)
    # Make sure output is flushed
    # (happens automatically to `StreamHandlers`, but not `FileHandlers`)
    for handle in log.handlers:
        handle.flush()
    # Raise given error
    raise err_type(err_str)


def log_memory(log, pref=None, lvl=logging.DEBUG, raise_flag=True):
    """Log the current memory usage.
    """
    import os
    import sys
    cyc_str = ""
    KB = 1024.0
    if pref is not None:
        cyc_str += "{}: ".format(pref)

    # Linux returns units in Bytes; OSX in kilobytes
    UNIT = KB*KB if sys.platform == 'darwin' else KB

    good = False
    # Use the `resource` module to check the maximum memory usage of this process
    try:
        import resource
        max_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        max_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        _str = "RSS Max Self: {:7.2f} [MB], Child: {:7.2f} [MB]".format(
            max_self/UNIT, max_child/UNIT)
        cyc_str += _str
    except Exception as err:
        log.log(lvl, "resource.getrusage failed.  '{}'".format(str(err)))
        if raise_flag:
            raise
    else:
        good = True

    # Use the `psutil` module to check the current memory/cpu usage of this process
    try:
        import psutil
        process = psutil.Process(os.getpid())
        rss = process.memory_info().rss
        cpu_perc = process.cpu_percent()
        mem_perc = process.memory_percent()
        num_thr = process.num_threads()
        _str = "; RSS: {:7.2f} [MB], {:7.2f}%; Threads: {:3d}, CPU: {:7.2f}%".format(
            rss/UNIT, mem_perc, num_thr, cpu_perc)
        cyc_str += _str
    except Exception as err:
        log.log(lvl, "psutil.Process failed.  '{}'".format(str(err)))
        if raise_flag:
            raise
    else:
        good = True

    if good:
        log.log(lvl, cyc_str)

    return
