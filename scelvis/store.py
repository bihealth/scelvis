"""Access to all data outside of the session.

This module uses the functions from ``data`` to load the data from the appropriate location
and memoizes the loaded data in the Flask cacche.  This module is aware of the Flask cache and
thus also of the Dash app.
"""

import contextlib

import fs.path
from logzero import logger

from . import data, settings
from .exceptions import ScelVisException
from .cache import cache


@contextlib.contextmanager
def list_identifiers(data_source):
    """Context manager for listing all identifiers inside the given ``data_source`` specification."""


@cache.memoize()
def does_exist(url, path, *more_components):
    """Return whether the given path exists behind the given URL."""
    if url.scheme in data.PYFS_SCHEMES:
        curr_fs = data.make_fs(url)
        return curr_fs.exists(fs.path.join(path, *more_components))
    elif url.scheme == "irods":
        pass  # TODO: implement me!
    else:
        raise ScelVisException("Invalid URL scheme: %s" % url.scheme)


@cache.memoize()
def glob_data_sets(url):
    """Return list of all data sets behind the given ``url``."""
    result = []
    if url.scheme in data.PYFS_SCHEMES:
        curr_fs = data.make_fs(url)
        for match in curr_fs.glob(fs.path.join("*", settings.ABOUT_FILENAME)):
            match_path = fs.path.basename(match.path)
            logger.info("Found data set %s at %s" % (match_path, data.redacted_urlunparse(url)))
            result.append(url._replace(path=fs.path.join(url.path, match.path)))
    elif url.scheme == "irods":
        pass  # TODO: implement me!
    else:
        raise ScelVisException("Invalid URL scheme: %s" % url.scheme)
    return result


@cache.memoize()
def load_all_metadata():
    """Load all meta data information from ``settings.DATA_SOURCES``.

    If ``data_source`` itself contains a file called ``about.md``, it is assumed that only one dataset is available.
    Otherwise, all sub directories will be scanned for ``about.md`` and one entry is returned for each dataset.
    """
    result = []
    for data_source in settings.DATA_SOURCES:
        if does_exist(data_source, settings.ABOUT_FILENAME):
            logger.info(
                "Loading single dataset from data source %s", data.redacted_urlunparse(data_source)
            )
            result.append(data.load_metadata(data_source))
        else:
            lst = []
            for match in glob_data_sets(data_source):
                identifier = fs.path.basename(fs.path.dirname(match.path))
                lst.append(data.load_metadata(data_source, identifier))
            if lst:
                logger.info("Loaded %d data sets from data directory.", len(lst))
            else:
                logger.warn(
                    "No data sets found in data directory %s", data.redacted_urlunparse(data_source)
                )
            result += lst
    return result


@cache.memoize()
def load_metadata(identifier):
    """Load metadata for the given identifier from data or upload directory."""
    for data_source in settings.DATA_SOURCES + [settings.UPLOAD_DIR]:
        if data_source and does_exist(data_source, identifier, settings.ABOUT_FILENAME):
            return data.load_metadata(data_source, identifier)


@cache.memoize()
def load_data(identifier):
    """Load data for the given identifier from data or upload directory."""
    for data_source in settings.DATA_SOURCES + [settings.UPLOAD_DIR]:
        if data_source and does_exist(data_source, identifier, settings.ABOUT_FILENAME):
            return data.load_data(data_source, identifier)
