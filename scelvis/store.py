"""Access to all data outside of the session.

This module uses the functions from ``data`` to load the data from the appropriate location
and memoizes the loaded data in the Flask cacche.  This module is aware of the Flask cache and
thus also of the Dash app.

Note well that the behaviour of iRODS is special because of how ticket access is implemented.
We have to get direct access to the collection for which we have a ticket.  Then we need to
perform a linear search for the collection's data objects and only THEN can we open them.
"""

import urllib.parse

import fs.path
from irods.exception import CAT_SQL_ERR, DoesNotExist
from logzero import logger

from scelvis.data import redacted_urlunparse
from . import data, settings
from .exceptions import ScelVisException
from .cache import cache


@cache.memoize()
def does_exist(url, path, *more_components):
    """Return whether the given path exists behind the given URL."""
    path_full = fs.path.join(url.path, path, *more_components)
    logger.info("Checking whether %s exists in %s", path_full, redacted_urlunparse(url))
    if url.scheme in data.PYFS_SCHEMES:
        curr_fs = data.make_fs(url)
        result = curr_fs.exists(path_full)
        return result
    elif url.scheme.startswith("irods"):
        with data.create_irods_session(url) as irods_session:
            path_collection = fs.path.dirname(path_full)
            name = fs.path.basename(path_full)
            try:
                collection = irods_session.collections.get(path_collection)
                for data_object in collection.data_objects:
                    if data_object.name == name:
                        return True
            except (CAT_SQL_ERR, DoesNotExist) as e:
                pass
            logger.info("=> False")
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
    elif url.scheme.startswith("irods"):
        with data.create_irods_session(url) as irods_session:
            # Get pointed-to collection.
            collection = irods_session.collections.get(url.path)
            for sub_coll in collection.subcollections:
                for data_obj in sub_coll.data_objects:
                    if data_obj.name == settings.ABOUT_FILENAME:
                        result.append(url._replace(path=fs.path.join(url.path, sub_coll.name, data_obj.name)))
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
    for url in settings.DATA_SOURCES:
        if does_exist(url, settings.ABOUT_FILENAME):
            logger.info(
                "Loading single dataset from data source %s", data.redacted_urlunparse(url)
            )
            identifier = fs.path.basename(url.path)
            collection = url._replace(path=fs.path.dirname(url.path))
            result.append(data.load_metadata(collection, identifier))
        else:
            lst = []
            for match in glob_data_sets(url):
                identifier = fs.path.basename(fs.path.dirname(match.path))
                lst.append(data.load_metadata(url, identifier))
            if lst:
                logger.info("Loaded %d data sets from data directory.", len(lst))
            else:
                logger.warn(
                    "No data sets found in data directory %s", data.redacted_urlunparse(url)
                )
            result += lst
    return result


def _load(identifier, load_func):
    urls = settings.DATA_SOURCES
    if settings.UPLOAD_DIR:
        urls += [urllib.parse.urlparse("file://%s" % settings.UPLOAD_DIR)]
    for url in urls:
        if fs.path.basename(url.path) == identifier:
            url = url._replace(path=fs.path.dirname(url.path))
            return load_func(url, identifier)
        else:
            if does_exist(url, identifier, settings.ABOUT_FILENAME):
                return load_func(url, identifier)


@cache.memoize()
def load_metadata(identifier):
    """Load metadata for the given identifier from data or upload directory."""
    return _load(identifier, data.load_metadata)


@cache.memoize()
def load_data(identifier):
    """Load data for the given identifier from data or upload directory."""
    return _load(identifier, data.load_data)
