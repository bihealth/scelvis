"""Access to all data outside of the session.

This module uses the functions from ``data`` to load the data from the appropriate location
and memoizes the loaded data in the Flask cacche.  This module is aware of the Flask cache and
thus also of the Dash app.

Note well that the behaviour of iRODS is special because of how ticket access is implemented.
We have to get direct access to the collection for which we have a ticket.  Then we need to
perform a linear search for the collection's data objects and only THEN can we open them.
"""

import urllib.parse
from urllib.parse import urlunparse

import fs.path
from fs.errors import ResourceInvalid
import s3fs
from irods.exception import CAT_SQL_ERR, DoesNotExist
import htmllistparse
from logzero import logger
import requests
from requests.exceptions import HTTPError

from scelvis.data import redacted_urlunparse
from . import data, settings
from .exceptions import ScelVisException
from .cache import cache


@cache.memoize()
def does_exist(url, path, *more_components):
    """Return whether the given path exists behind the given URL."""
    logger.info(
        "Checking whether %s exists in %s",
        fs.path.join(path, *more_components),
        redacted_urlunparse(url),
    )
    if url.scheme in data.PYFS_SCHEMES:
        path_full = fs.path.join(url.path, path, *more_components)
        path_dir = fs.path.dirname(path_full)
        path_basename = fs.path.basename(path_full)
        try:
            curr_fs = data.make_fs(url._replace(path=path_dir))
        except ResourceInvalid:
            return False
        result = curr_fs.exists(path_basename)
        return result
    elif url.scheme == "s3":
        anon = url.username is None and url.password is None
        s3 = s3fs.S3FileSystem(anon=anon, key=url.username, secret=url.password)
        return s3.exists("%s/%s" % (url.hostname, path))
    elif url.scheme.startswith("http"):
        path_full = fs.path.join(url.path, path, *more_components)
        try:
            res = requests.get(urlunparse(url._replace(path=path_full)))
            res.raise_for_status()
        except HTTPError:
            return False
        return res.ok
    elif url.scheme.startswith("irods"):
        path_full = fs.path.join(url.path, path, *more_components)
        with data.create_irods_session(url) as irods_session:
            path_collection = fs.path.dirname(path_full)
            name = fs.path.basename(path_full)
            try:
                collection = irods_session.collections.get(path_collection)
                for data_object in collection.data_objects:
                    if data_object.name == name:
                        return True
            except (CAT_SQL_ERR, DoesNotExist):
                pass  # swallow
            logger.info("=> False")
    else:
        raise ScelVisException("Invalid URL scheme: %s" % url.scheme)


@cache.memoize()
def glob_data_sets(url):
    """Return list of all data sets behind the given ``url``."""
    result = []
    if url.scheme in data.PYFS_SCHEMES:
        curr_fs = data.make_fs(url)
        for match in curr_fs.glob("*.h5ad"):
            match_path = fs.path.basename(match.path)
            logger.info("Found data set %s at %s" % (match_path, data.redacted_urlunparse(url)))
            result.append(url._replace(path=fs.path.join(url.path, match.path[1:])))
    elif url.scheme == "s3":
        anon = url.username is None and url.password is None
        s3 = s3fs.S3FileSystem(anon=anon, key=url.username, secret=url.password)
        if url.path:
            pattern = "%s/%s/*.h5ad" % (url.hostname, url.path)
        else:
            pattern = "%s/*.h5ad" % (url.hostname,)
        for match in s3.glob(pattern):
            result.append(url._replace(path=match.split("/", 1)[1]))
    elif url.scheme.startswith("http"):
        cwd, listing = htmllistparse.fetch_listing(urlunparse(url), timeout=30)
        for entry in listing:
            if entry.name.endswith(".h5ad"):
                result.append(url._replace(path=fs.path.join(cwd, entry.name)))
    elif url.scheme.startswith("irods"):
        with data.create_irods_session(url) as irods_session:
            # Get pointed-to collection.
            collection = irods_session.collections.get(url.path)
            for data_obj in collection.data_objects:
                if data_obj.name.endswith(".h5ad"):
                    result.append(url._replace(path=fs.path.join(url.path, data_obj.name)))
    else:
        raise ScelVisException("Invalid URL scheme: %s" % url.scheme)
    return result


@cache.memoize()
def _load_data_cached(url, identifier):
    return data.load_data(url, identifier)


@cache.memoize()
def load_all_metadata():
    """Load all meta data information from ``settings.DATA_SOURCES``.

    A data source can either be a URL to a file ending on ``.h5ad`` or a directory that contains ``.h5ad`` files.
    """
    result = []

    if settings.FAKE_DATA:
        result = [data.fake_data().metadata]

    for url in settings.DATA_SOURCES:
        if url.path.endswith(".h5ad"):
            logger.info("Loading single dataset from data source %s", data.redacted_urlunparse(url))
            identifier = fs.path.basename(url.path)[: -len(".h5ad")]
            result.append(_load_data_cached(url, identifier).metadata)
        else:
            lst = []
            for match in glob_data_sets(url):
                identifier = fs.path.basename(match.path)[: -len(".h5ad")]
                lst.append(_load_data_cached(match, identifier).metadata)
            if lst:
                logger.info("Loaded %d data sets from data directory.", len(lst))
            else:
                logger.warn(
                    "No data sets found in data directory %s", data.redacted_urlunparse(url)
                )
            result += lst
    return result


@cache.memoize()
def load_data(identifier):
    """Load data for the given identifier from data or upload directory."""
    if identifier == data.FAKE_DATA_ID:
        return data.fake_data()
    else:
        urls = settings.DATA_SOURCES
        if settings.UPLOAD_DIR:
            urls += [urllib.parse.urlparse("file://%s" % settings.UPLOAD_DIR)]
        for url in urls:
            if fs.path.basename(url.path)[: -len(".h5ad")] == identifier:
                return data.load_data(url, identifier)
            else:
                if identifier and does_exist(url, identifier + ".h5ad"):
                    return _load_data_cached(
                        url._replace(path=fs.path.join(url.path, identifier + ".h5ad")), identifier
                    )


@cache.memoize()
def load_metadata(identifier):
    """Load metadata for the given identifier from data or upload directory."""
    data = load_data(identifier)
    if data:
        return data.metadata
