"""Access to all data outside of the session.

This module uses the functions from ``data`` to load the data from the appropriate location
and memoizes the loaded data in the Flask cacche.  This module is aware of the Flask cache and
thus also of the Dash app.
"""

import os.path
import urllib

import fs
from logzero import logger

from . import data, settings
from .cache import cache


def get_file_system(url):
    """Return file system from the given ``url``."""
    url = urllib.parse.urlparse(url)
    if url.scheme and url.scheme != "file":
        base_url = "%s://%s" % (url.scheme, url.netloc)
    else:
        base_url = url.netloc
    return url, fs.open_fs(base_url)


@cache.memoize()
def load_all_metadata():
    """Load all meta data information ``data_dir``.

    If ``data_dir`` itself contains a file called ``about.md``, it is assumed that only one dataset is available.
    Otherwise, all sub directories will be scanned for ``about.md`` and one entry is returned for each dataset.
    """
    result = []
    for data_dir in settings.DATA_DIRS:
        url, the_fs = get_file_system(data_dir)
        if the_fs.exists(fs.path.join(url.path, settings.ABOUT_FILENAME)):
            logger.info("Loading single dataset from data directory %s", data_dir)
            result.append(data.load_metadata(data_dir))
        else:
            for match in the_fs.glob(fs.path.join(url.path, "*", settings.ABOUT_FILENAME)):
                logger.info("Found data set %s", fs.path.dirname(match.path))
                result.append(data.load_metadata(the_fs, fs.path.dirname(match.path)))
            if result:
                logger.info("Loaded %d data sets from data directory.", len(result))
            else:
                logger.warn("No data sets found in data directory %s", data_dir)
    return result


@cache.memoize()
def load_metadata(identifier):
    """Load metadata for the given identifier from data or upload directory."""
    for base_dir in settings.DATA_DIRS + [settings.UPLOAD_DIR]:
        if not base_dir:
            continue
        url, the_fs = get_file_system(base_dir)
        full_path = fs.path.join(url.path, identifier)
        if the_fs.exists(full_path):
            return data.load_metadata(the_fs, full_path)


@cache.memoize()
def load_data(identifier):
    """Load data for the given identifier from data or upload directory."""
    for base_dir in settings.DATA_DIRS + [settings.UPLOAD_DIR]:
        if not base_dir:
            continue
        url, the_fs = get_file_system(base_dir)
        full_path = fs.path.join(url.path, identifier)
        if the_fs.exists(full_path):
            return data.load_data(the_fs, full_path)
