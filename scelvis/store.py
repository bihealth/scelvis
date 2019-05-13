"""Access to all data outside of the session.

This module uses the functions from ``data`` to load the data from the appropriate location
and memoizes the loaded data in the Flask cacche.  This module is aware of the Flask cache and
thus also of the Dash app.
"""

import glob
import os.path

from logzero import logger

from . import data, settings
from .cache import cache


@cache.memoize()
def load_all_metadata():
    """Load all meta data information ``data_dir``.

    If ``data_dir`` itself contains a file called ``about.md``, it is assumed that only one dataset is available.
    Otherwise, all sub directories will be scanned for ``about.md`` and one entry is returned for each dataset.
    """
    data_dir = settings.DATA_DIR
    if os.path.exists(os.path.join(data_dir, settings.ABOUT_FILENAME)):
        logger.info("Loading single dataset from data directory %s", data_dir)
        return [data.load_metadata(data_dir)]
    else:
        result = []
        for path in glob.glob(os.path.join(data_dir, "*", settings.ABOUT_FILENAME)):
            result.append(data.load_metadata(os.path.dirname(path)))
        if result:
            logger.info("Loaded %d data sets from data directory.", len(result))
        else:
            logger.warn("No data sets found in data directory %s", data_dir)
        return result


@cache.memoize()
def load_metadata(identifier):
    """Load metadata for the given identifier from data or upload directory."""
    for base_dir in (settings.DATA_DIR, settings.UPLOAD_DIR):
        if not base_dir:
            continue
        full_path = os.path.join(base_dir, identifier)
        if os.path.exists(full_path):
            return data.load_metadata(full_path)


@cache.memoize()
def load_data(identifier):
    """Load data for the given identifier from data or upload directory."""
    for base_dir in (settings.DATA_DIR, settings.UPLOAD_DIR):
        if not base_dir:
            continue
        full_path = os.path.join(base_dir, identifier)
        if os.path.exists(full_path):
            return data.load_data(full_path)
