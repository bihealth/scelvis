"""Access to the datasets"""

import glob
import os.path

from logzero import logger

from . import data, settings
from .cache import cache


#: The file name to use for identifying data directories.
ABOUT_FILENAME = "about.md"


@cache.memoize()
def load_all_metadata():
    """Load all meta data information ``data_dir``.

    If ``data_dir`` itself contains a file called ``about.md``, it is assumed that only one dataset is available.
    Otherwise, all sub directories will be scanned for ``about.md`` and one entry is returned for each dataset.
    """
    data_dir = settings.DATA_DIR
    if os.path.exists(os.path.join(data_dir, ABOUT_FILENAME)):
        logger.info("Loading single dataset from data directory %s", data_dir)
        return [data.load_metadata(data_dir)]
    else:
        result = []
        for path in glob.glob(os.path.join(data_dir, "*", ABOUT_FILENAME)):
            result.append(data.load_metadata(os.path.dirname(path)))
        if result:
            logger.info("Loaded %d data sets from data directory.", len(result))
        else:
            logger.warn("No data sets found in data directory %s", data_dir)
        return result


@cache.memoize()
def load_metadata(identifier):
    """Load metadata for the given identifier."""
    return data.load_metadata(os.path.join(settings.DATA_DIR, identifier))


@cache.memoize()
def load_data(identifier):
    """Load data for the given identifier."""
    return data.load_data(os.path.join(settings.DATA_DIR, identifier))
