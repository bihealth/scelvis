"""Setup of the Flask cache."""

from flask_caching import Cache
from logzero import logger

from . import settings

#: The cache configuration, built from ``.settings``.
CACHE_CONFIG = {
    "DEBUG": settings.DEBUG,
    "CACHE_TYPE": settings.CACHE_TYPE,
    "CACHE_DEFAULT_TIMEOUT": settings.CACHE_DEFAULT_TIMEOUT,
    "CACHE_DIR": settings.CACHE_DIR,
    "CACHE_REDIS_URL": settings.CACHE_REDIS_URL,
}

#: The global cache instance.
cache = Cache()


def setup_cache(app):
    """Setup the Dash app's Flask app with the cache."""
    logger.info("Using cache configuration %s", CACHE_CONFIG)
    cache.init_app(app.server, config=CACHE_CONFIG)
