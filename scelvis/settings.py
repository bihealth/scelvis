"""Configuration settings.

Configuration is kept as global state.  This is not how Python apps usually do it but so far we have not found a better
way to do it with Dash.
"""

#: Whether or not to enable debugging
DEBUG = False

#: Path to the data directory.
DATA_DIR = None

#: The type of the cache to use from {"filesystem", "redis"}
CACHE_TYPE = "filesystem"
#: Default cache timeout.
CACHE_DEFAULT_TIMEOUT = None
#: For "filesystem" cache: directory to store the cache in.
CACHE_DIR = None
#: For "redis" cache: the URL to use for connecting to the cache.
CACHE_REDIS_URL = None

#: The height of the plot.
PLOT_HEIGHT = 500
