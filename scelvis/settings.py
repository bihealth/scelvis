"""Configuration settings.

Configuration is kept as global state.  This is not how Python apps usually do it but so far we have not found a better
way to do it with Dash.
"""

# Currently constant settings.

#: String to use for the Bootstrap "brand" at home.
HOME_BRAND = "SCelVis Home"

#: Maximal upload data files.
MAX_UPLOAD_DATA_SIZE = 1_000_000_000
#: Maximal upload text files.
MAX_UPLOAD_TEXT_SIZE = 1_000_000

# Currently configurable settings.

#: The prefix that this app will be served with.  This has to be properly set into Flask and Dash.
PUBLIC_URL_PREFIX = ""

#: Whether or not to enable debugging
DEBUG = False

#: Whether or not to show the fake data.
FAKE_DATA = False

#: Paths/URLs with data sources.
DATA_SOURCES = []

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

#: Temporary directory, e.g., as first upload directory.
TEMP_DIR = None
#: Maximal file upload size.
MAX_UPLOAD_SIZE = 512_000_000

#: Whether or not to enable file upload for visualization, default is ``True```.
UPLOAD_ENABLED = False
#: A specific download directory, by default a temporary directory is created.
UPLOAD_DIR = None

#: Whether or not to enable file upload for conversion, default is ``True``.
CONVERSION_ENABLED = False
