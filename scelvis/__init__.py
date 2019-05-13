from ._version import get_versions as _get_versions

#: The application's version from versioneer.
__version__ = _get_versions()["version"]
del _get_versions  # noqa
