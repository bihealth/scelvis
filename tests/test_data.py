"""Tests for the data loading module."""

from urllib import parse
import os.path

import pytest

from scelvis import settings
from scelvis import data


@pytest.yield_fixture(scope="function")
def scelvis_settings(tmpdir):
    old_tmpdir = settings.TEMP_DIR
    old_cachedir = settings.CACHE_DIR
    old_cachetype = settings.CACHE_TYPE
    old_cachetimeout = settings.CACHE_DEFAULT_TIMEOUT
    old_uploadenabled = settings.UPLOAD_ENABLED
    old_conversionenabled = settings.CONVERSION_ENABLED
    old_fakedata = settings.FAKE_DATA

    settings.TEMP_DIR = str(tmpdir)
    settings.CACHE_DIR = str(tmpdir)
    settings.CACHE_TYPE = "simple"
    settings.CACHE_DEFAULT_TIMEOUT = 60
    settings.UPLOAD_ENABLED = True
    settings.CONVERSION_ENABLED = True
    settings.FAKE_DATA = True

    yield

    settings.TEMP_DIR = old_tmpdir
    settings.CACHE_DIR = old_cachedir
    settings.CACHE_TYPE = old_cachetype
    settings.CACHE_DEFAULT_TIMEOUT = old_cachetimeout
    settings.UPLOAD_ENABLED = old_uploadenabled
    settings.CONVERSION_ENABLED = old_conversionenabled
    settings.FAKE_DATA = old_fakedata


def test_load_data_from_file(scelvis_settings):
    url = "file://%s" % os.path.join(os.path.dirname(__file__), "..", "examples", "dummy.h5ad")
    data_source = parse.urlparse(url)
    result = data.load_data(data_source, "dummy")
    assert result.metadata
