"""Tests for the Dash application."""

import pytest
from pytest_dash import wait_for
from pytest_dash.application_runners import import_app

from scelvis import settings


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


def test_render_upload(dash_threaded, scelvis_settings):
    driver = dash_threaded.driver
    dash_threaded(import_app("scelvis.app"))

    # Click "Go To" menu.
    item_goto = wait_for.wait_for_element_by_css_selector(driver, "#page-goto")
    item_goto.click()
    # Click contained "Upload" item.
    item_upload = wait_for.wait_for_element_by_css_selector(driver, "#menu-item-upload")
    item_upload.click()
    # Make sure that the site navigates to the "upload data" page.
    wait_for.wait_for_text_to_equal(driver, "#page-brand", "Upload Data")
    # TODO: actually test the upload page


def test_nagivate_to_data(dash_threaded, scelvis_settings):
    driver = dash_threaded.driver
    dash_threaded(import_app("scelvis.app"))

    # Click "Go To" menu.
    item_goto = wait_for.wait_for_element_by_css_selector(driver, "#page-goto")
    item_goto.click()
    # Click contained "fake data" item.
    item_upload = wait_for.wait_for_element_by_css_selector(driver, "#menu-item-builtin-fake-data")
    item_upload.click()
    # Make sure that the site navigates to the "fake data" page.
    wait_for.wait_for_text_to_equal(driver, "#page-brand", "fake data")


def test_render_cell_annotation(dash_threaded, scelvis_settings):
    """Click through the cell annotation"""
    driver = dash_threaded.driver
    dash_threaded(import_app("scelvis.app"))

    # Click "Go To" menu.
    item_goto = wait_for.wait_for_element_by_css_selector(driver, "#page-goto")
    item_goto.click()

    # Click contained "fake data" item.
    item_upload = wait_for.wait_for_element_by_css_selector(driver, "#menu-item-builtin-fake-data")
    item_upload.click()

    # Make sure that the site navigates to the "fake data" page.
    wait_for.wait_for_text_to_equal(driver, "#page-brand", "fake data")

    # Click through the different plot types.
    item = wait_for.wait_for_element_by_xpath(driver, "//label[contains(text(), 'scatter plot')]")
    item.click()
    item = wait_for.wait_for_element_by_xpath(driver, "//label[contains(text(), 'violin plot')]")
    item.click()
    item = wait_for.wait_for_element_by_xpath(driver, "//label[contains(text(), 'bar plot')]")
    item.click()


def test_render_gene_annotation(dash_threaded, scelvis_settings):
    """Click through the cell annotation"""
    driver = dash_threaded.driver
    dash_threaded(import_app("scelvis.app"))

    # Click "Go To" menu.
    item_goto = wait_for.wait_for_element_by_css_selector(driver, "#page-goto")
    item_goto.click()

    # Click contained "fake data" item.
    item_upload = wait_for.wait_for_element_by_css_selector(driver, "#menu-item-builtin-fake-data")
    item_upload.click()

    # Make sure that the site navigates to the "fake data" page.
    wait_for.wait_for_text_to_equal(driver, "#page-brand", "fake data")

    # Navigate to gene expression.
    tab = wait_for.wait_for_element_by_xpath(driver, "//a[contains(text(), 'Gene Expression')]")
    tab.click()

    # Click through the different plot types.
    item = wait_for.wait_for_element_by_xpath(
        driver, "(//label[contains(text(), 'scatter plot')])[2]"
    )
    item.click()
    item = wait_for.wait_for_element_by_xpath(
        driver, "(//label[contains(text(), 'violin plot')])[2]"
    )
    item.click()
    item = wait_for.wait_for_element_by_xpath(driver, "//label[contains(text(), 'dot plot')]")
    item.click()
