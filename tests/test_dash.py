"""Tests for the Dash application."""

import dash
import pytest
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


def test_render_upload(dash_duo, scelvis_settings):
    """test upload functionality"""

    app = dash.testing.application_runners.import_app("scelvis.app")
    dash_duo.start_server(app)

    # Click "Go To" menu.
    item_goto = dash_duo.wait_for_element_by_css_selector("#page-goto")
    item_goto.click()
    # Click contained "Upload" item.
    item_upload = dash_duo.wait_for_element_by_css_selector("#menu-item-upload")
    item_upload.click()
    # Make sure that the site navigates to the "upload data" page.
    dash_duo.wait_for_text_to_equal("#page-brand", "Upload Data")
    # TODO: actually test the upload page


def test_render_cell_annotation(dash_duo, scelvis_settings):
    """Click through the cell annotation"""

    app = dash.testing.application_runners.import_app("scelvis.app")
    dash_duo.start_server(app)

    # Click "Go To" menu.
    item = dash_duo.wait_for_element_by_css_selector("#page-goto")
    item.click()

    # Click contained "fake data" item.
    item = dash_duo.wait_for_element_by_css_selector("#menu-item-builtin-fake-data")
    item.click()

    # Make sure that the site navigates to the "fake data" page.
    dash_duo.wait_for_text_to_equal("#page-brand", "fake data")

    # Click through the different plot types.
    item = dash_duo.wait_for_element_by_css_selector("#meta_plot_type *:nth-child(1)")
    item.click()
    # plot = dash_duo.wait_for_element_by_css_selector('#meta_scatter_plot')
    # dash_duo.take_snapshot('meta_scatter')

    item = dash_duo.wait_for_element_by_css_selector("#meta_filter_cells_button")
    item.click()

    dropdown = dash_duo.wait_for_element_by_css_selector("#meta_filter_cells_attribute")
    dropdown.click()
    menu = dropdown.find_element_by_css_selector("div.Select-menu-outer")
    options = menu.find_elements_by_css_selector("div.VirtualizedSelectOption")
    options[1].click()

    checklist = dash_duo.wait_for_element_by_css_selector("#meta_filter_cells_choice")
    options = checklist.find_elements_by_css_selector("*")
    options[0].click()

    button = dash_duo.wait_for_element_by_css_selector("#meta_filter_cells_reset")
    button.click()

    item = dash_duo.wait_for_element_by_css_selector("#meta_plot_type *:nth-child(2)")
    item.click()
    # select sth from the dropdown
    dropdown = dash_duo.wait_for_element_by_css_selector("#meta_violin_select_vars")
    dropdown.click()
    menu = dropdown.find_element_by_css_selector("div.Select-menu-outer")
    options = menu.find_elements_by_css_selector("div.VirtualizedSelectOption")
    options[1].click()
    # plot = dash_duo.wait_for_element_by_css_selector('#meta_violin_plot')
    # dash_duo.take_snapshot('meta_violin')

    item = dash_duo.wait_for_element_by_css_selector("#meta_plot_type *:nth-child(3)")
    item.click()
    # plot = dash_duo.wait_for_element_by_css_selector('#meta_bar_plot')
    # dash_duo.take_snapshot('meta_bar')


def test_render_gene_annotation(dash_duo, scelvis_settings):
    """Click through the gene annotation"""

    app = dash.testing.application_runners.import_app("scelvis.app")
    dash_duo.start_server(app)

    # uncomment this to call up interactive ipython session
    # from IPython import embed
    # embed()

    # Click "Go To" menu.
    item_goto = dash_duo.wait_for_element_by_css_selector("#page-goto")
    item_goto.click()

    # Click contained "fake data" item.
    item = dash_duo.wait_for_element_by_css_selector("#menu-item-builtin-fake-data")
    item.click()

    # Make sure that the site navigates to the "fake data" page.
    dash_duo.wait_for_text_to_equal("#page-brand", "fake data")

    # Navigate to gene expression.
    tab = dash_duo.wait_for_element_by_css_selector(".nav-tabs li:nth-child(3)")
    tab.click()

    # Click through the different plot types.
    item = dash_duo.wait_for_element_by_css_selector("#expression_plot_type *:nth-child(1)")
    item.click()

    # select a gene from the dropdown
    dropdown = dash_duo.wait_for_element_by_css_selector("#expression_select_genes")
    dropdown.click()
    menu = dropdown.find_element_by_css_selector("div.Select-menu-outer")
    options = menu.find_elements_by_css_selector("div.VirtualizedSelectOption")
    options[1].click()
    # plot = dash_duo.wait_for_element_by_css_selector('#expression_scatter_plot')
    # dash_duo.take_snapshot('expression_scatter')

    item = dash_duo.wait_for_element_by_css_selector("#expression_plot_type *:nth-child(2)")
    item.click()
    # plot = dash_duo.wait_for_element_by_css_selector('#expression_violin_plot')
    # dash_duo.take_snapshot('expression_violin')

    item = dash_duo.wait_for_element_by_css_selector("#expression_plot_type *:nth-child(3)")
    item.click()
    # plot = dash_duo.wait_for_element_by_css_selector('#expression_dot_plot')
    # dash_duo.take_snapshot('expression_dot')
