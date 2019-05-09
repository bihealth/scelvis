"""Implementation of the app callbacks."""

import os.path

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from logzero import logger

from .layout import HOME_BRAND, render_dataset
from .ui import cells, PLOT_HEIGHT
from . import store


def get_page(pathname):
    """Helper function for parsing the URL path."""
    pathname = pathname or "/"
    if pathname in ("/", "/home"):
        return {"page": "home"}
    else:
        tokens = pathname.split("/")
        if len(tokens) < 3:
            return {"page": None}
        else:
            return {"page": "viz", "dataset": tokens[2]}


def register_page_content(app):
    """Register the display of the page content with the app."""

    @app.callback(
        dash.dependencies.Output("page-content", "children"),
        [dash.dependencies.Input("url", "pathname")],
    )
    def render_page_content(pathname):
        page = get_page(pathname)
        if page["page"] == "home":
            return display_home()
        elif page["page"] == "viz":
            return display_dataset(page["dataset"])
        else:
            return display_not_found()


def register_page_brand(app):
    """Register the display of the page brand with the app."""

    @app.callback(
        dash.dependencies.Output("page-navbar", "brand"),
        [dash.dependencies.Input("url", "pathname")],
    )
    def render_page_brand(pathname):
        page = get_page(pathname)
        if page["page"] in ("/", "/home"):
            return HOME_BRAND
        elif page["page"] == "viz":
            try:
                return "Dataset: %s" % store.load_metadata(page["dataset"]).title
            except FileNotFoundError:
                return "Dataset Not Found"
        else:
            return "Dataset Not Found"


def register_select_cell_plot_type(app):
    """Register callback for changing the controls on updating "CellAnnotation" plot type."""

    @app.callback(
        [dash.dependencies.Output("meta_plot_controls", "children")],
        [
            dash.dependencies.Input("url", "pathname"),
            dash.dependencies.Input("meta_plot_type", "value"),
        ],
    )
    def update_meta_plot_controls(pathname, plot_type):
        data = store.load_data(get_page(pathname)["dataset"])
        plots = {
            "scatter": cells.render_controls_scatter,
            "violin": cells.render_controls_violin,
            "bar": cells.render_controls_bars,
        }
        return plots[plot_type](data)

    @app.callback(
        [
            dash.dependencies.Output("meta_plot", "children"),
            dash.dependencies.Output("meta_select_cell_sample", "disabled"),
        ],
        [dash.dependencies.Input("meta_plot_type", "value")],
    )
    def update_meta_plots(plot_type):
        return cells.render_plot(plot_type)


def register_update_cell_scatter_plot_params(app):
    """Register handlers on updating scatter plot controls."""

    @app.callback(
        [
            dash.dependencies.Output("meta_scatter_plot", "figure"),
            dash.dependencies.Output("meta_scatter_download", "href"),
            dash.dependencies.Output("meta_scatter_download", "hidden"),
        ],
        [
            dash.dependencies.Input("url", "pathname"),
            dash.dependencies.Input("meta_scatter_select_x", "value"),
            dash.dependencies.Input("meta_scatter_select_y", "value"),
            dash.dependencies.Input("meta_scatter_select_color", "value"),
            dash.dependencies.Input("meta_select_cell_sample", "value"),
        ],
    )
    def get_meta_scatterplot(pathname, xc, yc, col, sample_size):
        data = store.load_data(get_page(pathname)["dataset"])
        return cells.render_plot_scatter(data, xc, yc, col, sample_size)


def register_update_cell_violin_plot_params(app):
    """Register handlers on updating violin plot controls."""

    @app.callback(
        [
            dash.dependencies.Output("meta_violin_plot", "figure"),
            dash.dependencies.Output("meta_violin_download", "href"),
            dash.dependencies.Output("meta_violin_download", "hidden"),
        ],
        [
            dash.dependencies.Input("url", "pathname"),
            dash.dependencies.Input("meta_violin_select_vars", "value"),
            dash.dependencies.Input("meta_violin_select_group", "value"),
            dash.dependencies.Input("meta_violin_select_split", "value"),
            dash.dependencies.Input("meta_select_cell_sample", "value"),
        ],
    )
    def get_meta_violinplot(pathname, variables, group, split, sample_size):
        data = store.load_data(get_page(pathname)["dataset"])
        return cells.render_plot_violin(data, variables, group, split, sample_size)


def register_update_cell_bar_chart_params(app):
    """Register handlers on updating bar chart controls."""

    @app.callback(
        dash.dependencies.Output("meta_bar_options", "options"),
        [dash.dependencies.Input("meta_bar_select_split", "value")],
    )
    def toggle_meta_bar_options(split):
        if split is None:
            return [dict(label="normalized", value="normalized")]
        else:
            return [dict(label=c, value=c) for c in ["normalized", "stacked"]]


def register_select_cell_plot_type_update_plots(app):
    """Register callback for changing the controls on updating "CellAnnotation" plot type."""


def display_home():
    """Return site content for the home screen."""
    with open(os.path.join(os.path.dirname(__file__), "static", "home.md")) as inputf:
        home_md = inputf.read()
    return dbc.Row(dbc.Col(html.Div(dcc.Markdown(home_md))))


def display_not_found():
    """Return site content in the case that the dataset could not be found."""
    return html.Div(
        children=[
            html.H3("Dataset Not Found"),
            html.P("The dataset that you specified could not be found!"),
        ]
    )


def display_dataset(identifier):
    """Display the dataset."""
    try:
        data = store.load_data(identifier)
    except FileNotFoundError:
        return ""
    else:
        return render_dataset(data)
