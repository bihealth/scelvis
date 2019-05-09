"""Implementation of the app callbacks."""

import os.path

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

from .layout import HOME_BRAND


def register_page_content(app):
    """Register the display of the page content with the app."""

    @app.callback(
        dash.dependencies.Output("page-content", "children"),
        [dash.dependencies.Input("url", "pathname")],
    )
    def render_page_content(pathname):
        pathname = pathname or "/"
        if pathname in ("/", "/home"):
            return display_home()
        else:
            tokens = pathname.split("/")
            if len(tokens) < 3:
                return display_not_found()
            else:
                return display_dataset(tokens[2])


def register_page_brand(app):
    """Register the display of the page brand with the app."""

    @app.callback(
        dash.dependencies.Output("page-navbar", "brand"),
        [dash.dependencies.Input("url", "pathname")],
    )
    def render_page_brand(pathname):
        pathname = pathname or "/"
        if pathname in ("/", "/home"):
            return HOME_BRAND
        else:
            tokens = pathname.split("/")
            if len(tokens) < 3:
                return "Dataset: Not Found"
            else:
                return "Dataset: %s" % tokens[2]


def display_home():
    """Return site content for the home screen."""
    with open(os.path.join(os.path.dirname(__file__), "static", "home.md")) as inputf:
        home_md = inputf.read()
    return dbc.Row(
        dbc.Col(
            html.Div(
                dcc.Markdown(home_md),
            )
        )
    )

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
    return html.Div(children=[html.H3("Dataset: %s" % identifier)])



