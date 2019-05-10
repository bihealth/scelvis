"""Definitions of Dash layout."""

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

from .__init__ import __version__
from .store import load_all_metadata
from .ui import render_tab_cells, render_tab_genes

#: String to use for the Bootstrap "brand" at home.
HOME_BRAND = "SCelVis Home"


def navbar():
    """Render the site navbar"""
    return dbc.NavbarSimple(
        children=[
            dbc.DropdownMenu(
                children=children_goto(), nav=True, in_navbar=True, label="Go To", id="page-goto"
            )
        ],
        brand=HOME_BRAND,
        color="primary",
        dark=True,
        id="page-navbar",
    )


def children_goto():
    """Render the "Go To" menu."""
    result = [
        dbc.DropdownMenuItem("Home", href="/home"),
        dbc.DropdownMenuItem(divider=True),
        dbc.DropdownMenuItem("Data Sets", header=True),
    ]
    metas = load_all_metadata()
    for meta in metas:
        result.append(dbc.DropdownMenuItem(meta.short_title, href="/viz/%s" % meta.id))
    if not metas:
        result.append(
            dbc.DropdownMenuItem(html.Span("no dataset available", className="text-muted"))
        )
    return result


def main_content():
    """Render page main content"""
    return html.Div(
        children=[
            dbc.Row(
                dbc.Col(
                    # content will be rendered in this element
                    html.Div(id="page-content")
                )
            )
        ],
        className="container pt-3",
    )


def footer():
    """Render page footer"""
    return html.Footer(
        html.Div(
            children=[
                html.Div(
                    children=[
                        html.Div(
                            children=[
                                html.Span(
                                    "SVViz v%s by BIH CUBI" % __version__, className="text-muted"
                                )
                            ],
                            className="col-6",
                        ),
                        html.Div(
                            children=[
                                html.A(
                                    children=[
                                        html.I(className="fas fa-globe-europe mr-1"),
                                        "CUBI Homepage",
                                    ],
                                    href="https://www.cubi.bihealth.org",
                                    className="text-muted mr-3",
                                ),
                                html.A(
                                    children=[
                                        html.I(className="fab fa-github mr-1"),
                                        "GitHub Project",
                                    ],
                                    href="https://github.com/bihealth/scelvis",
                                    className="text-muted",
                                ),
                            ],
                            className="col-6 text-right",
                        ),
                    ],
                    className="row",
                )
            ],
            className="container",
        ),
        className="footer",
    )


def render_dataset(data):
    """Render the page main content for dataset visualization."""
    return dbc.Tabs(
        children=[
            dbc.Tab(
                html.Div(
                    dbc.Row(dbc.Col(dcc.Markdown(data.metadata.readme))), className="mx-2 mt-2"
                ),
                label="About",
                tab_id="tab-about",
            ),
            dbc.Tab(
                html.Div(render_tab_cells(data), className="mx-2 mt-2"),
                label="Cell Annotation",
                tab_id="tab-cells",
            ),
            dbc.Tab(
                html.Div(render_tab_genes(data), className="mx-2 mt-2"),
                label="Gene Expression",
                tab_id="tab-genes",
            ),
        ],
        active_tab="tab-cells",
    )


def layout():
    """Overall layout"""
    return html.Div(
        children=[
            # Represents the URL bar, doesn't render anything.
            dcc.Location(id="url", refresh=False),
            # Navbar, content, footer.
            navbar(),
            main_content(),
            footer(),
        ]
    )
