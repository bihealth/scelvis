"""Definitions of Dash layout."""

import os.path

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

from .. import settings, store
from ..__init__ import __version__
from . import genes, cells


def render_not_found():
    """Return site content in the case that the dataset could not be found."""
    return html.Div(
        children=[
            html.H3("Dataset Not Found"),
            html.P("The dataset that you specified could not be found!"),
        ]
    )


def render_home():
    """Return site content for the home screen."""
    with open(os.path.join(os.path.dirname(__file__), "..", "static", "home.md")) as inputf:
        home_md = inputf.read()
        if settings.PUBLIC_URL_PREFIX:
            home_md = home_md.replace("](assets", "](%s/dash/assets" % settings.PUBLIC_URL_PREFIX)
    return dbc.Row(dbc.Col(html.Div(dcc.Markdown(home_md), className="home-text")))


def render_dataset(identifier):
    """Display the dataset."""
    data = store.load_data(identifier)
    if data:
        return render_dataset_page(data)
    else:
        return ""


def render_upload():
    """Display the upload form."""
    return dbc.Row(
        dbc.Col(
            [
                html.P(
                    children=[
                        (
                            "Use the following form to upload your .h5ad file.  You can convert "
                            "your pipeline output to such a file using "
                        ),
                        html.A(
                            children=[html.I(className="fas fa-redo mr-1"), "Convert Data"],
                            href="%s/convert/" % settings.PUBLIC_URL_PREFIX,
                        ),
                        ".",
                    ]
                ),
                dcc.Upload(
                    html.Div(
                        children=[
                            html.I(className="fas fa-file-upload mx-1"),
                            "Drag and Drop or ",
                            html.A(
                                children=[
                                    html.I(className="fas fa-hand-pointer mr-1"),
                                    "Click and Select files",
                                ]
                            ),
                            ".",
                        ],
                        className="text-center p-3 m-5 bg-light border border-primary rounded",
                    ),
                    id="file-upload",
                ),
            ]
        )
    )


def render_navbar():
    """Render the site navbar"""
    return dbc.Navbar(
        dbc.Container(
            children=[
                # Use row and col to control vertical alignment of logo / brand
                dbc.NavbarBrand(
                    [html.I(className="fas fa-home mr-1"), settings.HOME_BRAND], id="page-brand"
                ),
                dbc.Nav(
                    dbc.DropdownMenu(
                        children=render_children_goto(),
                        nav=True,
                        in_navbar=True,
                        label="Go To",
                        id="page-goto",
                    ),
                    navbar=True,
                ),
            ]
        ),
        # className="mb-5",
        color="primary",
        dark=True,
        id="page-navbar",
    )


def render_children_goto():
    """Render the "Go To" menu."""
    result = [
        dbc.DropdownMenuItem("Home", href="/dash/home"),
        dbc.DropdownMenuItem(divider=True),
        dbc.DropdownMenuItem("Data Sets", header=True),
    ]
    metas = store.load_all_metadata()
    for meta in metas:
        result.append(
            dbc.DropdownMenuItem(
                meta.short_title, id="menu-item-%s" % meta.id, href="/dash/viz/%s" % meta.id
            )
        )
    if not metas:
        result.append(
            dbc.DropdownMenuItem(html.Span("no dataset available", className="text-muted"))
        )
    if settings.UPLOAD_ENABLED or settings.CONVERSION_ENABLED:
        result.append(dbc.DropdownMenuItem(divider=True))
    if settings.UPLOAD_ENABLED:
        result.append(
            dbc.DropdownMenuItem(
                html.Span(
                    children=[html.I(className="fas fa-cloud-upload-alt mr-1"), "Upload Data"]
                ),
                id="menu-item-upload",
                href="/dash/upload",
            )
        )
    if settings.CONVERSION_ENABLED:
        result.append(
            dbc.DropdownMenuItem(
                html.Span(children=[html.I(className="fas fa-redo mr-1"), "Convert Data"]),
                id="menu-item-convert",
                href="%s/convert/" % settings.PUBLIC_URL_PREFIX,
                external_link=True,
            )
        )
    return result


def render_main_content():
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


def render_footer():
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


def render_dataset_page(data):
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
                html.Div(cells.render(data), className="mx-2 mt-2"),
                label="Cell Annotation",
                tab_id="tab-cells",
            ),
            dbc.Tab(
                html.Div(genes.render(data), className="mx-2 mt-2"),
                label="Gene Expression",
                tab_id="tab-genes",
            ),
        ],
        active_tab="tab-cells",
    )


def build_layout():
    """Build the overall Dash app layout"""
    return html.Div(
        children=[
            # Represents the URL bar, doesn't render anything.
            dcc.Location(id="url", refresh=False),
            # Navbar, content, footer.
            render_navbar(),
            render_main_content(),
            render_footer(),
        ],
        id="_dash-app-content",  # make pytest-dash happy
    )
