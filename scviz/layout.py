"""Definitions of Dash layout."""

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

from .__init__ import __version__

#: String to use for the Bootstrap "brand" at home.
HOME_BRAND = 'SCViz Home'

#: Site navbar.
navbar = dbc.NavbarSimple(
    children=[
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem("Home", href="/home"),
                dbc.DropdownMenuItem(divider=True),
                dbc.DropdownMenuItem("Data Sets", header=True),
                dbc.DropdownMenuItem("Page 2", href="/viz/2"),
                dbc.DropdownMenuItem("Page 3", href="/viz/3"),
            ],
            nav=True,
            in_navbar=True,
            label="Go To",
        )
    ],
    brand=HOME_BRAND,
    color="primary",
    dark=True,
    id="page-navbar"
)

#: Main page content.
main_content = html.Div(
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


#: Site footer.
footer = html.Footer(
    html.Div(
        children=[
            html.Div(
                children=[
                    html.Div(
                        children=[
                            html.Span("SVViz v%s by BIH CUBI" % __version__, className="text-muted")
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
                                children=[html.I(className="fab fa-github mr-1"), "GitHub Project"],
                                href="https://github.com/bihealth/scviz",
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


#: Overall layout.
layout = html.Div(
    [
        # Represents the URL bar, doesn't render anything.
        dcc.Location(id="url", refresh=False),
        # Navbar, content, footer.
        navbar,
        main_content,
        footer,
    ]
)
