"""Rendering of common Dash components."""

import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html

from ..exceptions import ScelVisException


def render_plot(data_type, plot_type):
    """Render the plot component."""
    if plot_type not in ("scatter", "violin", "bar", "dot"):
        raise ScelVisException("Invalid plot type: %s" % plot_type)
    elif data_type not in ("expression", "meta"):
        raise ScelVisException("Invalid data type: %s" % data_type)
    else:
        return [
                dcc.Graph(id="%s_%s_plot" % (data_type, plot_type)),
                html.A(
                    children=[
                        html.I(className="fas fa-cloud-download-alt pr-1"),
                        "download data for this plot",
                    ],
                    id="%s_%s_download" % (data_type, plot_type),
                    download="plot_data.csv",
                    href="",
                    hidden=True,
                    target="_blank",
                ),
        ]


def render_subsampling_dropdown(data, token):

    sample_choices = [
        {"label": g, "value": g} for g in [1000, 5000] if g < data.ad.obs.shape[0]
    ] + [{"label": "all", "value": "all"}]
    return html.Div(
        [
            html.Label("select cell sample"),
            dcc.Dropdown(
                id="%s_select_cell_sample" % token,
                options=sample_choices,
                value="all",
                disabled=False,
            ),
        ],
        title="Use random sample of cells for faster rendering.",
    )

def render_filter_cells_collapse (data, token):

    return html.Div(
        [
            dbc.Button(
                "filter cells",
                id="%s_filter_cells_button" % token,
                className='mb-3',
                color='secondary',
            ),
            dbc.Collapse(
                dbc.Card(dbc.CardBody(render_filter_cells_controls(data,token))),
                id="%s_filter_cells_collapse" % token,
            ),
        ]
    )

def render_filter_cells_controls (data, token):

    output = [
        html.Div(
            [
                html.Label("filter cells by "),
                dcc.Dropdown(
                    id="%s_filter_cells_attribute" % token,
                    options=[{"label": c, "value": c} for c in data.categorical_meta],
                    value="None",
                    multi=True,
                ),
            ],
            title=("choose attribute by which to filter cells"),
        ),
    ]
    for attribute in data.categorical_meta:
        output.append(
            dcc.Checklist(
                id = "%s_filter_cells_%s" % (token, attribute),
                options = [
                    {"label": c, "value": c} for c in data.ad.obs[attribute].cat.categories
                    ],
                value=data.ad.obs[attribute].cat.categories,
                className="mt-2",
                inputClassName="mr-1",
                style = {"display": "none"},
                )
            )

    return output

def render_filter_cells_choices (data, token):

    return [
        html.Div(
            [
                html.Label("filter cells by "),
                dcc.Dropdown(
                    id="%s_filter_cells_attribute" % token,
                    options=[{"label": c, "value": c} for c in data.ad.obs_keys()],
                    value="None",
                ),
            ],
            title=("choose attribute by which to filter cells"),
        ),
        # Placeholder for the attribute-specific controls.
        dcc.Loading(id="filter_select_cell_controls", type="circle")
    ]

