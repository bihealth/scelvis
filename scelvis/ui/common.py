"""Rendering of common Dash components."""

import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html

import json

from ..exceptions import ScelVisException
from .. import callbacks,settings


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
                    id="%s_filter_cells_options" % token,
                    options=[{"label": c, "value": k} for k,c in enumerate(data.categorical_meta)],
                    value="None"
                ),
            ],
            title=("choose attribute by which to filter cells"),
        ),
    ]
    choices = {}
    for n,attribute in enumerate(data.categorical_meta):
        output.append(
            dcc.Checklist(
                id = "%s_filter_cells_%i" % (token, n),
                options = [
                    {"label": c, "value": k} for k,c in enumerate(data.ad.obs[attribute].cat.categories)
                ],
                value=list(range(len(data.ad.obs[attribute].cat.categories))),
                className="mt-2",
                inputClassName="mr-1",
                style = {"display": "none"},
                )
            )
        choices[attribute] = list(data.ad.obs[attribute].cat.categories)

    output.append(html.Div(id="%s_filter_cells_choices" % token,
                           style={"display": "block"},
                           children=json.dumps(choices)))

    return output

