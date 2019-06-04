"""Rendering of common Dash components."""

import dash_core_components as dcc
import dash_html_components as html

from ..exceptions import ScelVisException


def render_plot(data_type, plot_type):
    """Render the plot component."""
    if plot_type not in ("scatter", "violin", "bar", "dot"):
        raise ScelVisException("Invalid plot type: %s" % plot_type)
    elif data_type not in ("expression", "meta"):
        raise ScelVisException("Invalid data type: %s" % data_type)
    else:
        return (
            [
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
            ],
            plot_type == "scatter",
        )


def render_subsampling_dropdown(data, token):
    sample_choices = [{"label": g, "value": g} for g in [1000, 5000] if g < data.meta.shape[0]] + [
        {"label": "all", "value": "all"}
    ]
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
