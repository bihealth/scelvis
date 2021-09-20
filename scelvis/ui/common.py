"""Rendering of common Dash components."""

import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html

import json
import numpy as np

from ..exceptions import ScelVisException


def render_plot(data_type, plot_type):
    """Render the plot component."""
    if plot_type not in ("scatter", "violin", "box", "bar", "dot"):
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


def render_filter_cells_collapse(data, token):

    children = [
        dbc.Button(
            "filter cells",
            id="%s_filter_cells_button" % token,
            className="text-left",
            color="primary",
            outline=True,
            size="md",
        ),
        dbc.Collapse(
            dbc.Card(dbc.CardBody(render_filter_cells_controls(data, token))),
            id="%s_filter_cells_collapse" % token,
        ),
    ]

    if token == "meta":
        # store list of choices in hidden div
        filters = {}
        children.append(
            html.Div(
                id="filter_cells_filters", style={"display": "none"}, children=json.dumps(filters)
            )
        )

    return html.Div(children=children, id="%s_filter_cells_div" % token)


def render_filter_cells_controls(data, token):

    options = (
        [{"label": "ncells", "value": "ncells"}]
        + [{"label": c, "value": c} for c in data.categorical_meta]
        + [{"label": c, "value": c} for c in data.numerical_meta]
        + [{"label": g, "value": g} for g in data.genes]
    )

    return [
        html.Div(
            [dcc.Dropdown(id="%s_filter_cells_attribute" % token, options=options, value="None")]
        ),
        # use Slider for ncells, CheckList for categorical filters and RangeSlider for numerical ones
        # and show the ones that's appropriate
        html.Div(
            id="%s_filter_cells_ncells_div" % token,
            children=[
                html.P(),
                dcc.Slider(
                    id=("%s_filter_cells_ncells" % token),
                    min=0,
                    max=1,
                    marks={0: "0", 1: "1"},
                    value=1,
                ),
                html.P(),
            ],
            style={"display": "none"},
        ),
        html.Div(
            id="%s_filter_cells_choice_div" % token,
            children=[
                dcc.Checklist(
                    id=("%s_filter_cells_choice" % token),
                    options=[{"label": "a", "value": "a"}],
                    value=["a"],
                    className="mt-2",
                    inputClassName="mr-1",
                )
            ],
            style={"display": "none"},
        ),
        html.Div(
            id="%s_filter_cells_range_div" % token,
            children=[
                html.P(),
                dcc.RangeSlider(
                    id=("%s_filter_cells_range" % token),
                    min=0,
                    max=1,
                    marks={0: "0", 1: "1"},
                    value=[0, 1],
                ),
                html.P(),
            ],
            style={"display": "none"},
        ),
        html.P(children="", id="%s_filter_cells_status" % token, style={"fontSize": 14}),
        dbc.Row(
            [
                dbc.Button(
                    "update plot",
                    id="%s_filter_cells_update" % token,
                    color="primary",
                    className="mr-1",
                    size="sm",
                ),
                dbc.Button(
                    "reset filters",
                    id="%s_filter_cells_reset" % token,
                    className="mr-1",
                    color="secondary",
                    size="sm",
                ),
            ]
        ),
    ]


def apply_filter_cells_filters(data, filters_json):

    ncells_tot = data.ad.obs.shape[0]
    take = np.ones(ncells_tot, dtype=bool)
    filters = json.loads(filters_json)
    for col, selected in filters.items():
        if selected is None:
            continue
        if col == "ncells":
            take = take & (np.random.random(ncells_tot) < selected / ncells_tot)
        elif col in data.categorical_meta:
            take = take & data.ad.obs_vector(col).isin(selected)
        else:
            take = (
                take
                & (data.ad.obs_vector(col) >= float(selected[0]))
                & (data.ad.obs_vector(col) <= float(selected[1]))
            )
    return data.ad[take, :]


def auto_tick(data_range, max_tick=10, tf_inside=False):
    """
    https://stackoverflow.com/questions/4947682/intelligently-calculating-chart-tick-positions
    tool function that automatically calculate optimal ticks based on range and the max number of ticks
    :param data_range:   range of data, e.g. [-0.1, 0.5]
    :param max_tick:     max number of ticks, an interger, default to 10
    :param tf_inside:    True/False if only allow ticks to be inside
    :return:             list of ticks
    """
    data_span = data_range[1] - data_range[0]
    scale = 10.0 ** np.floor(
        np.log10(data_span)
    )  # scale of data as the order of 10, e.g. 1, 10, 100, 0.1, 0.01, ...
    list_tick_size_nmlz = [
        5.0,
        2.0,
        1.0,
        0.5,
        0.2,
        0.1,
        0.05,
        0.02,
        0.01,
    ]  # possible tick sizes for normalized data in range [1, 10]
    tick_size_nmlz = 1.0  # initial tick size for normalized data
    for i in range(
        len(list_tick_size_nmlz)
    ):  # every loop reduces tick size thus increases tick number
        num_tick = (
            data_span / scale / list_tick_size_nmlz[i]
        )  # number of ticks for the current tick size
        if num_tick > max_tick:  # if too many ticks, break loop
            tick_size_nmlz = list_tick_size_nmlz[i - 1]
            break
    tick_size = tick_size_nmlz * scale  # tick sizse for the original data
    if tick_size > 0:
        ticks = (
            np.unique(np.arange(data_range[0] / tick_size, data_range[1] / tick_size).round())
            * tick_size
        )  # list of ticks
    else:
        ticks = data_range[0].round()

    if tf_inside:  # if only allow ticks within the given range
        ticks = ticks[(ticks >= data_range[0]) * (ticks <= data_range[1])]

    return ticks
