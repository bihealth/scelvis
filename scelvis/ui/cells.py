"""Render the plots and controls for the gene tab pane."""

import urllib.parse

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import plotly.subplots as subplots

import numpy as np
import pandas as pd
import scanpy as sc
import json


from .. import settings
from . import colors, common


def render_controls_scatter(data):
    """Render "top left" controls for the scatter plot for the given ``data``."""
    return [
        html.Div(
            children=[
                html.Label("select x axis"),
                dcc.Dropdown(
                    id="meta_scatter_select_x",
                    options=[{"label": c, "value": c} for c in data.coords + data.numerical_meta],
                    value=data.coords[0] if len(data.coords) > 1 else data.numerical_meta[0],
                ),
                html.Label("select y axis"),
                dcc.Dropdown(
                    id="meta_scatter_select_y",
                    options=[{"label": c, "value": c} for c in data.coords + data.numerical_meta],
                    value=data.coords[1] if len(data.coords) > 1 else data.numerical_meta[1],
                ),
                html.Label("select coloring"),
                dcc.Dropdown(
                    id="meta_scatter_select_color",
                    options=[
                        {"label": c, "value": c}
                        for c in data.categorical_meta + data.numerical_meta
                    ],
                    value=data.categorical_meta[0],
                ),
                render_select_cells(data),
            ],
            title="Select x- and y-coordinates for embedding (TSNE or UMAP); "
            "color points according to cell annotation (e.g., cluster identity or n_genes); "
            "select groups of cells for differential expression",
        )
    ]


def render_controls_violin(data):
    """Render "top left" controls for the violin plot for the given ``data``."""
    return [
        html.Div(
            children=[
                html.Label("select variable(s) and scaling"),
                dcc.Dropdown(
                    id="meta_violin_select_vars",
                    options=[{"label": c, "value": c} for c in data.numerical_meta],
                    value=None,
                    multi=True,
                ),
                html.Label("select grouping"),
                dcc.Dropdown(
                    id="meta_violin_select_group",
                    options=[{"label": c, "value": c} for c in data.categorical_meta],
                    value=data.categorical_meta[0],
                ),
                html.Label("select split"),
                dcc.Dropdown(
                    id="meta_violin_select_split",
                    options=[{"label": c, "value": c} for c in data.categorical_meta],
                    value=None,
                ),
            ],
            title=(
                "Select one or more numerical variables (e.g., n_genes or n_counts) to display in violin plots; "
                "choose how to group the cells (e.g., by cluster); and optionally also split these groups by another "
                "variable (e.g., by genotype)."
            ),
        )
    ]


def render_controls_bars(data):
    """Render "top left" controls for the bar chart for the given ``data``."""
    return [
        html.Div(
            children=[
                html.Label("select grouping"),
                dcc.Dropdown(
                    id="meta_bar_select_group",
                    options=[{"label": c, "value": c} for c in data.categorical_meta],
                    value=data.categorical_meta[0],
                ),
                html.Label("select split"),
                dcc.Dropdown(
                    id="meta_bar_select_split",
                    options=[{"label": c, "value": c} for c in data.categorical_meta],
                    value=None,
                ),
                html.Label("other properties"),
                dcc.Checklist(
                    id="meta_bar_options",
                    options=[
                        {"label": "normalized", "value": "normalized"},
                        {"label": "stacked", "value": "stacked"},
                    ],
                    value=[],
                ),
            ],
            title=(
                "Select group for which cell numbers are tallied; "
                "these groups can optionally be split by another variable for sub-tallies; "
                "bars can be normalized and/or stacked."
            ),
        )
    ]


def render_select_cells(data):
    """Render the collapse for selecting cells for differential expression"""
    return html.Div(
        [
            html.P(),
            dbc.Button(
                "differential expression",
                id="select_cells_button",
                className="text-left",
                color="primary",
                outline=True,
                size="md",
            ),
            dbc.Collapse(
                dbc.Card(dbc.CardBody(render_select_cells_controls(data))),
                id="select_cells_collapse",
            ),
        ],
        title='define two groups of cells in the plot with "Box select" or "Lasso select" and hit "Run" to perform differential expression',
    )


def render_select_cells_controls(data):
    """render the controls for selecting cells for differential expression"""
    return html.Div(
        [
            dbc.Row(
                [
                    dbc.Button(
                        "group A",
                        id="select_cells_group_A",
                        className="mr-1",
                        color="secondary",
                        size="sm",
                    ),
                    dbc.Button(
                        "group B",
                        id="select_cells_group_B",
                        className="mr-1",
                        color="secondary",
                        size="sm",
                    ),
                    dbc.Button(
                        "reset",
                        id="select_cells_reset",
                        color="secondary",
                        className="mr-1",
                        size="sm",
                    ),
                    dbc.Button(
                        "run", id="select_cells_run", color="primary", className="mr-1", size="sm"
                    ),
                ]
            ),
            html.P(),
            dbc.Row(
                [
                    dbc.Button(
                        "view table",
                        id="select_cells_view",
                        color="link",
                        style={
                            "padding-left": 0,
                            "padding-right": 0,
                            "padding-top": 0,
                            "padding-bottom": 2,
                        },
                    ),
                    " or get ",
                    html.A(
                        children=[html.I(className="fas fa-cloud-download-alt pr-1"), "results"],
                        download="results.csv",
                        id="select_cells_results_download",
                        href="",
                        target="_blank",
                    ),
                    " or ",
                    html.A(
                        children=[html.I(className="fas fa-cloud-download-alt pr-1"), "parameters"],
                        download="parameters.csv",
                        id="select_cells_parameters_download",
                        href="",
                        target="_blank",
                    ),
                ],
                id="select_cells_get_results",
                style={"display": "none"},
            ),
            html.Div(id="select_cells_choices", style={"display": "none"}),
            html.Div(id="select_cells_results", style={"display": "none"}),
        ]
    )


def render_controls(data):
    """Render the (left) controls column for the given ``data``."""
    return [
        # Select plot type.
        html.Div(
            children=[
                html.Label("select plot type"),
                dcc.RadioItems(
                    id="meta_plot_type",
                    options=[
                        {"label": "scatter plot", "value": "scatter"},
                        {"label": "violin plot", "value": "violin"},
                        {"label": "bar plot", "value": "bar"},
                    ],
                    value="scatter",
                    labelStyle={"display": "block"},
                    inputClassName="mr-1",
                ),
            ],
            title=(
                "SCATTER: plot annotation on two-dimensional embedding, e.g., TSNE or UMAP; "
                "VIOLIN: plot distributions of numerical variables (n_cells, n_counts, ...) in cell groups; "
                "BAR: plot cell numbers/fractions per (sub-)group"
            ),
        ),
        html.Hr(),
        # Control for filtering of cells.
        common.render_filter_cells_collapse(data, "meta"),
        html.Hr(),
        # Placeholder for the plot-specific controls.
        html.Div(id="meta_plot_controls"),
    ]


def render(data):
    """Render the "Cell Annotation" content."""
    return dbc.Row(
        children=[
            dbc.Col(children=render_controls(data), className="col-3"),
            # Placeholder for the plot.
            dbc.Col(children=[dcc.Loading(id="meta_plot", type="circle")], className="col-9"),
        ]
    )


def render_plot_scatter(data, xc, yc, col, filters_json):
    """Render the scatter plot figure."""

    if xc is None or yc is None or col is None:
        return {}, "", True

    ad_here = common.apply_filter_cells_filters(data, filters_json)

    if col in data.categorical_meta:
        # select color palette
        colvals = ad_here.obs[col].unique()
        cm = colors.get_cm(colvals)

        # plot scatter for each category separately
        traces = []
        for n, cv in enumerate(colvals):
            traces.append(
                go.Scattergl(
                    x=ad_here.obs[ad_here.obs[col] == cv][xc],
                    y=ad_here.obs[ad_here.obs[col] == cv][yc],
                    text=ad_here.obs[ad_here.obs[col] == cv].index,
                    mode="markers",
                    opacity=0.7,
                    marker={
                        "size": 5,
                        "color": cm[n % 40],
                        "line": {"width": 0.1, "color": "gray"},
                    },
                    showlegend=True,
                    name=cv,
                )
            )
    else:
        # for numerical data, plot scatter all at once
        traces = [
            go.Scattergl(
                x=ad_here.obs[xc],
                y=ad_here.obs[yc],
                text=ad_here.obs.index,
                mode="markers",
                opacity=0.7,
                marker={
                    "size": 5,
                    "color": ad_here.obs[col],
                    "colorscale": "Viridis",
                    "colorbar": {"title": col, "titleside": "right"},
                    "showscale": True,
                },
                name=col,
            )
        ]

    plot_data = ad_here.obs[[xc, yc, col]]
    csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(
        plot_data.to_csv(index=True, header=True, encoding="utf-8")
    )

    fig = {
        "data": traces,
        "layout": go.Layout(
            xaxis={"title": xc},
            yaxis={"title": yc},
            margin={"l": 40, "b": 40, "t": 10, "r": 10},
            legend={"x": 1.05, "y": 1},
            plot_bgcolor="rgb(255,255,255)",
            hovermode="closest",
            height=settings.PLOT_HEIGHT,
        ),
    }

    return fig, csv_string, False


def render_plot_violin(data, variables, group, split, filters_json):
    """Render the violin plot figure."""

    if variables is None or len(variables) == 0 or group is None:
        return {}, "", True

    ad_here = common.apply_filter_cells_filters(data, filters_json)

    # select color palette
    if split is None:
        groupvals = ad_here.obs[group].unique()
        cm = colors.get_cm(groupvals)
    else:
        splitvals = ad_here.obs[split].unique()
        cm = colors.get_cm(splitvals)

    nvar = len(variables)

    fig = subplots.make_subplots(
        rows=nvar,
        cols=1,
        specs=[[{}] for var in variables],
        shared_xaxes=True,
        vertical_spacing=0.001,
    )

    sg = 0
    for nv, var in enumerate(variables):
        if split is None:
            for n, cv in enumerate(groupvals):
                y = ad_here.obs[ad_here.obs[group] == cv][var]
                tr = go.Violin(
                    y=y,
                    name=cv,
                    fillcolor=cm[n % 40],
                    line={"color": "gray", "width": 0.5},
                    marker={"size": 1},
                    spanmode="manual",
                    span=[y.min(), None],
                    showlegend=nv == 0,
                    scalegroup=sg,
                )
                fig.append_trace(tr, nvar - nv, 1)
                sg += 1
            plot_data = ad_here.obs[variables + [group]]
        else:
            for n, sv in enumerate(splitvals):
                y = ad_here.obs[ad_here.obs[split] == sv][var]
                tr = go.Violin(
                    x=ad_here.obs[ad_here.obs[split] == sv][group],
                    y=y,
                    name=sv,
                    offsetgroup=n,
                    scalegroup=sg,
                    showlegend=nv == 0,
                    fillcolor=cm[n % 40],
                    line={"color": "gray", "width": 0.5},
                    marker={"size": 1},
                    spanmode="manual",
                    span=[y.min(), None],
                )
                fig.append_trace(tr, nvar - nv, 1)
                sg += 1
            plot_data = ad_here.obs[variables + [group, split]]

    for nv, var in enumerate(variables):
        if len(variables) - nv == 1:
            fig["layout"]["yaxis"].update(title=var)
        else:
            fig["layout"]["yaxis" + str(nvar - nv)].update(title=var)

    fig["layout"].update(
        xaxis={"title": group, "tickangle": -45},
        margin={"l": 50, "b": 80, "t": 10, "r": 10},
        legend={"x": 1.05, "y": 1},
        hovermode="closest",
        plot_bgcolor="rgb(255,255,255)",
        height=settings.PLOT_HEIGHT,
    )

    if split is not None:
        fig["layout"].update(violinmode="group")

    csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(
        plot_data.to_csv(index=True, header=True, encoding="utf-8")
    )

    return fig, csv_string, False


def render_plot_bars(data, group, split, options, filters_json):
    """render the bar chart plot."""

    if group is None:
        return {}, "", True

    if options is None:
        options = []

    ad_here = common.apply_filter_cells_filters(data, filters_json)

    if split is None:
        groupvals = ad_here.obs[group].unique()
        cm = colors.get_cm(groupvals)
        tally = ad_here.obs.groupby(group).size()
        if "normalized" in options:
            tally = tally.divide(tally.sum())
        traces = []
        for n, gv in enumerate(tally.index):
            tr = go.Bar(
                x=[gv],
                y=[tally[gv]],
                name=gv,
                marker=dict(line=dict(color="gray", width=0.5), color=cm[n % 40]),
            )
            traces.append(tr)

    else:
        splitvals = ad_here.obs[split].cat.categories
        cm = colors.get_cm(splitvals)

        tally = ad_here.obs.groupby([group, split]).size()
        if "normalized" in options:
            tally = tally.divide(tally.sum(level=0).astype(float), level=0)
        traces = []
        for n, sv in enumerate(splitvals):
            tr = go.Bar(
                x=tally.xs(sv, level=1).index,
                y=tally.xs(sv, level=1).values,
                name=sv,
                marker=dict(color=cm[n % 40], line=dict(color="gray", width=0.5)),
            )
            traces.append(tr)

    fig = {
        "data": traces,
        "layout": go.Layout(
            barmode="stack" if "stacked" in options else "group",
            xaxis={"title": group, "tickangle": -45},
            yaxis={"title": "cell frequency" if "normalized" in options else "cell number"},
            margin={"l": 50, "b": 100, "t": 10, "r": 10},
            legend={"x": 1.05, "y": 1},
            plot_bgcolor="rgb(255,255,255)",
            hovermode="closest",
            height=settings.PLOT_HEIGHT,
        ),
    }

    plot_data = tally
    csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(
        plot_data.to_csv(index=True, header=True, encoding="utf-8")
    )

    return fig, csv_string, False


def run_differential_expression(data, select_json):

    ad_here = data.ad

    selected = json.loads(select_json)

    ad_here.obs["group"] = np.nan
    ad_here.obs.loc[selected["group_A"], "group"] = "A"
    ad_here.obs.loc[selected["group_B"], "group"] = "B"

    ad_here = ad_here[~ad_here.obs["group"].isnull(), :]

    res = sc.tl.rank_genes_groups(ad_here, "group", copy=True).uns["rank_genes_groups"]

    res_df = pd.DataFrame(
        {
            "gene": pd.DataFrame(res["names"]).stack(),
            "logFC": pd.DataFrame(res["logfoldchanges"]).stack(),
            "pval": pd.DataFrame(res["pvals"]).stack(),
            "adjp": pd.DataFrame(res["pvals_adj"]).stack(),
        }
    )

    res_df.index = res_df.index.set_names(["n", "group"])
    res_df = res_df.reset_index().drop("n", axis=1)
    res_df = res_df[res_df["adjp"] < 0.05]

    results_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(
        res_df.to_csv(index=True, header=True, encoding="utf-8")
    )

    params = res["params"]
    params["group_A"] = ",".join(selected["group_A"])
    params["group_B"] = ",".join(selected["group_B"])

    params_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(
        pd.Series(params).to_csv(index=True, header=False, encoding="utf-8")
    )

    return res_df.to_json(), results_string, params_string, {"display": "block"}
