"""Render the plots and controls for the gene tab pane."""

import urllib.parse

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import plotly.tools as tools

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
                    options=[{"label": c, "value": c} for c in data.numerical_meta],
                    value=data.meta.columns[0],
                ),
                html.Label("select y axis"),
                dcc.Dropdown(
                    id="meta_scatter_select_y",
                    options=[{"label": c, "value": c} for c in data.numerical_meta],
                    value=data.meta.columns[1],
                ),
                html.Label("select coloring"),
                dcc.Dropdown(
                    id="meta_scatter_select_color",
                    options=[{"label": c, "value": c} for c in data.meta.columns],
                    value=data.categorical_meta[0],
                ),
            ],
            title="Select x- and y-coordinates for embedding (TSNE or UMAP) "
            "and color points according to cell annotation (e.g., cluster identity or n_genes)",
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
                        {
                            "label": "normalized",
                            "value": "normalized",
                            "title": "plot fractions instead of cell numbers",
                        },
                        {
                            "label": "stacked",
                            "value": "stacked",
                            "title": "stack bars instead of side-by-side",
                        },
                    ],
                    values=[],
                ),
            ],
            title=(
                "Select group for which cell numbers are tallied; "
                "these groups can optionally be split by another variable for sub-tallies; "
                "bars can be normalized and/or stacked."
            ),
        )
    ]


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
        # Placeholder for the plot-specific controls.
        dcc.Loading(id="meta_plot_controls", type="circle"),
        html.Hr(),
        # Control for sub-sampling of cells.
        common.render_subsampling_dropdown(data, "meta"),
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


def render_plot_scatter(data, xc, yc, col, sample_size):
    """Render the scatter plot figure."""
    if xc is None or yc is None or col is None:
        return {}, "", True

    if sample_size != "all":
        meta_here = data.meta.sample(n=sample_size)
    else:
        meta_here = data.meta

    if col in data.categorical_meta:
        # select color palette
        colvals = meta_here[col].unique()
        cm = colors.get_cm(colvals)

        # plot scatter for each category separately
        traces = []
        for n, cv in enumerate(colvals):
            traces.append(
                go.Scattergl(
                    x=meta_here[meta_here[col] == cv][xc],
                    y=meta_here[meta_here[col] == cv][yc],
                    text=meta_here[meta_here[col] == cv].index,
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
            # extra trace for legend with bigger marker
    #            traces.append(
    #                go.Scatter(
    #                    x=[None],
    #                    y=[None],
    #                    mode="markers",
    #                    marker={
    #                        "size": 20,
    #                        "color": cm[n % 40],
    #                        "line": {"width": 0.1, "color": "gray"},
    #                    },
    #                    showlegend=True,
    #                    visible="legendonly",
    #                    name=cv,
    #                )
    #            )
    else:
        # for numerical data, plot scatter all at once
        traces = [
            go.Scattergl(
                x=meta_here[xc],
                y=meta_here[yc],
                text=meta_here.index,
                mode="markers",
                opacity=0.7,
                marker={
                    "size": 5,
                    "color": meta_here[col],
                    "colorscale": "Viridis",
                    "colorbar": {"title": col, "titleside": "right"},
                    "showscale": True,
                },
                name=col,
            )
        ]

    plot_data = data.meta[[xc, yc, col]]
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
            hovermode="closest",
            height=settings.PLOT_HEIGHT,
        ),
    }
    return fig, csv_string, False


def render_plot_violin(data, variables, group, split, sample_size):
    """Render the violin plot figure."""
    if variables is None or len(variables) == 0:
        return {}, "", True

    if sample_size != "all":
        meta_here = data.meta.sample(n=sample_size)
    else:
        meta_here = data.meta

    # select color palette
    if split is None:
        groupvals = meta_here[group].unique()
        cm = colors.get_cm(groupvals)
    else:
        splitvals = meta_here[split].unique()
        cm = colors.get_cm(splitvals)

    nvar = len(variables)

    fig = tools.make_subplots(
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
                y = meta_here[meta_here[group] == cv][var]
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
            plot_data = data.meta[variables + [group]]
        else:
            for n, sv in enumerate(splitvals):
                y = meta_here[meta_here[split] == sv][var]
                tr = go.Violin(
                    x=meta_here[meta_here[split] == sv][group],
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
            plot_data = data.meta[variables + [group, split]]

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
        height=settings.PLOT_HEIGHT,
    )

    if split is not None:
        fig["layout"].update(violinmode="group")

    csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(
        plot_data.to_csv(index=True, header=True, encoding="utf-8")
    )

    return fig, csv_string, False


def render_plot_bars(data, group, split, options):
    """render the bar chart plot."""

    if group is None:
        return {}, "", True

    if split is None:
        groupvals = data.meta[group].unique()
        cm = colors.get_cm(groupvals)
        tally = data.meta.groupby(group).size()
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
        splitvals = data.meta[split].cat.categories
        cm = colors.get_cm(splitvals)

        tally = data.meta.groupby([group, split]).size()
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
            hovermode="closest",
            height=settings.PLOT_HEIGHT,
        ),
    }

    plot_data = tally
    csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(
        plot_data.to_csv(index=True, header=True, encoding="utf-8")
    )

    return fig, csv_string, False
