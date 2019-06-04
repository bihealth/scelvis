"""Render the plots and controls for the gene tab pane."""

import urllib.parse

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.tools as tools

from .. import settings
from . import colors, common


def render_controls_scatter(data):
    """Render "top left" controls for the scatter plot for the given ``data``."""
    return [
        html.Div(
            [
                html.Label("select x axis"),
                dcc.Dropdown(
                    id="expression_scatter_select_x",
                    options=[{"label": c, "value": c} for c in data.numerical_meta],
                    value=data.meta.columns[0],
                ),
                html.Label("select y axis"),
                dcc.Dropdown(
                    id="expression_scatter_select_y",
                    options=[{"label": c, "value": c} for c in data.numerical_meta],
                    value=data.meta.columns[1],
                ),
            ],
            title="""select x- and y-coordinates of embedding (e.g., TSNE or UMAP)""",
        )
    ]


def render_controls_violin(data):
    """Render "top left" controls for the violin plot for the given ``data``."""
    return html.Div(
        [
            html.Label("select grouping"),
            dcc.Dropdown(
                id="expression_violin_select_group",
                options=[{"label": c, "value": c} for c in data.categorical_meta],
                value=data.categorical_meta[0],
            ),
            html.Label("select split"),
            dcc.Dropdown(
                id="expression_violin_select_split",
                options=[{"label": c, "value": c} for c in data.categorical_meta],
                value=None,
            ),
        ],
        title=(
            "Choose how to group cells for plotting gene expression distributions (e.g., by cluster); "
            "optionally how to split these groups (e.g., by genotype)"
        ),
    )


def render_controls_dot(data):
    """Render "top left" controls for the bubble plot for the given ``data``."""
    return html.Div(
        [
            html.Label("select grouping"),
            dcc.Dropdown(
                id="expression_dot_select_group",
                options=[{"label": c, "value": c} for c in data.categorical_meta],
                value=data.categorical_meta[0],
            ),
            html.Label("select split"),
            dcc.Dropdown(
                id="expression_dot_select_split",
                options=[
                    {"label": c, "value": c}
                    for c in data.categorical_meta
                    if len(data.meta[c].unique()) < 4
                ],
                value=None,
            ),
        ],
        title=(
            "Choose grouping of cells for y-coordinate of dot plot (e.g., by cluster); "
            "optionally split these groups (e.g., by genotype)"
        ),
    )


def render_controls(data):
    """Render the (left) controls column for the given ``data``."""
    return [
        # Select plot type.
        html.Div(
            [
                html.Label("select plot type"),
                dcc.RadioItems(
                    id="expression_plot_type",
                    options=[
                        {"label": "scatter plot", "value": "scatter"},
                        {"label": "violin plot", "value": "violin"},
                        {"label": "dot plot", "value": "dot"},
                    ],
                    value="scatter",
                    labelStyle={"display": "block"},
                    inputClassName="mr-1",
                ),
            ],
            title=(
                "SCATTER: gene expression on 2-dimensional embedding, one plot per gene; "
                "VIOLIN: gene expression distributions in different (sub)groups, one row per gene; "
                "DOT: summarized gene expression in (sub)groups for multiple genes"
            ),
        ),
        html.Hr(),
        # Placeholder for the plot-specific controls.
        dcc.Loading(id="expression_plot_controls", type="circle"),
        html.Hr(),
        # Selection of genes and markers.
        html.Div(
            [
                html.Label("select gene(s)"),
                dcc.Dropdown(
                    id="expression_select_genes",
                    options=[{"label": g, "value": g} for g in data.genes],
                    value=[],
                    multi=True,
                ),
                dcc.Checklist(
                    id="expression_toggle_marker_list",
                    options=[
                        {
                            "label": "use marker table for selection",
                            "value": "markers",
                            "disabled": (data.markers is None),
                        }
                    ],
                    values=[],
                    className="mt-2",
                    inputClassName="mr-1",
                ),
            ],
            title=(
                "Select one or more genes (remove from list by clicking on x) "
                "or select from marker gene list by checking the box"
            ),
        ),
        html.Hr(),
        # Control for sub-sampling of cells.
        common.render_subsampling_dropdown(data, "expression"),
    ]


def render(data):
    """Render the "Gene Expression" content."""
    return html.Div(
        children=[
            dbc.Row(
                children=[
                    dbc.Col(children=render_controls(data), className="col-3"),
                    # Placeholder for the plot.
                    dbc.Col(
                        children=[dcc.Loading(id="expression_plot", type="circle")],
                        className="col-9",
                    ),
                ]
            ),
            dbc.Row(
                children=[
                    dbc.Col(children=[dcc.Loading(id="expression_marker_list", type="circle")])
                ]
            ),
        ]
    )


def render_marker_list(data, values):
    if "markers" in values:
        return [
            dbc.Row([dbc.Col(html.Hr())]),
            dbc.Row(
                [
                    dbc.Col(
                        [dbc.Button("use selected genes", color="primary", id="marker_selection")],
                        className="col-3",
                    ),
                    dbc.Col(
                        [
                            dash_table.DataTable(
                                id="marker_list",
                                columns=[{"name": i, "id": i} for i in data.markers.columns],
                                data=data.markers.round(3).to_dict("rows"),
                                n_fixed_rows=1,
                                style_as_list_view=True,
                                sorting=True,
                                sorting_type="single",
                                row_selectable="multi",
                                pagination_mode=False,
                                selected_rows=[],
                                style_table={"maxHeight": 200, "overflowY": "scroll"},
                            )
                        ],
                        className="col-9",
                    ),
                ]
            ),
        ]

    else:
        return [html.P()]


def render_plot_scatter(data, xc, yc, genelist, sample_size):

    gl = [g for g in genelist if g in data.DGE.index]

    if gl is None or len(gl) == 0 or xc is None or yc is None:
        return {}, "", True

    if sample_size != "all":
        DGE_here = data.DGE.sample(n=sample_size, axis=1)
        meta_here = data.meta.loc[DGE_here.columns]
    else:
        DGE_here = data.DGE
        meta_here = data.meta

    ngenes = len(gl)
    if ngenes > 1:
        nrow = int(np.floor(np.sqrt(ngenes)))
        ncol = int(np.ceil(ngenes / nrow))
    else:
        nrow = 1
        ncol = 1

    fig = tools.make_subplots(
        rows=nrow,
        cols=ncol,
        specs=[[{} for c in range(ncol)] for r in range(nrow)],
        shared_xaxes=True,
        shared_yaxes=True,
        vertical_spacing=0.01,
        horizontal_spacing=0.01,
        subplot_titles=gl,
    )

    # rescale all expression values to the same min and max
    maxval = max(1, DGE_here.loc[gl].max().max())

    for ng, gene in enumerate(gl):
        # for numerical data, plot scatter all at once
        trace = go.Scatter(
            x=meta_here[xc],
            y=meta_here[yc],
            text=data.cells,
            mode="markers",
            marker={
                "size": 3,
                "color": maxval * DGE_here.loc[gene] / max(1, DGE_here.loc[gene].max()),
                "colorscale": "Viridis",
                "colorbar": {"title": "expression", "titleside": "right"},
                "showscale": (ng == 0),
            },
        )

        nc = ng % ncol + 1
        nr = ng // ncol + 1
        fig.append_trace(trace, nr, nc)

    # fix axes labels and scales
    for k in dir(fig["layout"]):
        if k.startswith("yaxis"):
            fig["layout"][k].update(title=yc)
        elif k.startswith("xaxis"):
            fig["layout"][k].update(title=xc)

    fig["layout"].update(
        margin={"l": 40, "b": 40, "t": 40, "r": 40},
        showlegend=False,
        hovermode="closest",
        height=settings.PLOT_HEIGHT,
    )

    plot_data = data.DGE.loc[gl].T.to_csv(index=True, header=True, encoding="utf-8")
    csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(plot_data)

    return fig, csv_string, False


def render_plot_violin(data, pathname, genelist, sample_size, group, split):

    gl = [g for g in genelist if g in data.DGE.index]

    if gl is None or len(gl) == 0 or group is None:
        return {}, "", True

    if sample_size != "all":
        DGE_here = data.DGE.sample(n=sample_size, axis=1)
        meta_here = data.meta.loc[DGE_here.columns]
    else:
        DGE_here = data.DGE
        meta_here = data.meta

    # select color palette
    if split is None:
        groupvals = meta_here[group].cat.categories
        cm = colors.get_cm(groupvals)
    else:
        splitvals = meta_here[split].cat.categories
        cm = colors.get_cm(splitvals)

    ngenes = len(gl)

    fig = tools.make_subplots(
        rows=ngenes, cols=1, specs=[[{}] for gene in gl], shared_xaxes=True, vertical_spacing=0.01
    )
    sg = 0
    for ng, gene in enumerate(gl):

        if split is None:
            for n, cv in enumerate(groupvals):
                y = DGE_here.loc[gene, meta_here[group] == cv]
                tr = go.Violin(
                    y=y,
                    name=cv,
                    fillcolor=cm[n % 40],
                    line={"color": "gray", "width": 0.5},
                    marker={"size": 1},
                    spanmode="manual",
                    span=[y.min(), None],
                    showlegend=ng == 0,
                    scalegroup=sg,
                )
                fig.append_trace(tr, ngenes - ng, 1)
                sg += 1
            plot_data = data.DGE.loc[gl].T.join(data.meta[group])
        else:
            for n, sv in enumerate(splitvals):
                y = DGE_here.loc[gene, meta_here[split] == sv]
                tr = go.Violin(
                    x=meta_here[meta_here[split] == sv][group],
                    y=y,
                    name=sv,
                    offsetgroup=n,
                    scalegroup=sg,
                    showlegend=ng == 0,
                    fillcolor=cm[n % 40],
                    line={"color": "gray", "width": 0.5},
                    spanmode="manual",
                    marker={"size": 1},
                    span=[y.min(), None],
                )
                fig.append_trace(tr, ngenes - ng, 1)
                sg += 1
            plot_data = data.DGE.loc[gl].T.join(data.meta[[group, split]])

    for ng, gene in enumerate(gl):
        if ngenes - ng == 1:
            fig["layout"]["yaxis"].update(title=gene)
        else:
            fig["layout"]["yaxis" + str(ngenes - ng)].update(title=gene)

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


def render_plot_dot(data, pathname, genelist, group, split):

    gl = [g for g in genelist if g in data.DGE.index]

    if gl is None or len(gl) == 0 or group is None:
        return {}, "", True

    groupvals = data.meta[group].cat.categories
    ngroup = len(groupvals)
    ngenes = len(gl)

    def my_agg(x):
        return pd.DataFrame([x.mean(), (x > 0).mean()], index=["expression", "pct_cells"])

    if split is None:
        plot_data = (
            data.DGE.loc[gl].T.join(data.meta[group]).groupby(group).apply(my_agg).unstack(level=1)
        )
        means = plot_data.xs("expression", axis=1, level=1).values
        sizes = plot_data.xs("pct_cells", axis=1, level=1).values
        xv, yv = np.meshgrid(np.arange(ngenes), np.arange(ngroup))
        traces = [
            go.Scatter(
                x=xv.ravel(),
                y=yv.ravel(),
                mode="markers",
                showlegend=False,
                opacity=1,
                marker={
                    "size": 20 * sizes.ravel(),
                    "color": means.ravel(),
                    "colorscale": colors.my_gradients[0],
                },
            )
        ]
        # extra invisible traces for legend
        traces += [
            go.Scatter(
                x=[0],
                y=[0],
                mode="markers",
                marker={"size": 20 * s, "color": "black"},
                name=p,
                visible="legendonly",
                legendgroup="size",
                showlegend=True,
            )
            for s, p in zip([0, 0.5, 1], ["0% cells", "50% cells", "100% cells"])
        ]
        traces += [
            go.Scatter(
                x=[0],
                y=[0],
                mode="markers",
                marker={"size": 20, "color": colors.my_gradients[0][c][1], "symbol": "square"},
                name=p,
                visible="legendonly",
                legendgroup="color",
                showlegend=True,
            )
            for c, p in zip([0, 1], ["low", "high"])
        ]
        layout = go.Layout(
            xaxis=go.layout.XAxis(
                showgrid=False,
                zeroline=False,
                showline=False,
                tickvals=np.arange(ngenes),
                ticktext=gl,
                tickangle=-90,
            ),
            yaxis=go.layout.YAxis(
                showgrid=False,
                zeroline=False,
                showline=False,
                tickvals=np.arange(ngroup),
                ticktext=groupvals,
            ),
            margin={"l": 150, "b": 100, "t": 40, "r": 40},
            showlegend=True,
            hovermode="closest",
            height=settings.PLOT_HEIGHT,
        )

    else:
        splitvals = data.meta[split].cat.categories
        nsplit = len(splitvals)
        traces = []
        plot_data = (
            data.DGE.loc[gl]
            .T.join(data.meta[[group, split]])
            .groupby([group, split])
            .apply(my_agg)
            .unstack(level=2)
        )
        for n, sv in enumerate(splitvals):
            means = plot_data.xs("expression", axis=1, level=1).xs(sv, axis=0, level=1).values
            sizes = plot_data.xs("pct_cells", axis=1, level=1).xs(sv, axis=0, level=1).values
            xv, yv = np.meshgrid(np.arange(ngenes), (nsplit + 1) * np.arange(ngroup))
            tr = go.Scatter(
                x=xv.ravel(),
                y=yv.ravel() + n,
                mode="markers",
                showlegend=False,
                opacity=1,
                marker={
                    "size": 20 * sizes.ravel(),
                    "color": means.ravel(),
                    "colorscale": colors.my_gradients[n],
                },
            )
            traces.append(tr)

        # extra invisible traces for legend
        traces += [
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker={"size": 20 * s, "color": "black"},
                name=p,
                visible="legendonly",
                legendgroup="size",
                showlegend=True,
            )
            for s, p in zip([0, 0.5, 1], ["0% cells", "50% cells", "100% cells"])
        ]
        traces += [
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker={"size": 20, "symbol": "square", "color": colors.my_gradients[-1][c][1]},
                name=p,
                visible="legendonly",
                legendgroup="color",
                showlegend=True,
            )
            for c, p in zip([0, 1], ["low", "high"])
        ]
        traces += [
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker={"size": 20, "symbol": "square", "color": colors.my_gradients[n][1][1]},
                name=sv,
                visible="legendonly",
                legendgroup="split",
                showlegend=True,
            )
            for n, sv in enumerate(splitvals)
        ]

        layout = go.Layout(
            xaxis=go.layout.XAxis(
                showgrid=False,
                zeroline=False,
                showline=False,
                tickvals=np.arange(ngenes),
                ticktext=gl,
                tickangle=-90,
            ),
            yaxis=go.layout.YAxis(
                showgrid=False,
                zeroline=False,
                showline=False,
                tickvals=(nsplit + 1) * np.arange(ngroup) + (nsplit - 1) / 2.0,
                ticktext=groupvals,
            ),
            margin={"l": 150, "b": 100, "t": 40, "r": 40},
            showlegend=True,
            hovermode="closest",
            height=settings.PLOT_HEIGHT,
        )

    csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(
        plot_data.to_csv(index=True, header=True, encoding="utf-8")
    )

    fig = {"data": traces, "layout": layout}
    return fig, csv_string, False
