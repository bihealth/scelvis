"""Render the plots and controls for the gene tab pane."""

import urllib.parse

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import pandas as pd
import numpy as np
import plotly.graph_objs as go
import plotly.subplots as subplots


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
                    options=[{"label": c, "value": c} for c in data.coords + data.numerical_meta],
                    value=data.coords[0] if len(data.coords) > 1 else data.numerical_meta[0],
                ),
                html.Label("select y axis"),
                dcc.Dropdown(
                    id="expression_scatter_select_y",
                    options=[{"label": c, "value": c} for c in data.coords + data.numerical_meta],
                    value=data.coords[1] if len(data.coords) > 1 else data.numerical_meta[1],
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


def render_controls_box(data):
    """Render "top left" controls for the box plot for the given ``data``."""
    return html.Div(
        [
            html.Label("select grouping"),
            dcc.Dropdown(
                id="expression_box_select_group",
                options=[{"label": c, "value": c} for c in data.categorical_meta],
                value=data.categorical_meta[0],
            ),
            html.Label("select split"),
            dcc.Dropdown(
                id="expression_box_select_split",
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
                    if len(data.ad.obs[c].unique()) < 4
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
                        {"label": "box plot", "value": "box"},
                        {"label": "dot plot", "value": "dot"},
                    ],
                    value="scatter",
                    labelStyle={"display": "block"},
                    inputClassName="mr-1",
                ),
            ],
            title=(
                "SCATTER: gene expression on 2-dimensional embedding, one plot per gene; "
                "VIOLIN/BOX: gene expression distributions in different (sub)groups, one row per gene; "
                "DOT: summarized gene expression in (sub)groups for multiple genes"
            ),
        ),
        html.Hr(),
        # Control for filtering of cells.
        common.render_filter_cells_collapse(data, "expression"),
        html.Hr(),
        # Placeholder for the plot-specific controls.
        html.Div(id="expression_plot_controls"),
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
                    id="expression_toggle_gene_list",
                    options=[
                        {
                            "label": "show marker table",
                            "value": "markers",
                            "disabled": (data.markers is None),
                        },
                        {"label": "show diffexp results", "value": "diffexp", "disabled": False},
                    ],
                    value=[],
                    className="mt-2",
                    inputClassName="mr-1",
                ),
            ],
            title=(
                "Select one or more genes (remove from list by clicking on x) "
                "or select from gene lists by checking the boxes"
            ),
        ),
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
                    dbc.Col(
                        children=[
                            html.Div(id="expression_marker_list"),
                            html.Div(id="expression_diffexp_list"),
                        ]
                    )
                ]
            ),
        ]
    )


def render_marker_list(data, value):
    show = value is not None and "markers" in value
    if show:
        columns = [{"name": i, "id": i} for i in data.markers.columns]
        table_data = data.markers.round(3).to_dict("rows")
    else:
        columns = []
        table_data = []
    return [
        dbc.Row([dbc.Col(html.Hr())], style={"display": "block" if show else "none"}),
        dbc.Row(
            [
                dbc.Col(
                    [dbc.Button("use selected genes", color="primary", id="marker_selection")],
                    className="col-3",
                    style={"display": "block" if show else "none"},
                ),
                dbc.Col(
                    [
                        dash_table.DataTable(
                            id="marker_list",
                            columns=columns,
                            data=table_data,
                            fixed_rows={"headers": True},
                            style_as_list_view=True,
                            sort_action="native",
                            sort_mode="single",
                            row_selectable="multi",
                            selected_rows=[],
                            style_table={"maxHeight": 200, "overflowY": "scroll"},
                        )
                    ],
                    className="col-9",
                    style={"display": "block" if show else "none"},
                ),
            ]
        ),
    ]


def render_diffexp_list(value, diffexp_json):
    show = value is not None and "diffexp" in value
    if show:
        diffexp = pd.read_json(diffexp_json)
        columns = [{"name": i, "id": i} for i in diffexp.columns]
        table_data = diffexp.round(3).to_dict("rows")
    else:
        columns = []
        table_data = []
    return [
        dbc.Row([dbc.Col(html.Hr())], style={"display": "block" if show else "none"}),
        dbc.Row(
            [
                dbc.Col(
                    [dbc.Button("use selected genes", color="primary", id="diffexp_selection")],
                    className="col-3",
                    style={"display": "block" if show else "none"},
                ),
                dbc.Col(
                    [
                        dash_table.DataTable(
                            id="diffexp_list",
                            columns=columns,
                            data=table_data,
                            fixed_rows={"headers": True},
                            style_as_list_view=True,
                            sort_action="native",
                            sort_mode="single",
                            row_selectable="multi",
                            selected_rows=[],
                            style_table={"maxHeight": 200, "overflowY": "scroll"},
                        )
                    ],
                    className="col-9",
                    style={"display": "block" if show else "none"},
                ),
            ]
        ),
    ]


def render_plot_scatter(data, xc, yc, genelist, filters_json):

    gl = [g for g in genelist if g in data.ad.var_names]

    if gl is None or len(gl) == 0 or xc is None or yc is None:
        return {}, "", True

    ad_here = common.apply_filter_cells_filters(data, filters_json)

    ngenes = len(gl)
    if ngenes > 1:
        nrow = int(np.floor(np.sqrt(ngenes)))
        ncol = int(np.ceil(ngenes / nrow))
    else:
        nrow = 1
        ncol = 1

    fig = subplots.make_subplots(
        rows=nrow,
        cols=ncol,
        specs=[[{} for c in range(ncol)] for r in range(nrow)],
        shared_xaxes=True,
        shared_yaxes=True,
        vertical_spacing=0.05,
        horizontal_spacing=0.05,
        subplot_titles=gl,
    )

    # rescale all expression values to the same min and max
    maxval = max(1, ad_here[:, gl].X.max().max())

    for ng, gene in enumerate(gl):
        # for numerical data, plot scatter all at once
        trace = go.Scatter(
            x=ad_here.obs[xc],
            y=ad_here.obs[yc],
            text=ad_here.obs_names,
            mode="markers",
            marker={
                "size": 3,
                "color": maxval * ad_here.obs_vector(gene) / max(1, ad_here.obs_vector(gene).max()),
                "colorscale": "deep",
                "colorbar": {"title": "expression", "titleside": "right"},
                "showscale": (ng == 0),
            },
        )

        nc = ng % ncol + 1
        nr = ng // ncol + 1
        fig.append_trace(trace, nr, nc)
        if nc == 1:
            fig["layout"]["yaxis" + (str(ng + 1) if ng > 0 else "")].update(title=yc)
        if nr == nrow:
            fig["layout"]["xaxis" + (str(ng + 1) if ng > 0 else "")].update(title=xc)

    fig["layout"].update(
        margin={"l": 40, "b": 40, "t": 40, "r": 40},
        plot_bgcolor="rgb(255,255,255)",
        showlegend=False,
        hovermode="closest",
        height=settings.PLOT_HEIGHT,
    )

    plot_data = ad_here[:, gl].to_df().to_csv(index=True, header=True, encoding="utf-8")
    csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(plot_data)

    return fig, csv_string, False


def render_plot_violin_box(data, pathname, plot_type, genelist, group, split, filters_json):

    gl = [g for g in genelist if g in data.ad.var_names]

    if gl is None or len(gl) == 0 or group is None:
        return {}, "", True

    ad_here = common.apply_filter_cells_filters(data, filters_json)

    # select color palette
    if split is None:
        groupvals = ad_here.obs[group].unique()
        cm = colors.get_cm(groupvals)
    else:
        splitvals = ad_here.obs[split].unique()
        cm = colors.get_cm(splitvals)

    ngenes = len(gl)

    fig = subplots.make_subplots(
        rows=ngenes, cols=1, specs=[[{}] for gene in gl], shared_xaxes=True, vertical_spacing=0.01
    )
    sg = 0
    for ng, gene in enumerate(gl):

        if split is None:
            for n, cv in enumerate(groupvals):
                y = ad_here[ad_here.obs[group] == cv, :].obs_vector(gene)
                if plot_type == "violin":
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
                else:
                    tr = go.Box(
                        y=y,
                        name=cv,
                        fillcolor=cm[n % 40],
                        line={"color": "gray", "width": 0.5},
                        marker={"size": 1},
                        showlegend=ng == 0,
                    )
                fig.append_trace(tr, ngenes - ng, 1)
                sg += 1
            plot_data = ad_here[:, gl].to_df().join(ad_here.obs[group])
        else:
            for n, sv in enumerate(splitvals):
                y = ad_here[ad_here.obs[split] == sv, :].obs_vector(gene)
                if plot_type == "violin":
                    tr = go.Violin(
                        x=ad_here[ad_here.obs[split] == sv].obs[group],
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
                else:
                    tr = go.Box(
                        x=ad_here[ad_here.obs[split] == sv].obs[group],
                        y=y,
                        name=sv,
                        offsetgroup=n,
                        showlegend=ng == 0,
                        fillcolor=cm[n % 40],
                        line={"color": "gray", "width": 0.5},
                        marker={"size": 1},
                    )
                fig.append_trace(tr, ngenes - ng, 1)
                sg += 1
            plot_data = ad_here[:, gl].to_df().join(ad_here.obs[[group, split]])

    for ng, gene in enumerate(gl):
        if ngenes - ng == 1:
            fig["layout"]["yaxis"].update(title=gene)
        else:
            fig["layout"]["yaxis" + str(ngenes - ng)].update(title=gene)

    fig["layout"]["xaxis" + (str(ngenes) if ngenes > 1 else "")].update(title=group)

    fig["layout"].update(
        xaxis={"tickangle": -45},
        plot_bgcolor="rgb(255,255,255)",
        margin={"l": 50, "b": 80, "t": 10, "r": 10},
        legend={"x": 1.05, "y": 1},
        hovermode="closest",
        height=settings.PLOT_HEIGHT,
    )

    if split is not None:
        if plot_type == "violin":
            fig["layout"].update(violinmode="group")
        else:
            fig["layout"].update(boxmode="group")

    csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(
        plot_data.to_csv(index=True, header=True, encoding="utf-8")
    )

    return fig, csv_string, False


def render_plot_dot(data, pathname, genelist, group, split, filters_json):

    gl = [g for g in genelist if g in data.ad.var_names]

    if gl is None or len(gl) == 0 or group is None:
        return {}, "", True

    ad_here = common.apply_filter_cells_filters(data, filters_json)

    groupvals = ad_here.obs[group].unique()
    ngroup = len(groupvals)
    ngenes = len(gl)

    def my_agg(x):
        return pd.DataFrame([x.mean(), (x > 0).mean()], index=["expression", "pct_cells"])

    if split is None:
        plot_data = (
            ad_here[:, gl]
            .to_df()
            .join(ad_here.obs[group])
            .groupby(group)
            .apply(my_agg)
            .unstack(level=1)
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
                x=[None],
                y=[None],
                mode="markers",
                marker={"size": 20 * s, "color": "black"},
                name=p,
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
                marker={"size": 20, "color": colors.my_gradients[0][c][1], "symbol": "square"},
                name=p,
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
        splitvals = ad_here.obs[split].unique()
        nsplit = len(splitvals)
        traces = []
        # fix bug in pandas groupby when having more than one grouping variable
        ii = pd.MultiIndex.from_arrays(
            pd.core.reshape.util.cartesian_product([groupvals, splitvals])
        )
        plot_data = (
            ad_here[:, gl]
            .to_df()
            .join(ad_here.obs[[group, split]])
            .groupby([group, split])
            .apply(my_agg)
            .unstack(level=2)
            .reindex(ii)
            .fillna(0)
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
            plot_bgcolor="rgb(255,255,255)",
        )

    csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(
        plot_data.to_csv(index=True, header=True, encoding="utf-8")
    )

    fig = {"data": traces, "layout": layout}
    return fig, csv_string, False
