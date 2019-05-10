"""Implementation of the app callbacks."""

import os.path
import urllib.parse

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.tools as tools

from .ui.colors import get_cm, my_gradients
from .layout import HOME_BRAND, render_dataset
from .ui import cells, genes
from .settings import PLOT_HEIGHT
from . import store


def get_page(pathname):
    """Helper function for parsing the URL path."""
    pathname = pathname or "/"
    if pathname in ("/", "/home"):
        return {"page": "home"}
    else:
        tokens = pathname.split("/")
        if len(tokens) < 3:
            return {"page": None}
        else:
            return {"page": "viz", "dataset": tokens[2]}


def display_home():
    """Return site content for the home screen."""
    with open(os.path.join(os.path.dirname(__file__), "static", "home.md")) as inputf:
        home_md = inputf.read()
    return dbc.Row(dbc.Col(html.Div(dcc.Markdown(home_md))))


def display_not_found():
    """Return site content in the case that the dataset could not be found."""
    return html.Div(
        children=[
            html.H3("Dataset Not Found"),
            html.P("The dataset that you specified could not be found!"),
        ]
    )


def display_dataset(identifier):
    """Display the dataset."""
    try:
        data = store.load_data(identifier)
    except FileNotFoundError:
        return ""
    else:
        return render_dataset(data)


def register_page_content(app):
    """Register the display of the page content with the app."""

    @app.callback(
        dash.dependencies.Output("page-content", "children"),
        [dash.dependencies.Input("url", "pathname")],
    )
    def render_page_content(pathname):
        page = get_page(pathname)
        if page["page"] == "home":
            return display_home()
        elif page["page"] == "viz":
            return display_dataset(page["dataset"])
        else:
            return display_not_found()


def register_page_brand(app):
    """Register the display of the page brand with the app."""

    @app.callback(
        dash.dependencies.Output("page-navbar", "brand"),
        [dash.dependencies.Input("url", "pathname")],
    )
    def render_page_brand(pathname):
        page = get_page(pathname)
        if page["page"] in ("/", "/home"):
            return HOME_BRAND
        elif page["page"] == "viz":
            try:
                return "Dataset: %s" % store.load_metadata(page["dataset"]).title
            except FileNotFoundError:
                return "Dataset Not Found"
        else:
            return "Dataset Not Found"


def register_select_cell_plot_type(app):
    """Register callback for changing the controls on updating "CellAnnotation" plot type."""

    @app.callback(
        [dash.dependencies.Output("meta_plot_controls", "children")],
        [
            dash.dependencies.Input("url", "pathname"),
            dash.dependencies.Input("meta_plot_type", "value"),
        ],
    )
    def update_meta_plot_controls(pathname, plot_type):
        data = store.load_data(get_page(pathname)["dataset"])
        plots = {
            "scatter": cells.render_controls_scatter,
            "violin": cells.render_controls_violin,
            "bar": cells.render_controls_bars,
        }
        return plots[plot_type](data)

    @app.callback(
        [
            dash.dependencies.Output("meta_plot", "children"),
            dash.dependencies.Output("meta_select_cell_sample", "disabled"),
        ],
        [dash.dependencies.Input("meta_plot_type", "value")],
    )
    def update_meta_plots(plot_type):
        return cells.render_plot(plot_type)


def register_update_cell_scatter_plot_params(app):
    """Register handlers on updating scatter plot controls."""

    @app.callback(
        [
            dash.dependencies.Output("meta_scatter_plot", "figure"),
            dash.dependencies.Output("meta_scatter_download", "href"),
            dash.dependencies.Output("meta_scatter_download", "hidden"),
        ],
        [
            dash.dependencies.Input("url", "pathname"),
            dash.dependencies.Input("meta_scatter_select_x", "value"),
            dash.dependencies.Input("meta_scatter_select_y", "value"),
            dash.dependencies.Input("meta_scatter_select_color", "value"),
            dash.dependencies.Input("meta_select_cell_sample", "value"),
        ],
    )
    def get_meta_scatterplot(pathname, xc, yc, col, sample_size):
        data = store.load_data(get_page(pathname)["dataset"])
        return cells.render_plot_scatter(data, xc, yc, col, sample_size)


def register_update_cell_violin_plot_params(app):
    """Register handlers on updating violin plot controls."""

    @app.callback(
        [
            dash.dependencies.Output("meta_violin_plot", "figure"),
            dash.dependencies.Output("meta_violin_download", "href"),
            dash.dependencies.Output("meta_violin_download", "hidden"),
        ],
        [
            dash.dependencies.Input("url", "pathname"),
            dash.dependencies.Input("meta_violin_select_vars", "value"),
            dash.dependencies.Input("meta_violin_select_group", "value"),
            dash.dependencies.Input("meta_violin_select_split", "value"),
            dash.dependencies.Input("meta_select_cell_sample", "value"),
        ],
    )
    def get_meta_violinplot(pathname, variables, group, split, sample_size):
        data = store.load_data(get_page(pathname)["dataset"])
        return cells.render_plot_violin(data, variables, group, split, sample_size)


def register_update_cell_bar_chart_params(app):
    """Register handlers on updating bar chart controls."""

    @app.callback(
        dash.dependencies.Output("meta_bar_options", "options"),
        [dash.dependencies.Input("meta_bar_select_split", "value")],
    )
    def toggle_meta_bar_options(split):
        if split is None:
            return [{"label": "normalized", "value": "normalized"}]
        else:
            return [{"label": c, "value": c} for c in ["normalized", "stacked"]]


def register_select_gene_plot_type(app):
    """Register callback for changing the controls on updating "Gene Expression" plot type."""

    @app.callback(
        [dash.dependencies.Output("expression_plot_controls", "children")],
        [
            dash.dependencies.Input("url", "pathname"),
            dash.dependencies.Input("expression_plot_type", "value"),
        ],
    )
    def update_expression_plot_controls(pathname, plot_type):
        data = store.load_data(get_page(pathname)["dataset"])
        plots = {
            "scatter": genes.render_controls_scatter,
            "violin": genes.render_controls_violin,
            "dot": genes.render_controls_dot,
        }
        return [plots[plot_type](data)]

    @app.callback(
        [
            dash.dependencies.Output("expression_plot", "children"),
            dash.dependencies.Output("expression_select_cell_sample", "disabled"),
        ],
        [dash.dependencies.Input("expression_plot_type", "value")],
    )
    def update_meta_plots(plot_type):
        return genes.render_plot(plot_type)


def register_select_gene_marker_list(app):
    """Register callbacks related to changing marker setting."""

    @app.callback(
        dash.dependencies.Output("expression_marker_list", "children"),
        [
            dash.dependencies.Input("url", "pathname"),
            dash.dependencies.Input("expression_toggle_marker_list", "values"),
        ],
    )
    def toggle_marker_list(pathname, values):
        data = store.load_data(get_page(pathname)["dataset"])
        if "markers" in values:
            return [
                dbc.Row([dbc.Col(html.Hr())]),
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                dbc.Button(
                                    "use selected genes", color="primary", id="marker_selection"
                                )
                            ],
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

    @app.callback(
        dash.dependencies.Output("expression_select_genes", "value"),
        [
            dash.dependencies.Input("url", "pathname"),
            dash.dependencies.Input("marker_selection", "n_clicks"),
        ],
        [
            dash.dependencies.State("marker_list", "selected_rows"),
            dash.dependencies.State("expression_toggle_marker_list", "values"),
        ],
    )
    def update_gene_selection(pathname, n_clicks, rows, values):
        data = store.load_data(get_page(pathname)["dataset"])
        if "markers" in values:
            return list(data.markers.iloc[rows]["gene"])
        else:
            return []


def register_select_gene_scatter_plot(app):
    """Register callbacks for updating scatter plot"""

    @app.callback(
        [
            dash.dependencies.Output("expression_scatter_plot", "figure"),
            dash.dependencies.Output("expression_scatter_download", "href"),
            dash.dependencies.Output("expression_scatter_download", "hidden"),
        ],
        [
            dash.dependencies.Input("url", "pathname"),
            dash.dependencies.Input("expression_scatter_select_x", "value"),
            dash.dependencies.Input("expression_scatter_select_y", "value"),
            dash.dependencies.Input("expression_select_genes", "value"),
            dash.dependencies.Input("expression_select_cell_sample", "value"),
        ],
    )
    def get_expression_scatterplot(pathname, xc, yc, genelist, sample_size):
        data = store.load_data(get_page(pathname)["dataset"])
        if genelist is None or len(genelist) == 0:
            return {}, "", True

        if sample_size != "all":
            DGE_here = data.DGE.sample(n=sample_size, axis=1)
            meta_here = data.meta.loc[DGE_here.columns]
        else:
            DGE_here = data.DGE
            meta_here = data.meta

        ngenes = len(genelist)
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
            subplot_titles=genelist,
        )

        # rescale all expression values to the same min and max
        maxval = max(1, DGE_here.loc[genelist].max().max())

        for ng, gene in enumerate(genelist):
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
            height=PLOT_HEIGHT,
        )

        data = data.DGE.loc[genelist].T.to_csv(index=True, header=True, encoding="utf-8")
        csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(data)

        return fig, csv_string, False


def register_select_gene_violin_plot(app):
    """Register callbacks for updating violin plot"""

    @app.callback(
        [
            dash.dependencies.Output("expression_violin_plot", "figure"),
            dash.dependencies.Output("expression_violin_download", "href"),
            dash.dependencies.Output("expression_violin_download", "hidden"),
        ],
        [
            dash.dependencies.Input("url", "pathname"),
            dash.dependencies.Input("expression_select_genes", "value"),
            dash.dependencies.Input("expression_select_cell_sample", "value"),
            dash.dependencies.Input("expression_violin_select_group", "value"),
            dash.dependencies.Input("expression_violin_select_split", "value"),
        ],
    )
    def get_expression_violinplot(pathname, genelist, sample_size, group, split):
        data = store.load_data(get_page(pathname)["dataset"])
        if genelist is None or len(genelist) == 0:
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
            cm = get_cm(groupvals)
        else:
            splitvals = meta_here[split].cat.categories
            cm = get_cm(splitvals)

        ngenes = len(genelist)

        fig = tools.make_subplots(
            rows=ngenes,
            cols=1,
            specs=[[{}] for gene in genelist],
            shared_xaxes=True,
            vertical_spacing=0.01,
        )
        sg = 0
        for ng, gene in enumerate(genelist):

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
                the_data = data.DGE.loc[genelist].T.join(data.meta[group])
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
                the_data = data.DGE.loc[genelist].T.join(data.meta[[group, split]])

        for ng, gene in enumerate(genelist):
            if ngenes - ng == 1:
                fig["layout"]["yaxis"].update(title=gene)
            else:
                fig["layout"]["yaxis" + str(ngenes - ng)].update(title=gene)

        fig["layout"].update(
            xaxis={"title": group, "tickangle": -45},
            margin={"l": 50, "b": 80, "t": 10, "r": 10},
            legend={"x": 1.05, "y": 1},
            hovermode="closest",
            height=PLOT_HEIGHT,
        )

        if split is not None:
            fig["layout"].update(violinmode="group")

        csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(
            the_data.to_csv(index=True, header=True, encoding="utf-8")
        )

        return fig, csv_string, False


def register_select_gene_dot_plot(app):
    """Register callbacks for updating dot plot"""

    @app.callback(
        [
            dash.dependencies.Output("expression_dot_plot", "figure"),
            dash.dependencies.Output("expression_dot_download", "href"),
            dash.dependencies.Output("expression_dot_download", "hidden"),
        ],
        [
            dash.dependencies.Input("url", "pathname"),
            dash.dependencies.Input("expression_select_genes", "value"),
            dash.dependencies.Input("expression_dot_select_group", "value"),
            dash.dependencies.Input("expression_dot_select_split", "value"),
        ],
    )
    def get_expression_dotplot(pathname, genelist, group, split):
        data = store.load_data(get_page(pathname)["dataset"])
        if genelist is None or len(genelist) == 0:
            return {}, "", True

        groupvals = data.meta[group].cat.categories
        ngroup = len(groupvals)
        ngenes = len(genelist)

        def my_agg(x):
            return pd.DataFrame([x.mean(), (x > 0).mean()], index=["expression", "pct_cells"])

        if split is None:
            data = (
                data.DGE.loc[genelist]
                .T.join(data.meta[group])
                .groupby(group)
                .apply(my_agg)
                .unstack(level=1)
            )
            means = data.xs("expression", axis=1, level=1).values
            sizes = data.xs("pct_cells", axis=1, level=1).values
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
                        "colorscale": my_gradients[0],
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
                    marker={"size": 20, "color": my_gradients[0][c][1], "symbol": "square"},
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
                    ticktext=genelist,
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
                height=PLOT_HEIGHT,
            )

        else:
            splitvals = data.eta[split].cat.categories
            nsplit = len(splitvals)
            traces = []
            data = (
                data.DGE.loc[genelist]
                .T.join(data.meta[[group, split]])
                .groupby([group, split])
                .apply(my_agg)
                .unstack(level=2)
            )
            for n, sv in enumerate(splitvals):
                means = data.xs("expression", axis=1, level=1).xs(sv, axis=0, level=1).values
                sizes = data.xs("pct_cells", axis=1, level=1).xs(sv, axis=0, level=1).values
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
                        "colorscale": my_gradients[n],
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
                    marker={"size": 20, "symbol": "square", "color": my_gradients[-1][c][1]},
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
                    marker={"size": 20, "symbol": "square", "color": my_gradients[n][1][1]},
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
                    ticktext=genelist,
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
                height=PLOT_HEIGHT,
            )

        csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(
            data.to_csv(index=True, header=True, encoding="utf-8")
        )

        fig = {"data": traces, "layout": layout}
        return fig, csv_string, False
