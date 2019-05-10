"""The UI for the gene annotation."""

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

from ..exceptions import ScelVisException


def render_controls_scatter(data):
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
            "Choose how to group cells for plotting gene expression distributions (e.g., by cluster) and optionally "
            "how to split these groups (e.g., by genotype)"
        ),
    )


def render_controls_dot(data):
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
            "Choose grouping of cells for y-coordinate of dot plot (e.g., by cluster); optionally split these "
            "groups (e.g., by genotype)"
        ),
    )


def render_controls(data):
    """Render the controls column for the given ``data``."""
    # Choices for sub-sampling.
    sample_choices = [{"label": g, "value": g} for g in [100, 1000, 5000] if g < data.meta.shape[0]]
    sample_choices.append({"label": "all", "value": "all"})

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
                ),
            ],
            title=(
                "scatter: gene expression on 2-dimensional embedding, one plot per gene "
                "violin: gene expression distributions in different (sub)groups, one row per gene "
                "dot: summarized gene expression in (sub)groups for multiple genes"
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
                "Select one or more genes (remove from list by clicking on x) or select from marker gene list "
                "by checking the box"
            ),
        ),
        html.Hr(),
        # Control for sub-sampling of cells.
        html.Div(
            [
                html.Label("select cell sample"),
                dcc.Dropdown(
                    id="expression_select_cell_sample",
                    options=sample_choices,
                    value="all",
                    disabled=False,
                ),
            ],
            title="use random sample of cells for faster rendering",
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
                    dbc.Col(children=[dcc.Loading(id="expression_marker_list", type="circle")])
                ]
            ),
        ]
    )


def render_plot(plot_type):
    """Render the "Gene Expression" plot (placeholders)."""

    if plot_type not in ("scatter", "violin", "dot"):
        raise ScelVisException("Invalid plot type: %s", plot_type)

    return (
        [
            dcc.Graph(id="expression_%s_plot" % plot_type),
            html.A(
                children=[
                    html.I(className="fas fa-cloud-download-alt pr-1"),
                    "download data for this plot",
                ],
                id="expression_%s_download" % plot_type,
                download="plot_data.csv",
                href="",
                hidden=True,
                target="_blank",
            ),
        ],
        plot_type == "dot",
    )
