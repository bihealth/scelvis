"""Registeration of Dash app callbacks.

Each ``register_*`` function sets up callbacks for one aspect of the app.  This module does not
any component building itself.  Instead, this is done in the module ``.ui``.
"""

import base64
import os.path
import uuid
import json

import dash
import dash_html_components as html
from dash.dependencies import Input, Output, State
import pandas as pd
from logzero import logger
from werkzeug.utils import secure_filename

from . import ui, settings, store


def get_route(pathname):
    """Helper function for routing URLs.

    >>> get_route("")
    "home", {}
    >>> get_route("/home")
    "home", {}
    >>> get_route("/")
    "home", {}
    >>> get_route("/viz/dataset")
    "home", {"dataset": "dataset"}
    >>> get_route("/bogus")
    None, {}

    Args:

    :pathname: The path to the current page.

    Returns:

    :page:   The name of the page to render.
    :kwargs: The arguments to render with.
    """
    pathname = pathname or "/dash/"
    suffixes = ("dash", "dash/", "dash/home", "dash/home/")
    prefixes = ("", settings.PUBLIC_URL_PREFIX + "/")
    home_urls = ["%s%s" % (prefix, x) for prefix in prefixes for x in suffixes]
    if pathname in home_urls:
        return "home", {}
    else:
        tokens = pathname.split("/")
        if len(tokens) < 3:
            return None, {}
        elif tokens[2] == "upload":
            return "upload", {}
        elif tokens[2] == "viz":
            if len(tokens) < 4:
                return None, {}
            else:
                return "viz", {"dataset": tokens[3]}
        else:
            return None, {}


def register_page_content(app):
    """Register the display of the page content with the app."""

    @app.callback(Output("page-content", "children"), [Input("url", "pathname")])
    def render_page_content(pathname):
        view, kwargs = get_route(pathname)
        if view == "home":
            return ui.main.render_home()
        elif view == "upload":
            return ui.main.render_upload()
        elif view == "viz":
            return ui.main.render_dataset(kwargs.get("dataset"))
        else:
            return ui.main.render_not_found()


def register_page_brand(app):
    """Register the display of the page brand with the app."""

    @app.callback(Output("page-brand", "children"), [Input("url", "pathname")])
    def render_page_brand(pathname):
        view, kwargs = get_route(pathname)
        if view == "home":
            return [html.I(className="fas fa-home mr-1"), settings.HOME_BRAND]
        if view == "upload":
            return [html.I(className="fas fa-cloud-upload-alt mr-1"), "Upload Data"]
        elif view == "viz":
            metadata = store.load_metadata(kwargs.get("dataset"))
            if metadata:
                return [html.I(className="fas fa-file-alt mr-1"), metadata.title]
            else:
                return [html.I(className="fas fa-file-alt mr-1"), "Not Found"]
        else:
            return [html.I(className="fas fa-file-alt mr-1"), "Not Found"]


def register_select_cell_plot_type(app):
    """Register callback for changing the controls on updating "CellAnnotation" plot type."""

    @app.callback(
        [Output("meta_plot_controls", "children")],
        [Input("url", "pathname"), Input("meta_plot_type", "value")],
    )
    def update_meta_plot_controls(pathname, plot_type):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        plots = {
            "scatter": ui.cells.render_controls_scatter,
            "violin": ui.cells.render_controls_violin,
            "bar": ui.cells.render_controls_bars,
        }
        return plots[plot_type](data)

    @app.callback(Output("meta_plot", "children"), [Input("meta_plot_type", "value")])
    def update_meta_plots(plot_type):
        return ui.common.render_plot("meta", plot_type)


def register_update_cell_scatter_plot_params(app):
    """Register handlers on scatter plot."""

    @app.callback(
        [
            Output("meta_scatter_plot", "figure"),
            Output("meta_scatter_download", "href"),
            Output("meta_scatter_download", "hidden"),
        ],
        [
            Input("url", "pathname"),
            Input("meta_scatter_select_x", "value"),
            Input("meta_scatter_select_y", "value"),
            Input("meta_scatter_select_color", "value"),
            Input("filter_cells_filters", "children"),
        ],
    )
    def get_meta_plot_scatter(pathname, xc, yc, col, filters_json):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        return ui.cells.render_plot_scatter(data, xc, yc, col, filters_json)


def register_toggle_select_cells_controls(app):
    @app.callback(
        Output("select_cells_collapse", "is_open"),
        [Input("select_cells_button", "n_clicks")],
        [State("select_cells_collapse", "is_open")],
    )
    def toggle_select_cells_controls(n, is_open):
        if n:
            return not is_open
        return is_open


def register_update_select_cells_choices(app):
    @app.callback(
        Output("select_cells_choices", "children"),
        [
            Input("meta_scatter_plot", "selectedData"),
            Input("select_cells_group_A", "n_clicks"),
            Input("select_cells_group_B", "n_clicks"),
            Input("select_cells_reset", "n_clicks"),
        ],
        [State("select_cells_choices", "children")],
    )
    def update_select_cells_choices(
        selectedData, n_clicks_A, n_clicks_B, n_clicks_reset, choices_json
    ):
        ctx = dash.callback_context
        if choices_json is not None:
            choices = json.loads(choices_json)
        else:
            choices = {}
        if ctx.triggered and "group_A" in ctx.triggered[0]["prop_id"] and selectedData is not None:
            choices["group_A"] = [p["text"] for p in selectedData["points"]]
        elif (
            ctx.triggered and "group_B" in ctx.triggered[0]["prop_id"] and selectedData is not None
        ):
            choices["group_B"] = [p["text"] for p in selectedData["points"]]
        elif ctx.triggered and "reset" in ctx.triggered[0]["prop_id"]:
            choices = {}
        return json.dumps(choices)


def register_activate_select_cells_buttons(app):
    @app.callback(
        [
            Output("select_cells_group_A", "disabled"),
            Output("select_cells_group_B", "disabled"),
            Output("select_cells_reset", "disabled"),
            Output("select_cells_run", "disabled"),
        ],
        [Input("meta_scatter_plot", "selectedData"), Input("select_cells_choices", "children")],
    )
    def activate_select_cells_buttons(selectedData, choices_json):
        if choices_json is not None:
            choices = json.loads(choices_json)
        else:
            choices = {}
        disabled_A = (selectedData is None) | ("group_A" in choices)
        disabled_B = (selectedData is None) | ("group_B" in choices)
        disabled_reset = ("group_A" not in choices) & ("group_B" not in choices)
        disabled_run = ("group_A" not in choices) | ("group_B" not in choices)

        return (disabled_A, disabled_B, disabled_reset, disabled_run)


def register_run_differential_expression(app):
    @app.callback(
        [
            Output("select_cells_results", "children"),
            Output("select_cells_results_download", "href"),
            Output("select_cells_parameters_download", "href"),
            Output("select_cells_get_results", "style"),
        ],
        [
            Input("url", "pathname"),
            Input("select_cells_run", "n_clicks"),
            Input("select_cells_choices", "children"),
        ],
    )
    def run_differential_expression(pathname, n_clicks, select_json):
        ctx = dash.callback_context
        if ctx.triggered and "run" in ctx.triggered[0]["prop_id"] and n_clicks:
            _, kwargs = get_route(pathname)
            data = store.load_data(kwargs.get("dataset"))
            res = ui.cells.run_differential_expression(data, select_json)
            return res
        else:
            return "", "", "", {"display": "none"}

    @app.callback(
        [Output("main-tabs", "active_tab"), Output("expression_toggle_gene_list", "value")],
        [Input("select_cells_view", "n_clicks"), Input("select_cells_reset", "n_clicks")],
        [State("expression_toggle_gene_list", "value"), State("main-tabs", "active_tab")],
    )
    def switch_view(n1, n2, selected, at):
        ctx = dash.callback_context
        if ctx.triggered and "view" in ctx.triggered[0]["prop_id"]:
            if "diffexp" not in selected:
                selected.append("diffexp")
            return "tab-genes", selected
        elif ctx.triggered and "reset" in ctx.triggered[0]["prop_id"]:
            if "diffexp" in selected:
                selected.remove("diffexp")
            return at, selected


def register_update_cell_violin_plot_params(app):
    """Register handlers on violin plot."""

    @app.callback(
        [
            Output("meta_violin_plot", "figure"),
            Output("meta_violin_download", "href"),
            Output("meta_violin_download", "hidden"),
        ],
        [
            Input("url", "pathname"),
            Input("meta_violin_select_vars", "value"),
            Input("meta_violin_select_group", "value"),
            Input("meta_violin_select_split", "value"),
            Input("filter_cells_filters", "children"),
        ],
    )
    def get_meta_plot_violin(pathname, variables, group, split, choices_json):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        return ui.cells.render_plot_violin(data, variables, group, split, choices_json)


def register_update_cell_bar_chart_params(app):
    """Register handlers on bar chart and its controls."""

    @app.callback(Output("meta_bar_options", "options"), [Input("meta_bar_select_split", "value")])
    def toggle_meta_bar_options(split):
        if split is None:
            return [{"label": "normalized", "value": "normalized"}]
        else:
            return [
                {"label": "normalized", "value": "normalized"},
                {"label": "stacked", "value": "stacked"},
            ]

    @app.callback(
        [
            Output("meta_bar_plot", "figure"),
            Output("meta_bar_download", "href"),
            Output("meta_bar_download", "hidden"),
        ],
        [
            Input("url", "pathname"),
            Input("meta_bar_select_group", "value"),
            Input("meta_bar_select_split", "value"),
            Input("meta_bar_options", "value"),
            Input("filter_cells_filters", "children"),
        ],
    )
    def get_meta_plot_bars(pathname, group, split, options, choices_json):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        return ui.cells.render_plot_bars(data, group, split, options, choices_json)


def register_select_gene_plot_type(app):
    """Register callback for changing the controls on updating "Gene Expression" plot type."""

    @app.callback(
        [Output("expression_plot_controls", "children")],
        [Input("url", "pathname"), Input("expression_plot_type", "value")],
    )
    def update_expression_plot_controls(pathname, plot_type):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        plots = {
            "scatter": ui.genes.render_controls_scatter,
            "violin": ui.genes.render_controls_violin,
            "dot": ui.genes.render_controls_dot,
        }
        return [plots[plot_type](data)]

    @app.callback(Output("expression_plot", "children"), [Input("expression_plot_type", "value")])
    def update_expression_plots(plot_type):
        return ui.common.render_plot("expression", plot_type)


def register_select_gene_list(app):
    """Register callbacks related to changing gene list settings."""

    @app.callback(
        Output("expression_marker_list", "children"),
        [Input("url", "pathname"), Input("expression_toggle_gene_list", "value")],
    )
    def toggle_marker_list(pathname, value):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        return ui.genes.render_marker_list(data, value)

    @app.callback(
        Output("expression_toggle_gene_list", "options"),
        [Input("select_cells_results", "children")],
        [State("expression_toggle_gene_list", "options")],
    )
    def activate_toggle_gene_list(diffexp_json, options):
        options[1]["disabled"] = diffexp_json is None or len(diffexp_json) == 0
        return options

    @app.callback(
        Output("expression_diffexp_list", "children"),
        [Input("expression_toggle_gene_list", "value")],
        [State("select_cells_results", "children")],
    )
    def toggle_gene_list(value, diffexp_json):
        return ui.genes.render_diffexp_list(value, diffexp_json)

    @app.callback(
        Output("expression_select_genes", "value"),
        [
            Input("url", "pathname"),
            Input("marker_selection", "n_clicks"),
            Input("diffexp_selection", "n_clicks"),
        ],
        [
            State("marker_list", "selected_rows"),
            State("diffexp_list", "selected_rows"),
            State("expression_toggle_gene_list", "value"),
            State("expression_select_genes", "value"),
            State("select_cells_results", "children"),
        ],
    )
    def update_gene_selection(
        pathname,
        n_clicks_markers,
        n_clicks_diffexp,
        rows_markers,
        rows_diffexp,
        selected_tables,
        selected_genes,
        diffexp_json,
    ):
        genelist = selected_genes
        if "markers" in selected_tables:
            _, kwargs = get_route(pathname)
            data = store.load_data(kwargs.get("dataset"))
            genelist += list(data.markers.iloc[rows_markers]["gene"])
        if "diffexp" in selected_tables and diffexp_json is not None:
            diffexp = pd.read_json(diffexp_json)
            genelist += list(diffexp.iloc[rows_diffexp]["gene"])
        return list(set(genelist))


def register_select_gene_scatter_plot(app):
    """Register callbacks for updating scatter plot"""

    @app.callback(
        [
            Output("expression_scatter_plot", "figure"),
            Output("expression_scatter_download", "href"),
            Output("expression_scatter_download", "hidden"),
        ],
        [
            Input("url", "pathname"),
            Input("expression_scatter_select_x", "value"),
            Input("expression_scatter_select_y", "value"),
            Input("expression_select_genes", "value"),
            Input("filter_cells_filters", "children"),
        ],
    )
    def get_expression_plot_scatter(pathname, xc, yc, genelist, choices_json):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        return ui.genes.render_plot_scatter(data, xc, yc, genelist, choices_json)


def register_select_gene_violin_plot(app):
    """Register callbacks for updating violin plot"""

    @app.callback(
        [
            Output("expression_violin_plot", "figure"),
            Output("expression_violin_download", "href"),
            Output("expression_violin_download", "hidden"),
        ],
        [
            Input("url", "pathname"),
            Input("expression_select_genes", "value"),
            Input("expression_violin_select_group", "value"),
            Input("expression_violin_select_split", "value"),
            Input("filter_cells_filters", "children"),
        ],
    )
    def get_expression_plot_violin(pathname, genelist, group, split, choices_json):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        return ui.genes.render_plot_violin(data, pathname, genelist, group, split, choices_json)


def register_select_gene_dot_plot(app):
    """Register callbacks for updating dot plot"""

    @app.callback(
        [
            Output("expression_dot_plot", "figure"),
            Output("expression_dot_download", "href"),
            Output("expression_dot_download", "hidden"),
        ],
        [
            Input("url", "pathname"),
            Input("expression_select_genes", "value"),
            Input("expression_dot_select_group", "value"),
            Input("expression_dot_select_split", "value"),
            Input("filter_cells_filters", "children"),
        ],
    )
    def get_expression_plot_dot(pathname, genelist, group, split, choices_json):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        return ui.genes.render_plot_dot(data, pathname, genelist, group, split, choices_json)


def register_toggle_filter_cells_controls(app, token):
    @app.callback(
        Output("%s_filter_cells_collapse" % token, "is_open"),
        [Input("%s_filter_cells_button" % token, "n_clicks")],
        [State("%s_filter_cells_collapse" % token, "is_open")],
    )
    def toggle_filter_cells_controls(n, is_open):
        if n:
            return not is_open
        return is_open


def register_update_filter_cells_controls(app, token):
    @app.callback(
        [
            Output("%s_filter_cells_choice_div" % token, "style"),
            Output("%s_filter_cells_choice" % token, "options"),
            Output("%s_filter_cells_choice" % token, "value"),
            Output("%s_filter_cells_range_div" % token, "style"),
            Output("%s_filter_cells_range" % token, "marks"),
            Output("%s_filter_cells_range" % token, "min"),
            Output("%s_filter_cells_range" % token, "max"),
            Output("%s_filter_cells_range" % token, "value"),
            Output("%s_filter_cells_range" % token, "step"),
        ],
        [Input("url", "pathname"), Input("%s_filter_cells_attribute" % token, "value")],
        [State("filter_cells_filters", "children")],
    )
    def update_filter_cells_controls(pathname, attribute, filters_json):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        if attribute is None or attribute == "None":
            return (
                {"display": "none"},
                [],
                None,
                {"display": "none"},
                {0: "0", 1: "1"},
                0,
                1,
                [0, 1],
                0,
            )
        filters = json.loads(filters_json)
        values = data.ad.obs_vector(attribute)
        if not pd.api.types.is_numeric_dtype(values):
            categories = list(data.ad.obs[attribute].cat.categories)
            return (
                {"display": "block"},
                [{"label": v, "value": v} for v in categories],
                filters[attribute] if attribute in filters else categories,
                {"display": "none"},
                {0: "0", 1: "1"},
                0,
                1,
                [0, 1],
                0,
            )
        else:
            range_min = values.min()
            range_max = values.max()
            if attribute in filters:
                val_min = float(filters[attribute][0])
                val_max = float(filters[attribute][1])
            else:
                val_min = range_min
                val_max = range_max
            return (
                {"display": "none"},
                [],
                None,
                {"display": "block"},
                dict(
                    (int(t) if t % 1 == 0 else t, "{0:g}".format(t))
                    for t in ui.common.auto_tick([range_min, range_max], max_tick=4, tf_inside=True)
                ),
                range_min,
                range_max,
                [val_min, val_max],
                (range_max - range_min) / 1000,
            )


def register_update_filter_cells_filters(app):
    @app.callback(
        [
            Output("filter_cells_filters", "children"),
            Output("meta_filter_cells_div", "title"),
            Output("expression_filter_cells_div", "title"),
        ],
        [
            Input("url", "pathname"),
            Input("meta_filter_cells_choice", "value"),
            Input("meta_filter_cells_range", "value"),
            Input("meta_filter_cells_reset", "n_clicks"),
            Input("expression_filter_cells_choice", "value"),
            Input("expression_filter_cells_range", "value"),
            Input("expression_filter_cells_reset", "n_clicks"),
        ],
        [
            State("meta_filter_cells_attribute", "value"),
            State("expression_filter_cells_attribute", "value"),
            State("filter_cells_filters", "children"),
        ],
    )
    def update_filter_cells_filters(
        pathname,
        meta_cat_value,
        meta_range_value,
        meta_reset_n,
        expression_cat_value,
        expression_range_value,
        expression_reset_n,
        meta_attribute,
        expression_attribute,
        filters_json,
    ):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        ctx = dash.callback_context

        filters = json.loads(filters_json)
        status = "active filters: "
        # if reset button was hit, check all boxes using stored values in filters_json
        attributes = filters.keys()
        if ctx.triggered and "reset" in ctx.triggered[0]["prop_id"]:
            for attribute in list(attributes):
                del filters[attribute]
            return (json.dumps(filters), status, status)

        for cat_value, range_value, attribute in [
            (meta_cat_value, meta_range_value, meta_attribute),
            (expression_cat_value, expression_range_value, expression_attribute),
        ]:
            if attribute is not None and attribute != "None":
                values = data.ad.obs_vector(attribute)
                if not pd.api.types.is_numeric_dtype(values):
                    filters[attribute] = sorted(cat_value)
                else:
                    filters[attribute] = ["{0:g}".format(v) for v in range_value]

        status += ", ".join(attributes)
        return (json.dumps(filters), status, status)


def register_activate_filter_cells_reset(app):
    @app.callback(
        [
            Output("meta_filter_cells_reset", "disabled"),
            Output("expression_filter_cells_reset", "disabled"),
        ],
        [Input("url", "pathname"), Input("filter_cells_filters", "children")],
    )
    def activate_filter_cells_reset(pathname, filters_json):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        if filters_json is not None:
            filters = json.loads(filters_json)
        else:
            filters = {}
        disabled = True
        attributes = filters.keys()
        for attribute in attributes:
            values = data.ad.obs_vector(attribute)
            if not pd.api.types.is_numeric_dtype(values):
                if filters[attribute] != list(data.ad.obs[attribute].cat.categories):
                    disabled = False
            else:
                range_min = values.min()
                range_max = values.max()
                val_min = float(filters[attribute][0])
                val_max = float(filters[attribute][1])
                if val_min > range_min or val_max < range_max:
                    disabled = False

        return (disabled, disabled)


def register_file_upload(app):
    """Register callbacks for the file upload"""

    logger.info("-> file_upload")

    # TODO: move main handler into "upload" module?
    @app.callback(
        Output("url", "pathname"),
        [Input("file-upload", "contents")],
        [State("file-upload", "filename")],
    )
    def file_uploaded(contents, filename):
        logger.info("Handling upload of %s", filename)
        filename = secure_filename(filename)
        _content_type, content_string = contents.split(",")
        data_uuid = str(uuid.uuid4())
        logger.info("Data will have UUID %s", data_uuid)
        # Decode base64 string and write out to final file.
        filepath = os.path.join(settings.UPLOAD_DIR, "%s.h5ad" % data_uuid)
        logger.info("Writing to %s", filepath)
        with open(filepath, "wb") as tmpf:
            tmpf.write(base64.b64decode(content_string))
        # Redirect to view the data set.
        return "/dash/viz/%s" % data_uuid
