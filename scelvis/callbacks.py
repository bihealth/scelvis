"""Registeration of Dash app callbacks.

Each ``register_*`` function sets up callbacks for one aspect of the app.  This module does not
any component building itself.  Instead, this is done in the module ``.ui``.
"""

import base64
import os.path
import shutil
import tempfile
import uuid

import fs
import dash
import dash_html_components as html
from logzero import logger
from werkzeug.utils import secure_filename

from . import ui, settings, store
from .exceptions import ScelVisException
from .ui import cells, common, genes


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

    @app.callback(
        dash.dependencies.Output("page-content", "children"),
        [dash.dependencies.Input("url", "pathname")],
    )
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

    @app.callback(
        dash.dependencies.Output("page-brand", "children"),
        [dash.dependencies.Input("url", "pathname")],
    )
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
        [dash.dependencies.Output("meta_plot_controls", "children")],
        [
            dash.dependencies.Input("url", "pathname"),
            dash.dependencies.Input("meta_plot_type", "value"),
        ],
    )
    def update_meta_plot_controls(pathname, plot_type):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
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
        return common.render_plot("meta", plot_type)


def register_update_cell_scatter_plot_params(app):
    """Register handlers on scatter plot."""

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
    def get_meta_plot_scatter(pathname, xc, yc, col, sample_size):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        return cells.render_plot_scatter(data, xc, yc, col, sample_size)


def register_update_cell_violin_plot_params(app):
    """Register handlers on violin plot."""

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
    def get_meta_plot_violin(pathname, variables, group, split, sample_size):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        return cells.render_plot_violin(data, variables, group, split, sample_size)


def register_update_cell_bar_chart_params(app):
    """Register handlers on bar chart and its controls."""

    @app.callback(
        dash.dependencies.Output("meta_bar_options", "options"),
        [dash.dependencies.Input("meta_bar_select_split", "value")],
    )
    def toggle_meta_bar_options(split):
        if split is None:
            return [
                {
                    "label": "normalized",
                    "value": "normalized",
                    "title": "plot fractions instead of cell numbers",
                }
            ]
        else:
            return [
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
            ]

    @app.callback(
        [
            dash.dependencies.Output("meta_bar_plot", "figure"),
            dash.dependencies.Output("meta_bar_download", "href"),
            dash.dependencies.Output("meta_bar_download", "hidden"),
        ],
        [
            dash.dependencies.Input("url", "pathname"),
            dash.dependencies.Input("meta_bar_select_group", "value"),
            dash.dependencies.Input("meta_bar_select_split", "value"),
            dash.dependencies.Input("meta_bar_options", "values"),
        ],
    )
    def get_meta_plot_bars(pathname, group, split, options):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        return cells.render_plot_bars(data, group, split, options)


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
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
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
        return common.render_plot("expression", plot_type)


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
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        return ui.genes.render_marker_list(data, values)

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
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
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
    def get_expression_plot_scatter(pathname, xc, yc, genelist, sample_size):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        return ui.genes.render_plot_scatter(data, xc, yc, genelist, sample_size)


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
    def get_expression_plot_violin(pathname, genelist, sample_size, group, split):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        return ui.genes.render_plot_violin(data, pathname, genelist, sample_size, group, split)


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
    def get_expression_plot_dot(pathname, genelist, group, split):
        _, kwargs = get_route(pathname)
        data = store.load_data(kwargs.get("dataset"))
        return ui.genes.render_plot_dot(data, pathname, genelist, group, split)


def register_file_upload(app):
    """Register callbacks for the file upload"""

    logger.info("-> file_upload")

    # TODO: move main handler into "upload" module?
    # TODO: upload still assumes about.md file; broken
    @app.callback(
        dash.dependencies.Output("url", "pathname"),
        [dash.dependencies.Input("file-upload", "contents")],
        [dash.dependencies.State("file-upload", "filename")],
    )
    def file_uploaded(contents, filename):
        logger.info("Handling upload of %s", filename)
        filename = secure_filename(filename)
        _content_type, content_string = contents.split(",")
        data_uuid = str(uuid.uuid4())
        logger.info("Data will have UUID %s", data_uuid)
        with tempfile.TemporaryDirectory() as tmpdir:
            logger.info("Using temporary directory %s", tmpdir)
            # Create UUID directory to be moved later.
            uuid_path = os.path.join(tmpdir, data_uuid)
            os.makedirs(uuid_path)
            # Decode base64 string and write out to temporary file.
            filepath = os.path.join(tmpdir, filename)
            with open(filepath, "wb") as tmpf:
                tmpf.write(base64.b64decode(content_string))
            # Try to open archive with PyFilesystem2
            if filename.endswith(".zip"):
                protocol = "zip"
            elif filename.endswith(".tar") or filename.endswith(".tar.gz"):
                protocol = "tar"
            else:
                raise ScelVisException(
                    "Does not have one of the valid extensions .zip, .tar, or .tar.gz: %s"
                    % filename
                )
            # Extract the files (up to a maximal size) from archive.
            with fs.open_fs("%s://%s" % (protocol, filepath)) as archive_fs:
                for path in archive_fs.walk(filter=["about.md"]):
                    if path.files:
                        base_path = path.path
                        break
                else:
                    raise ScelVisException("Could not find about.md in %s" % filename)
                # Extract about.md and data.h5ad.
                for fname, max_size in (
                    ("about.md", settings.MAX_UPLOAD_TEXT_SIZE),
                    ("data.h5ad", settings.MAX_UPLOAD_DATA_SIZE),
                ):
                    path_inf = os.path.join(base_path, fname)
                    path_outf = os.path.join(uuid_path, fname)
                    logger.info(
                        "Extracting %s from %s://%s to %s", path_inf, protocol, filepath, path_outf
                    )
                    with archive_fs.open(path_inf, "rb") as inf:
                        with open(path_outf, "wb") as outf:
                            shutil.copyfileobj(inf, outf, length=max_size)
            # Move into upload directory.
            dest_path = os.path.join(settings.UPLOAD_DIR, data_uuid)
            shutil.move(uuid_path, dest_path)
            # Redirect to view the data set.
            return "/dash/viz/%s" % data_uuid
