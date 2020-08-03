"""Setup the SCelVis Dash application.

When importing this module, the app is built and configurd.  Thus, it is important that before
this module is imported, the values in ``.settings`` must already have been setup.
"""

import os.path
import os
import tempfile
import tarfile
import zipfile

import dash
import flask
import uuid
import anndata
from flask import request, helpers, send_from_directory
from logzero import logger
from werkzeug.utils import secure_filename


from . import convert, cache, callbacks, settings
from .__init__ import __version__
from .ui.main import build_layout

#: Path to assets.
ASSETS_FOLDER = os.path.join(os.path.dirname(__file__), "assets")

#: The Flask application to use.
app_flask = flask.Flask(__name__)

# Setup temporary upload folder
app_flask.config["UPLOAD_FOLDER"] = settings.TEMP_DIR
# Setup maximal file upload size
app_flask.config["MAX_CONTENT_LENGTH"] = settings.MAX_UPLOAD_DATA_SIZE
# Setup URL prefix for Flask.
app_flask.config["APPLICATION_ROOT"] = "%s/" % settings.PUBLIC_URL_PREFIX

#: The Dash application to run.
app = dash.Dash(
    __name__,
    # Use our specific Flask app
    server=app_flask,
    # All files from "assets/" will be served as "/assets/*"
    assets_folder=ASSETS_FOLDER,
    # The visualization will be served below "/dash"
    routes_pathname_prefix="/dash/",
    requests_pathname_prefix="%s/dash/" % settings.PUBLIC_URL_PREFIX,
)

# Setup the cache.
cache.setup_cache(app)

# preload data if required
if settings.CACHE_PRELOAD_DATA:
    from . import store

    logger.info("preloading data ...")
    store.load_all_metadata()
    logger.info("preloading done.")

# Set app title
app.title = "SCelVis v%s" % __version__

# Serve assets locally
app.css.config.serve_locally = True
app.scripts.config.serve_locally = True

app.config.suppress_callback_exceptions = True

# Setup the application's main layout.
app.layout = build_layout()

# Register the callbacks with the app.
#
# Callbacks for title and main content.
callbacks.register_page_content(app)
callbacks.register_page_brand(app)

# Callbacks for the "cells" tab pane.
callbacks.register_select_cell_plot_type(app)
callbacks.register_update_cell_scatter_plot_params(app)
callbacks.register_toggle_select_cells_controls(app)
callbacks.register_update_select_cells_selected(app)
callbacks.register_activate_select_cells_buttons(app)
callbacks.register_update_cell_violin_plot_params(app)
callbacks.register_update_cell_box_plot_params(app)
callbacks.register_update_cell_bar_chart_params(app)
callbacks.register_toggle_filter_cells_controls(app, "meta")
callbacks.register_update_filter_cells_controls(app)

# Cellbacks for the "genes" tab pane.
callbacks.register_select_gene_plot_type(app)
callbacks.register_select_gene_list(app)
callbacks.register_select_gene_scatter_plot(app)
callbacks.register_select_gene_violin_plot(app)
callbacks.register_select_gene_box_plot(app)
callbacks.register_select_gene_dot_plot(app)
callbacks.register_toggle_filter_cells_controls(app, "expression")

# callbacks for the filter cells div (on "cells" tab pane, but required for both)
callbacks.register_update_filter_cells_filters(app)
callbacks.register_activate_filter_cells_reset(app)
callbacks.register_run_differential_expression(app)

# Callbacks for the file upload.
callbacks.register_file_upload(app)


# Add redirection for root.
@app_flask.route("/")
def redirect_root():
    return flask.redirect("%s/dash/" % settings.PUBLIC_URL_PREFIX)


# add redirection for custom static folder
if settings.CUSTOM_STATIC_FOLDER:

    @app.server.route(os.path.join("/static", "<filename>"))
    def custom_static_route(filename):
        return send_from_directory(settings.CUSTOM_STATIC_FOLDER, filename)


# Mount upload site.
@app_flask.route("/upload/", methods=("GET", "POST"))
def upload_route():
    """Perform file upload."""

    import faulthandler

    faulthandler.enable()

    if request.method == "POST":
        if (
            "file" not in request.files
            or not request.files["file"].filename
            or not request.files["file"].filename.endswith(".h5ad")
        ):
            return """
            <!doctype html>
            <p>no .h5ad file provided!</p>
            <p>
            <a href="%(application_root)s/dash/upload" target="_PARENT">try again</a>
            </p>""" % {
                "application_root": settings.PUBLIC_URL_PREFIX
            }
        file = request.files["file"]
        data_uuid = str(uuid.uuid4())
        logger.info("Data will have UUID %s", data_uuid)
        # Decode base64 string and write out to final file.
        filepath = os.path.join(settings.UPLOAD_DIR, "%s.h5ad" % data_uuid)
        logger.info("Writing to %s", filepath)
        file.save(filepath)
        try:
            anndata.read_h5ad(filepath)
            return """
            <!doctype html>
            <p>upload successful!</p>
            <p>
            <a href="%(application_root)s/dash/viz/%(data_uuid)s" target="_PARENT">view uploaded dataset</a>
            </p>""" % {
                "application_root": settings.PUBLIC_URL_PREFIX,
                "data_uuid": data_uuid,
            }
        except OSError as err:
            return """
            <!doctype html>
            <p>%(err)s</p>
            <p>
            <a href="%(application_root)s/dash/upload" target="_PARENT">try again</a>
            </p>""" % {
                "application_root": settings.PUBLIC_URL_PREFIX,
                "err": err,
            }
    else:
        return """
        <!doctype html>
        <form method=post enctype=multipart/form-data>
        <input type=file name=file>
        <input type=submit value=Upload>
        </form>
        """


# enable download of conversion results
@app_flask.route("/download/<string:data_uuid>", methods=("GET",))
def download_route(data_uuid):
    """Download converted file."""
    try:
        data_uuid = str(uuid.UUID(data_uuid))
    except ValueError:
        return """
        <!doctype html>
        <p>invalid identifier</p>
        <p>
        """
    out_file = os.path.join(settings.UPLOAD_DIR, data_uuid + ".h5ad")
    return helpers.send_file(out_file, mimetype="application/binary", as_attachment=True)


# Mount conversion site.
@app_flask.route("/convert/", methods=("GET", "POST"))
def convert_route():
    """Perform conversion file upload."""

    def find(name, path):
        for root, _dirs, files in os.walk(path):
            if name in files:
                return os.path.join(root, name)

    if request.method == "POST":
        if "file" not in request.files or not request.files["file"].filename:
            return """
            <!doctype html>
            <p>no valid file provided!</p>
            <p>
            <a href="%(application_root)s/dash/convert" target="_PARENT">try again</a>
            </p>""" % {
                "application_root": settings.PUBLIC_URL_PREFIX
            }
        file = request.files["file"]
        filename = secure_filename(file.filename)
        filepath = os.path.join(app_flask.config["UPLOAD_FOLDER"], filename)
        file.save(filepath)
        with tempfile.TemporaryDirectory() as tmpdir:
            if filepath.endswith(".zip"):
                logger.info("Extracting ZIP file %s to %s", filepath, tmpdir)
                with zipfile.ZipFile(filepath) as zipf:
                    zipf.extractall(tmpdir)
            elif filepath.endswith(".tar") or filepath.endswith(".tar.gz"):
                logger.info("Extracting TAR file %s to %s", filepath, tmpdir)
                with tarfile.open(filepath) as tarf:
                    tarf.extractall(tmpdir)
            else:
                return """
                <!doctype html>
                <p>file does not have .zip, .tar or .tar.gz extension!</p>
                <p>
                <a href="%(application_root)s/dash/convert" target="_PARENT">try again</a>
                </p>""" % {
                    "application_root": settings.PUBLIC_URL_PREFIX
                }
            cellranger_needle = "filtered_feature_bc_matrix.h5"
            logger.info("Looking for %s file", cellranger_needle)
            needle_path = find(cellranger_needle, tmpdir)
            if needle_path is None:
                text_needle = "coords.tsv"
                logger.info("Looking for %s file", text_needle)
                needle_path = find(text_needle, tmpdir)
                if needle_path is None:
                    loom_needle = "data.loom"
                    logger.info("Looking for %s file", loom_needle)
                    needle_path = find(loom_needle, tmpdir)
                    format_ = "loom"
                else:
                    format_ = "text"
            else:
                format_ = "cell-ranger"
            input_dir = os.path.dirname(needle_path)
            logger.info("Starting conversion...")
            with tempfile.TemporaryDirectory() as tmpdir2:
                logger.info("Writing about.md")
                about_md = os.path.join(tmpdir2, "about.md")
                with open(about_md, "wt") as aboutf:
                    title = request.form.get("title")
                    short_title = request.form.get("short_title")
                    description = request.form.get("description")
                    if title or short_title:
                        print("----", file=aboutf)
                        print("title: %s" % (title or "Untitled"), file=aboutf)
                        if short_title:
                            print("short_title: %s" % short_title, file=aboutf)
                        print("----", file=aboutf)
                    print(description or "No description", file=aboutf)
                data_uuid = str(uuid.uuid4())
                logger.info("Data will have UUID %s", data_uuid)
                out_file = os.path.join(settings.UPLOAD_DIR, data_uuid + ".h5ad")
                logger.info("Performing conversion (%s)", out_file)
                try:
                    convert.run(
                        convert.Config(
                            indir=input_dir, about_md=about_md, out_file=out_file, format=format_
                        )
                    )
                    return """
                    <!doctype html>
                    <p>conversion successful!</p>
                    <p>
                    <a href="%(application_root)s/dash/viz/%(data_uuid)s" target="_PARENT">view</a>
                    or
                    <a href="%(application_root)s/download/%(data_uuid)s" target="_PARENT">download</a>
                    converted dataset
                    </p>
                    """ % {
                        "application_root": settings.PUBLIC_URL_PREFIX,
                        "data_uuid": data_uuid,
                    }
                except Exception as err:
                    return """
                    <!doctype html>
                    <p>conversion failed (%(err)s)</p>
                    <p>
                    <a href="%(application_root)s/dash/convert" target="_PARENT">try again</a>
                    </p>""" % {
                        "application_root": settings.PUBLIC_URL_PREFIX,
                        "err": err,
                    }

    else:
        return """
        <!doctype html>
        <form method=post enctype=multipart/form-data>
        <label>Title</label> <input type=text name=title><br>
        <label>Short Title</label> <input type=text name=short_title><br>
        <label>Description</label><br>
        <textarea cols=40 rows=5 name=description></textarea><br>
        <input type=file name=file>
        <input type=submit value=Upload>
        </form>
        """
