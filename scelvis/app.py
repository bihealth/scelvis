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
from flask import flash, request, redirect, helpers
from logzero import logger
from werkzeug.utils import secure_filename


from . import convert, cache, callbacks, settings
from .__init__ import __version__
from .exceptions import ScelVisException
from .ui.main import build_layout

#: Path to assets.
ASSETS_FOLDER = os.path.join(os.path.dirname(__file__), "assets")

#: The Flask application to use.
app_flask = flask.Flask(__name__)

# Setup temporary upload folder
app_flask.config["UPLOAD_FOLDER"] = settings.TEMP_DIR
# Setup maximal file upload size
app_flask.config["MAX_CONTENT_LENGTH"] = settings.MAX_UPLOAD_SIZE

#: The Dash application to run.
app = dash.Dash(
    __name__,
    # Use our specific Flask app
    server=app_flask,
    # All files from "assets/" will be served as "/assets/*"
    assets_folder=ASSETS_FOLDER,
    # The visualization will be served below "/dash"
    url_base_pathname="/dash/",
)

# Setup the cache.
cache.setup_cache(app)

# Set app title
app.title = "SCelVis v%s" % __version__

# Serve assets locally
app.css.config.serve_locally = True
app.scripts.config.serve_locally = True

# TODO: Better use the approach from this URL:
# - https://community.plot.ly/t/dynamic-controls-and-dynamic-output-components/5519
app.config.supress_callback_exceptions = True

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
callbacks.register_update_cell_violin_plot_params(app)
callbacks.register_update_cell_bar_chart_params(app)

# Cellbacks for the "genes" tab pane.
callbacks.register_select_gene_plot_type(app)
callbacks.register_select_gene_marker_list(app)
callbacks.register_select_gene_scatter_plot(app)
callbacks.register_select_gene_violin_plot(app)
callbacks.register_select_gene_dot_plot(app)

# Callbacks for the file upload.
callbacks.register_file_upload(app)

# Add redirection for root.
@app_flask.route("/")
def redirect_root():
    return flask.redirect("/dash/")


# Mount conversion site.
@app_flask.route("/convert/", methods=("GET", "POST"))
def convert_route():
    """Perform conversion file upload."""

    def find(name, path):
        for root, dirs, files in os.walk(path):
            if name in files:
                return os.path.join(root, name)

    # TODO: error handling, error handling, error handling
    # TODO: prettify the form
    # TODO: protect against "zip bomb"
    if request.method == "POST":
        if "file" not in request.files or not request.files["file"].filename:
            flash("No file uploaded!")
            return redirect(request.url)
        file = request.files["file"]
        filename = secure_filename(file.filename)
        filepath = os.path.join(app_flask.config["UPLOAD_FOLDER"], filename)
        file.save(filepath)
        with tempfile.TemporaryDirectory() as tmpdir:
            if filepath.endswith(".zip"):
                logger.info("Extracting ZIP file %s to %s", filepath, tmpdir)
                with zipfile.ZipFile(filepath) as zipf:
                    zipf.extractall(tmpdir)
            elif filepath.endswith(".tar") or filepath.endwith(".tar.gz"):
                logger.info("Extracting TAR file %s to %s", filepath, tmpdir)
                with tarfile.TarFile(filepath) as tarf:
                    tarf.extractall(tmpdir)
            else:
                raise ScelVisException(
                    "Does not have one of the valid extensions .zip, .tar, or .tar.gz: %s"
                    % filepath
                )
            needle = "filtered_feature_bc_matrix.h5"
            logger.info("Looking for %s file", needle)
            needle_path = find(needle, tmpdir)
            # TODO: properly handle not found
            input_dir = os.path.dirname(needle_path)
            logger.info("Starting conversion...")
            with tempfile.TemporaryDirectory() as tmpdir2:
                output_dir = os.path.join(tmpdir2, filename)
                os.makedirs(output_dir)
                convert.run(convert.Config(indir=input_dir, outdir=output_dir))
                outzip = os.path.join(tmpdir2, "%s.zip" % output_dir)
                logger.info("Writing about.md")
                with open(os.path.join(output_dir, "about.md"), "wt") as aboutf:
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
                logger.info("Creating output ZIP file %s", outzip)
                with zipfile.ZipFile(outzip, "w") as zipf:
                    for member_name in ("about.md", "markers.csv", "data.h5ad"):
                        path_in = os.path.join(output_dir, member_name)
                        path_zip = os.path.join(filename, member_name)
                        logger.info("Adding %s as %s", path_in, path_zip)
                        zipf.write(path_in, path_zip)
                # Send file to the user.
                return helpers.send_file(outzip, mimetype="application/zip", as_attachment=True)
    else:
        return """
            <!doctype html>
            <title>Convert File</title>
            <h1>Upload ZIP or TAR.GZ of CellRanger Output</h1>
            <a href="/dash/">Back to Visualisation</a>
            <form method=post enctype=multipart/form-data>
            <label>Title</label> <input type=text name=title><br>
            <label>Short Title</label> <input type=text name=short_title><br>
            <label>Description</label><br>
            <textarea cols=40 rows=5 name=description></textarea><br>
            <label>CellRanger Output</label> <input type=file name=file>
            <input type=submit value=Upload>
            </form>
            """
