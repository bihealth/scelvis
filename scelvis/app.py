"""Setup the SCelVis Dash application.

Note that before importing this module, you will have to configure the settings in ``.settings``.  See this module's
docstring for more information.
"""

import os.path

import dash

from .__init__ import __version__
from .layout import layout
from . import cache, callbacks

#: Path to assets.
ASSETS_FOLDER = os.path.join(os.path.dirname(__file__), "assets")


#: The Dash application to run.
app = dash.Dash(
    __name__,
    assets_folder=ASSETS_FOLDER,
    # external_stylesheets=["/%s/bootstrap.min.css" % ASSETS_ROUTE, "/%s/scelvis.css" % ASSETS_ROUTE],
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
app.layout = layout()

# Register the callbacks with the app.
callbacks.register_page_content(app)
callbacks.register_page_brand(app)

callbacks.register_select_cell_plot_type(app)
callbacks.register_update_cell_scatter_plot_params(app)
callbacks.register_update_cell_violin_plot_params(app)
callbacks.register_update_cell_bar_chart_params(app)

callbacks.register_select_gene_plot_type(app)
callbacks.register_select_gene_marker_list(app)
callbacks.register_select_gene_scatter_plot(app)
callbacks.register_select_gene_violin_plot(app)
callbacks.register_select_gene_dot_plot(app)
