# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Workaround --------------------------------------------------------------
#
# Workaround to install and execute git-lfs on Read the Docs
#
# cf. https://github.com/readthedocs/readthedocs.org/issues/1846#issuecomment-477184259
import os

if not os.path.exists("./git-lfs"):
    os.system(
        "wget https://github.com/git-lfs/git-lfs/releases/download/v2.7.1/git-lfs-linux-amd64-v2.7.1.tar.gz"
    )
    os.system("tar xvfz git-lfs-linux-amd64-v2.7.1.tar.gz")
    os.system("./git-lfs install")  # make lfs available in current repository
    os.system("./git-lfs fetch")  # download content from remote
    os.system("./git-lfs checkout")  # make local files to have the real content on them
# ----------------------------------------------------------------------------


# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
# sys.path.insert(0, os.path.abspath('.'))

# -- Project information -----------------------------------------------------

project = "SCelVis"
copyright = "2019-2020, BIH Core Unit Bioinformatics"
author = "Benedikt Obermayer, Manuel Holtgrewe"

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = ""
# The full version, including alpha/beta/rc tags.
release = ""

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = []

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build"]

# The master toctree document.
master_doc = "index"


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
