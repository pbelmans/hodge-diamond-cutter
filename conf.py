# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

# from sage.env import SAGE_DOC_SRC, SAGE_DOC, SAGE_SRC
#
# try:
#    import sage.all
# except ImportError:
#    raise RuntimeError("to build the documentation you need to be inside a Sage shell (run first the command 'sage -sh' in a shell")

import os
import sys
from sage.env import SAGE_DOC_SRC
import sage.all

sys.path.insert(0, os.path.abspath("."))

# -- Project information -----------------------------------------------------

project = "Hodge diamond cutter"
copyright = "2021, Pieter Belmans"
author = "Pieter Belmans"
# The full version, including alpha/beta/rc tags
release = "v1.2"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
]

mathjax3_config = {
    "tex": {
        "inlineMath": [["$", "$"]],
        "displayMath": [["$$", "$$"], ["\\[", "\\]"]],
    },
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"
html_logo = "logo.png"
html_favicon = "favicon.ico"

html_js_files = [
    (
        "https://plausible.io/js/script.js",
        {"data-domain": "cutter.ncag.info", "defer": "defer"},
    ),
]
