# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import pathlib

# sys.path.insert(0, os.path.abspath('../..'))

CURRENT_PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))
SCORPIO_PATH = CURRENT_PATH.parent.parent


sys.path.insert(0, str(SCORPIO_PATH))


# -- Project information -----------------------------------------------------

project = "SCORPIO"
copyright = "2020, Jose Benavides"
author = "Jose Benavides"

# The full version, including alpha/beta/rc tags
with open(SCORPIO_PATH / "scorpio.py") as fp:
    for line in fp.readlines():
        if line.startswith("__version__ = "):
            release = line.split("=", 1)[-1].replace('"', '').strip()
            break


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    'nbsphinx'
]


# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = [".rst", ".md"]

# The master toctree document.
# master_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
html_theme = 'bootstrap-astropy'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]


# Custom

html_logo = "_static/logo.png"


html_theme_options = {
    'logotext1': 'Scorpio',  # white,  semi-bold
    'logotext2': '',  # orange, light
    'logotext3': ':docs',   # white,  light
    'astropy_project_menubar': False
    }



import m2r


with open(SCORPIO_PATH / "README.md") as fp:
    md = fp.read()


index = f"""
..
   Automatic created file. Don't edit

{m2r.convert(md)}

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   tutorial.ipynb
   api.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

"""

with open(CURRENT_PATH / "index.rst", "w") as fp:
    fp.write(index)
