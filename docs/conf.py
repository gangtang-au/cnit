# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# -- Path setup --------------------------------------------------------------
# Add the 'src' directory to sys.path, so Sphinx can find your code modules!
sys.path.insert(0, os.path.abspath('../src'))

# -- Project information -----------------------------------------------------
project = 'CNit'
copyright = '2025, Gang Tang'
author = 'Gang Tang'
release = '2.1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",         # Auto API docs from docstrings
    "sphinx.ext.napoleon",        # Support for NumPy/Google docstrings
    "sphinx_autodoc_typehints",   # Show type hints in docs
    "sphinxcontrib.bibtex",       # BibTeX citations
]

# Templates and static paths
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Prevent ugly default values (like <Quantity(...)>) from showing up
autodoc_class_signature = "separated"

# Autodoc default options
autodoc_default_options = {
    "members": True,
    "exclude-members": "__init__",
    "inherited-members": True,
    "undoc-members": True,
    "private-members": False,
    "show-inheritance": True,
    "member-order": "bysource",
}

# Napoleon extension settings (for NumPy/Google docstrings)
napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_use_rtype = False
napoleon_numpy_returns_no_rtype = True

# Type hints in API docs
typehints_fully_qualified = False
typehints_document_rtype = True
typehints_use_rtype = False

# BibTeX file(s) for references
bibtex_bibfiles = ["bibliography.bib"]

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'         # Use the popular ReadTheDocs theme
html_static_path = ['_static']
html_css_files = ['custom.css']          # Custom CSS file
html_js_files = ['custom.js']

# -- Misc configuration ------------------------------------------------------
# Make the command |CO2| automatically subscript everywhere
rst_epilog = """
.. |CO2| replace:: CO\ :sub:`2`
"""

# Optional: Left-align math equations for MathJax
mathjax3_config = {"chtml": {"displayAlign": "center"}}
