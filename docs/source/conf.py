# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'polyase'
copyright = '2025, Nadja Nolte'
author = 'Nadja Nolte'
release = '0.1.0'

import os
import sys
sys.path.insert(0, os.path.abspath('../../polyase/'))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.coverage', 'sphinx.ext.napoleon', 'nbsphinx' ]

templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']




def setup(app):
    app.add_css_file('my_theme.css')
# Auto-generate API documentation


# Configure nbsphinx for RequireJS (critical for Plotly)
nbsphinx_requirejs_path = ''  # Disable the default RequireJS
nbsphinx_requirejs_options = {
    "paths": {
        "plotly": "https://cdn.plot.ly/plotly-2.27.0.min"
    }
}

# Add plotly script (keep this)
html_js_files = [
    'https://cdn.plot.ly/plotly-2.27.0.min.js',
]

# Ensure notebooks are executed with outputs
nbsphinx_execute = 'never'  # Use 'always' if you want RTD to execute them

# If using sphinx_rtd_theme, you might need this
html_theme = 'sphinx_rtd_theme'

nbsphinx_execute = 'never'

