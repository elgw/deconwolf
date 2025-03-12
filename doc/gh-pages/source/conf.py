# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'deconwolf'
copyright = '2023-2024, Erik Wernersson'
author = 'Erik Wernersson'
release = '0.4.4'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []

# At some point it would make sense to switch to sphinx also for
# man pages. This is here just to test.
# build with `make man`
# the output will go to build/man
man_pages = [ ("man/dw-nuclei", "dw-nuclei", "", "", 1)]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
# https://github.com/readthedocs/sphinx_rtd_theme
html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
