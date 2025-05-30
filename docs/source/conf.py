# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
#from unittest.mock import MagicMock


sys.path.append(os.path.abspath("..")+'/..')
#sys.path.insert(os.path.abspath('../..')) 	
sys.path.insert(0, os.path.abspath("../../"))
 
# Fix matplotlib non import
MOCK_MODULES = ['yaml']
    

# -- Project information -----------------------------------------------------

project = 'Seismic Graph'
copyright = '2025, Rouskin Lab'
author = "Yves Martin des Taillades, Scott Grote, Matthew Allan, Alberic de Lajarte, Casper L'Esperance-Kerckhoff"

# The full version, including alpha/beta/rc tags
release = '0.0.1'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

from recommonmark.parser import CommonMarkParser

source_parsers = {
    '.md': CommonMarkParser,
}

source_suffix = ['.rst', '.md']

extensions = [
    'sphinx.ext.napoleon',
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    # 'sphinxcontrib.blockdiag',
    'sphinx.ext.autosectionlabel',
    'recommonmark',
    'sphinx_panels',
    'sphinx_click',
    ]

# add search bar
html_sidebars = {
    "**": ["search-field.html", "sidebar-nav-bs.html", "sidebar-ethical-ads.html"]
}

# Fontpath for blockdiag (truetype font)
blockdiag_fontpath = '/usr/share/fonts/truetype/ipafont/ipagp.ttf'

# Provide a GitHub API token:
# Pass the SPHINX_GITHUB_CHANGELOG_TOKEN environment variable to your build
# OR
sphinx_github_changelog_token = "..."

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True
napoleon_type_aliases = None
napoleon_attr_annotations = True

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output
html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'

autosectionlabel_prefix_document = True

# Options for build
# fail_on_warning = True
# autodoc_mock_imports = MOCK_MODULES
nitpick_ignore = [('py:class', 'type')]

# Generate the plots for the gallery
sys.path.append(os.path.abspath(""))
from scripts import gallery_generator
gallery_generator.main()
