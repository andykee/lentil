# Configuration file for the Sphinx documentation builder.

import os
import sys
import re

from datetime import datetime

# -- Path setup --------------------------------------------------------------

sys.path.insert(0, os.path.abspath('../'))
path = os.path.abspath(os.path.dirname(__file__))


# -- Project information -----------------------------------------------------

project = 'Lentil'
author = 'Andy Kee'
copyright = f'{datetime.now().year} California Institute of Technology'

with open(os.path.normpath(os.path.join(path, '..', 'lentil', '__init__.py'))) as f:
    version = release = re.search("__version__ = '(.*?)'", f.read()).group(1)


# -- General configuration ---------------------------------------------------

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              'sphinx.ext.mathjax',
              'sphinx.ext.napoleon',
              'sphinx.ext.viewcode',
              'sphinx_remove_toctrees']
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

html_theme = 'pydata_sphinx_theme'
html_theme_options = {
    'show_prev_next': False,
    'google_analytics_id': 'UA-180546240-1',
    'github_url': 'https://github.com/andykee/lentil'
}
html_logo = '_static/img/lentil.png'

html_additional_pages = {
    'index': 'indexcontent.html'
}

html_static_path = ['_static']
html_show_sphinx = False
html_show_sourcelink = False
html_scaled_image_link = False

html_js_files = ['js/copybutton.js']
html_css_files = ['css/lentil.css', 'css/syntax-highlighting.css']

pygments_style = 'default'

# if true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = True

autodoc_default_options = {
    'member-order': 'alphabetical',
    'exclude-members': '__init__, __weakref__, __dict__, __module__'
}

autosummary_generate = True

#remove_from_toctrees = ["generated/*"]
