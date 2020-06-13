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
              'sphinx.ext.viewcode']
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

html_theme = 'candle'
html_theme_options = {
    'prev_next_buttons_location': None,
    'collapse_navigation': False,
    'display_version': True,
    'logo_only': True,
    'navbar_links': [
        ('Get Started', 'https://lentil.readthedocs.io/en/latest/index.html#get-started'),
        ('User Guide', 'https://lentil.readthedocs.io/en/latest/index.html#user-guide'),
        ('API', 'https://lentil.readthedocs.io/en/latest/index.html#api'),
        ('Resources', 'https://lentil.readthedocs.io/en/latest/index.html#resources'),
        ('Github', 'https://github.com/andykee/lentil')
    ]
}
html_logo = '_static/img/lentil.png'

html_static_path = ['_static']
html_show_sphinx = False
html_show_sourcelink = False
html_scaled_image_link = False

html_js_files = ['js/copybutton.js']

pygments_style = 'default'

# if true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = True

autodoc_default_options = {
    'member-order': 'alphabetical',
    'exclude-members': '__init__, __weakref__, __dict__, __module__'
}

