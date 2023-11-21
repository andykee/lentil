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
copyright = f'2017-{datetime.now().year} California Institute of Technology'

with open(os.path.normpath(os.path.join(path, '..', 'lentil', '__init__.py'))) as f:
    version = release = re.search("__version__ = '(.*?)'", f.read()).group(1)


# -- General configuration ---------------------------------------------------

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              'sphinx.ext.mathjax',
              'sphinx.ext.napoleon',
              'sphinx.ext.viewcode',
              'sphinx_remove_toctrees',
              'sphinx_copybutton',
              'sphinx_design',
              'matplotlib.sphinxext.plot_directive',
              'numpydoc']
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'docs'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

html_theme = 'pydata_sphinx_theme'
html_theme_options = {
    'show_prev_next': False,
    'github_url': 'https://github.com/andykee/lentil',
    "logo": {
        "link": "docs",
        "image_light": "_static/logo-light.svg",
        "image_dark": "_static/logo-dark.svg",
   },
   "collapse_navigation": True,
   "navbar_persistent": ["search-button"],
   "favicons": [
      {
         "rel": "icon",
         "sizes": "16x16",
         "href": "favicon/favicon-16x16.png",
      },
      {
         "rel": "icon",
         "sizes": "32x32",
         "href": "favicon/favicon-32x32.png",
      },
      {
         "rel": "icon",
         "sizes": "48x48",
         "href": "favicon/favicon-48x48.png",
      },
      {
         "rel": "apple-touch-icon",
         "sizes": "180x180",
         "href": "favicon/apple-touch-icon-180x180.png",
         "color": "#000000",
      },
   ]
}

html_additional_pages = {
    'index': 'index.html'
}

html_sidebars = {
    'index': []
}

html_static_path = ['_static']
html_show_sphinx = False
html_show_sourcelink = False
html_scaled_image_link = False

html_css_files = ['css/lentil.css']

pygments_style = 'default'

# if true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = True

autodoc_default_options = {
    'member-order': 'bysource',
    'undoc-members': None
#    'exclude-members': '__init__, __weakref__, __dict__, __module__'
}

autosummary_generate = True

#remove_from_toctrees = ["generated/*"]

# -- Plot config -------------------------------------------------------------
dpi = 144

plot_rcparams = {}  # noqa
plot_rcparams['font.size'] = 12*72/dpi  # 12 pt
plot_rcparams['axes.titlesize'] = 14*72/dpi  # 14 pt
plot_rcparams['axes.labelsize'] = 12*72/dpi  # 12 pt
plot_rcparams['axes.linewidth'] = 0.5
plot_rcparams['lines.linewidth'] = 1
plot_rcparams['lines.markersize'] = 2
plot_rcparams['xtick.major.width'] = 0.5
plot_rcparams['xtick.major.size'] = 2
plot_rcparams['ytick.major.width'] = 0.5
plot_rcparams['ytick.major.size'] = 2
plot_rcparams['grid.linewidth'] = 0.5
plot_rcparams['xtick.labelsize'] = 12*72/dpi  # 12 pt
plot_rcparams['ytick.labelsize'] = 12*72/dpi  # 12 pt
plot_rcparams['legend.fontsize'] = 12*72/dpi  # 12 pt
plot_rcparams['figure.figsize'] = (2.5, 2.5)
plot_rcparams['figure.subplot.wspace'] = 0.2
plot_rcparams['figure.subplot.hspace'] = 0.2
plot_rcparams['savefig.bbox'] = 'tight'
plot_rcparams['savefig.transparent'] = True

plot_apply_rcparams = True
plot_html_show_source_link = False
plot_html_show_formats = False
plot_formats = [('png', dpi*2)]
plot_pre_code = """
import numpy as np
np.random.seed(12345)
"""




#def fix_attributes(app, pagename, templatename, context, doctree):
#   if 'generated' in pagename:
#      context['body'] = context['body'].replace('Variables', 'Attributes')

#def setup(app):
#     app.connect("html-page-context", fix_attributes)