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

today_fmt = '%B %d, %Y'

# -- General configuration ---------------------------------------------------

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              'sphinx.ext.mathjax',
              'sphinx.ext.napoleon',
              'sphinx.ext.viewcode',
              'sphinx_design',
              'matplotlib.sphinxext.plot_directive']
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
   "pygment_light_style": "tango",
   "pygment_dark_style": "nord",
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

# if true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = True

add_function_parentheses = True

autodoc_default_options = {
    'member-order': 'alphabetical',
    'exclude-members': '__init__, __weakref__, __dict__, __module__',
}

autosummary_generate = True

rst_prolog = """
.. currentmodule:: lentil

.. |Pupil| replace:: :class:`Pupil`
.. |Image| replace:: :class:`Image`
.. |Detector| replace:: :class:`Detector`
.. |Plane| replace:: :class:`Plane`
.. |Wavefront| replace:: :class:`Wavefront`
.. |Tilt| replace:: :class:`Tilt`

"""

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
import lentil
import matplotlib.pyplot as plt
import numpy as np
np.random.seed(12345)
"""

