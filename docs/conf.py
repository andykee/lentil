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

# Set by the docs workflow to the directory being published ('latest' or a
# release number like '0.8.10') so the version switcher selects the right entry
docs_version = os.environ.get('DOCS_VERSION', 'latest')

# -- General configuration ---------------------------------------------------

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              'sphinx.ext.mathjax',
              'sphinx.ext.napoleon',
              'sphinx.ext.viewcode',
              'sphinx_copybutton',
              'sphinx_design',
              'sphinx_favicon',
              'matplotlib.sphinxext.plot_directive'
              ]
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

html_theme = 'pydata_sphinx_theme'
html_theme_options = {
    'show_prev_next': False,
    'github_url': 'https://github.com/andykee/lentil',
    "logo": {
        "image_light": "_static/logo-light.svg",
        "image_dark": "_static/logo-dark.svg",
    },
    "collapse_navigation": True,
    "navbar_persistent": ["search-button"],
    "navbar_end": ["version-switcher", "theme-switcher", "navbar-icon-links"],
    "article_header_start": [],
    "secondary_sidebar_items": ["page-toc"],
    "switcher": {
       "json_url": "https://andykee.github.io/lentil/versions.json",
       "version_match": docs_version,
    },
    "pygments_light_style": "tango",
    "pygments_dark_style": "nord",
}

# html_additional_pages = {
#     'index': 'index.html'
# }

html_sidebars = {
    '**': ['sidebar-nav-root']
}

html_static_path = ['_static']
html_show_sphinx = False
html_show_sourcelink = False
html_scaled_image_link = False

html_css_files = ['css/lentil.css']

favicons = [
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

# if true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = True

add_function_parentheses = True

autodoc_default_options = {
    'member-order': 'alphabetical',
    'exclude-members': '__new__, __init__, __weakref__, __dict__, __module__',
    'undoc-members': False,
}

autosummary_generate = True

# Some class/function names differ only in their capitalisation. Since these
# resolve to autosummary filenames that also differ only in capitalisation,
# this causes problems when the documentation is built on an OS/filesystem
# that is enforcing case-insensitive semantics.  This setting defines some
# custom names to prevent the clash from happening.
autosummary_filename_map = {
    "lentil.ptype.PType": "lentil.ptype.PType_cls.rst"
}

# Temporary workaround for sphinx #13257
# autodoc: signature of __new__ overwrites signature of __init__
# https://github.com/sphinx-doc/sphinx/issues/13257#issuecomment-4182037875
import inspect


def _fix_class_signature(app, what, name, obj, options, signature, return_annotation):
    if what == 'class' and '__init__' in obj.__dict__:
        sig = inspect.signature(obj.__init__)
        params = [p for p in sig.parameters.values() if p.name != 'self']
        return (str(inspect.Signature(params)), None)


# Hide the secondary sidebar (page TOC) on generated API pages. This can't be
# done with a secondary_sidebar_items dict without pydata-sphinx-theme warning
# that each page matches two wildcard patterns. Runs after the theme's own
# html-page-context handler (priority 500), which sets secondary_sidebar_items.
def _remove_api_secondary_sidebar(app, pagename, templatename, context, doctree):
    if pagename.startswith('api/generated/'):
        context['secondary_sidebar_items'] = []


def setup(app):
    app.connect('autodoc-process-signature', _fix_class_signature)
    app.connect('html-page-context', _remove_api_secondary_sidebar, priority=501)


rst_prolog = """
.. currentmodule:: lentil

.. |Pupil| replace:: :class:`Pupil`
.. |Image| replace:: :class:`Image`
.. |Detector| replace:: :class:`Detector`
.. |Plane| replace:: :class:`Plane`
.. |Wavefront| replace:: :class:`Wavefront`
.. |Tilt| replace:: :class:`Tilt`
.. |Field| replace:: :class:`~lentil.field.Field`

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
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, PowerNorm, AsinhNorm
import numpy as np
np.random.seed(12345)

# set up the OPD colormap to display NaNs as light gray
opd_cmap = matplotlib.colormaps.get_cmap('RdBu_r')
opd_cmap.set_bad(color='#dddddd')

psf_cmap = matplotlib.colormaps.get_cmap('turbo')
"""
