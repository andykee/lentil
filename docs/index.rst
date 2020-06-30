.. title:: Lentil

********************
Lentil Documentation
********************

Lentil is a Python library for modeling the imaging chain of an optical system.
It was originally developed at NASA's Jet Propulsion Lab by the Wavefront Sensing and
Control group (383E) to provide an easy to use framework for simulating point spread
functions of segmented aperture telescopes.

.. note::

    Lentil is still under active development and new features continue to be added.
    Until Lentil reaches version 1.0, the API is not guaranteed to be stable, but
    changes breaking backwards compatibility will be noted.

.. toctree::
    :caption: Get Started
    :name: get_started
    :maxdepth: 1

    installation
    basics

.. toctree::
    :caption: User Guide
    :name: user-guide
    :maxdepth: 1

    user_guide/overview
    user_guide/plane
    user_guide/diffraction
    user_guide/wavefront_error
    user_guide/radiometry
    user_guide/image_sensors
    user_guide/imaging_artifacts
    patterns/index
    user_guide/performance
    user_guide/matlab

.. toctree::
    :caption: API
    :name: api
    :maxdepth: 1

    api/lentil
    api/convolvable
    api/detector
    api/modeltools
    api/radiometry
    api/util
    api/wfe
    api/zernike
    api/private

.. toctree::
    :caption: Resources
    :name: resources
    :maxdepth: 1

    changes
    Issues <https://github.com/andykee/lentil/issues>
    contributing
    License <https://github.com/andykee/lentil/blob/master/LICENSE.rst>
