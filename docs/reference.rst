.. _api:

.. currentmodule:: lentil

*************
API reference
*************

This page gives an overview of all public functions, modules, and objects included
in Lentil. All classes and functions exposed in the ``lentil.*`` namespace are public.

Some subpackages are public including ``lentil.radiometry`` and
``lentil.detector``.

.. _api.planes:

Planes
======
.. autosummary::
    :toctree: generated/

    lentil.Plane
    lentil.Pupil
    lentil.Image
    lentil.Detector
    lentil.Grism
    lentil.Tilt
    lentil.Rotate
    lentil.Flip

.. _api.wavefront:

Wavefront
=========

.. autosummary::
    :toctree: generated/

    lentil.Wavefront

.. _api.propagation:

Numerical diffraction propagation
=================================

.. autosummary::
    :toctree: generated/

    lentil.propagate

.. _api.zernike:

Zernike polynomials
===================
.. autosummary::
    :toctree: generated/

    lentil.zernike
    lentil.zernike_basis
    lentil.zernike_compose
    lentil.zernike_coordinates
    lentil.zernike_fit
    lentil.zernike_remove

.. _api.wfe:

Wavefront errors
================
.. autosummary::
    :toctree: generated/

    lentil.power_spectrum
    lentil.translation_defocus

.. _api.imaging:

Imaging artifacts
=================
.. autosummary::
    :toctree: generated/

    lentil.jitter
    lentil.smear

.. _api.detector:

Detector effects
================

Charge collection
-----------------
.. autosummary::
    :toctree: generated/

    lentil.detector.collect_charge
    lentil.detector.collect_charge_bayer

Pixel effects
-------------
.. autosummary::
    :toctree: generated/

    lentil.detector.pixel
    lentil.detector.pixelate

Noise
-----
.. autosummary::
    :toctree: generated/

    lentil.detector.shot_noise
    lentil.detector.read_noise
    lentil.detector.charge_diffusion
    lentil.detector.dark_current
    lentil.detector.rule07_dark_current

Readout
-------
.. autosummary::
    :toctree: generated/

    lentil.detector.adc

Cosmic rays
-----------
.. autosummary::
    :toctree: generated/

    lentil.detector.cosmic_rays

.. _api.radiometry:

Radiometry
==========

.. autosummary::
    :toctree: generated/

    lentil.radiometry.Spectrum
    lentil.radiometry.Blackbody
    lentil.radiometry.Material
    lentil.radiometry.planck_radiance
    lentil.radiometry.planck_exitance
    lentil.radiometry.vegaflux
    lentil.radiometry.path_transmission
    lentil.radiometry.path_emission

.. _api.util:

Utilities
=========

Shapes
------
.. autosummary::
    :toctree: generated/

    lentil.circle
    lentil.circlemask
    lentil.hexagon
    lentil.slit

Array manipulation
------------------
.. autosummary::
    :toctree: generated/

    lentil.boundary
    lentil.centroid
    lentil.pad
    lentil.window
    lentil.rebin
    lentil.rescale

Miscellaneous
-------------
.. autosummary::
    :toctree: generated/

    lentil.pixelscale_nyquist
    lentil.normalize_power

.. _apt.internal:

Internals
=========

.. warning::

    The ``lentil.field``, ``lentil.fourier``, ``lentil.helper``, and ``lentil.util``
    top-level modules are intended for internal use. Stable functionality in these
    modules is not guaranteed.

Field
-----
.. autosummary::
    :toctree: generated/

    lentil.field.Field
    lentil.field.NDField
    lentil.field.extent
    lentil.field.overlap
    lentil.field.boundary
    lentil.field.insert
    lentil.field.merge
    lentil.field.multiply
    lentil.field.reduce

Fourier transforms
------------------
.. autosummary::
    :toctree: generated/

    lentil.fourier.dft2
    lentil.fourier.idft2
    lentil.fourier.czt2
    lentil.fourier.iczt2

Helper functions
----------------
.. autosummary::
    :toctree: generated/


    lentil.helper.slice_offset

Utilities
---------

.. autosummary::
    :toctree: generated/

    lentil.util.mesh
    lentil.util.gaussian2d
    lentil.util.boundary_slice
    lentil.util.expc
    lentil.util.v2m
    lentil.util.m2v
    lentil.util.make_mask
    lentil.util.make_index
