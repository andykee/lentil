.. _api.internals:

*********
Internals
*********

.. currentmodule:: lentil

Planes
------
.. autosummary::
    :toctree: generated/

    plane._PlaneBase
    plane._TiltBase

Plane type
----------
.. autosummary::
    :toctree: generated/

    ptype.ptype
    ptype.PType

Field
-----
.. autosummary::
    :toctree: generated/

    field.Field
    field.boundary
    field.insert
    field.merge
    field.overlap
    field.reduce

Extent
------
.. autosummary::
    :toctree: generated/

    extent.array_extent
    extent.array_center
    extent.intersect
    extent.intersection_extent
    extent.intersection_shape
    extent.intersection_slices
    extent.intersection_shift

Fourier transforms
------------------
.. autosummary::
    :toctree: generated/

    fourier.dft2

Helper functions
----------------
.. autosummary::
    :toctree: generated/

    helper.mesh
    helper.boundary_slice
    helper.slice_offset
    helper.gaussian2d
    helper.get_rng
