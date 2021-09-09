.. _apt.internal:

*********
Internals
*********
.. currentmodule:: lentil

.. warning::

    The ``lentil.fourier`` and ``lentil.util`` top-level modules are intended 
    for internal use. Stable functionality is the goal, but is not explicitly
    guaranteed.

Fourier transforms
==================
.. autosummary::
    :toctree: ../generated/

    lentil.fourier.dft2
    lentil.fourier.idft2
    lentil.fourier.czt2
    lentil.fourier.iczt2

Utilities
=========

Shapes
------
.. autosummary::
    :toctree: ../generated/

    lentil.util.mesh
    lentil.util.gaussian2d

Array manipulation
------------------
.. autosummary::
    :toctree: ../generated/

    lentil.util.boundary_slice
    lentil.util.slice_offset

Performance
-----------
.. autosummary::
    :toctree: ../generated/

    lentil.util.expc

Sparse matrix tools
-------------------
.. autosummary::
    :toctree: ../generated/
    
    lentil.util.v2m
    lentil.util.m2v
    lentil.util.make_mask
    lentil.util.make_index
