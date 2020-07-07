.. _api-lentil:

******
lentil
******

.. currentmodule:: lentil

Optical Planes
==============
.. autosummary::

    lentil.Plane
    lentil.Pupil
    lentil.Image

Special Planes
==============
.. autosummary::

    lentil.Grism

Diffraction Modeling
====================
.. autosummary::

    lentil.propagate

Coordinate Transformations
==========================
.. autosummary::

    lentil.Tilt
    lentil.Rotate
    lentil.Flip


.. autoclass:: lentil.Plane
    :members:

.. autoclass:: lentil.Pupil
    :members:
    :inherited-members:

.. autoclass:: lentil.Image
    :members: pixelscale, shape, multiply, frame, mask

.. autoclass:: lentil.Grism

.. autofunction:: lentil.propagate

.. autoclass:: lentil.Tilt

.. autoclass:: lentil.Rotate

.. autoclass:: lentil.Flip
