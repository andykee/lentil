**********************
Optimizing Performance
**********************

.. _caching:

Caching
=======
A number of plane attributes are accessed with each propagation wavelength. This
behavior does does not impact performance unless any of these attributes are computed on
the fly or are otherwise expensive to retrieve. To mitigate any potential performance
impacts, Lentil's :func:`~lentil.propagate` method performs a pre-propagation caching
step by calling each plane's :func:`~lentil.Plane.cache_propagate` method, temporarily
storing a copy of possibly expensive to compute attributes for faster access during the
actual numerical propagation. When the propagation calculations are complete, a
post-propagation cleanup calls each plane's :func:`~lentil.Plane.clear_cache_propagate`
to clear any cached values.

The cached attributes are defined in a list in each Plane's
:attr:`~lentil.Plane.cache_attrs` attribute. This list is user-settable but the only
valid values are `'amplitude'` and `'phase'`. The default behavior is to cache both
amplitude and phase attributes.


.. _performance-image-simulation:

Image simulation
================



.. Speeding up the FFT
.. ===================


.. Multiprocessing
.. ===============

.. Other Performance Tweaks
.. ========================

.. Faster photon to electron (quantum efficiency) calculations
.. -----------------------------------------------------------
