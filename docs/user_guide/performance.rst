**********************
Optimizing Performance
**********************

Selecting sensible propagation parameters
=========================================
Calls to the the :func:`~lentil.propagate` method are almost always the most expensive
part of any model. Even the most efficient and streamlined models will suffer from poor
performance when working with large arrays and many wavelengths. Understanding how each
of the :func:`~lentil.propagate` method's parameters influence accuracy and speed will
allow you to identify appropriate values for your specific needs.

``propagate()``
---------------

.. list-table:: 
   :widths: 20 40 40
   :header-rows: 1

   * - Parameter
     - Usage
     - Performance considerations
   * - ``wave``, ``weight``
     - Arrays representing wavelength and the corresponding weight of each wavelength in
       the propagation.
     - Propagation time scales linearly with the number of wavelengths. 
   * - ``npix``
     - Shape of output plane.
     - The output plane size has minimal impact on propagation performance unless it is
       astronomically large or if ``flatten`` is ``False``. ``npix`` should be set to 
       ensure all data is adequately captured by the output plane.
   * - ``npix_chip``
     - Shape of propagation plane.
     - Propagation time scales quadtatically with ``npix_chip`` but because values are
       typically small there are unlikely to be any large opportunities for increasing
       performance by reducing ``npix_chip``.
   * - ``oversample``
     - Number of times to oversample the output and propagation planes.
     - Propagation time scales quadratically with ``oversample``. For accuracy, ``oversample``
       should be selected to ensure propagations are Nyquist sampled, but there is typically no 
       benefit in selecting larger values.
   * - ``flatten``
     - If ``True``, the cube of wavelength-dependent output planes is flattened into a single
       2D array before being returned. If ``False``, a cube of output planes is returned.
     - 

``radiometry.trim_spectrum()``
------------------------------

.. list-table:: 
   :widths: 20 40 40
   :header-rows: 1

   * - Parameter
     - Usage
     - Performance considerations
   * - ``trim_tol``
     - 
     - 

``radiometry.sample_spectrum()`` and ``radiometry.trim_spectrum()``
-------------------------------------------------------------------

.. list-table:: 
   :widths: 20 40 40
   :header-rows: 1

   * - Parameter
     - Usage
     - Performance considerations
   * - ``wave_sampling``
     - 
     - 

Multiprocessing
===============

.. _performance-image-simulation:


Using appropriately sized planes
================================
Planes should be sized to ensure the smallest spatial features of interest are
adequately sampled. 


Image simulation
================


.. _caching:

Caching
=======

Plane attributes
----------------
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

DFT matrices
------------

Profiling your code
===================
There are several approaches to finding bottlenecks and inefficiencies. To really 
understand what is happening, the code needs to be profiled. The Python standard
library includes several `profilers <https://docs.python.org/3/library/profile.html>`_.
``cProfile`` is simple to use and its profile results files can be visualized using
`snakeviz <https://jiffyclub.github.io/snakeviz/>`_.


.. Faster photon to electron (quantum efficiency) calculations
.. -----------------------------------------------------------
