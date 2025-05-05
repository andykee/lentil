Release notes
=============

v0.8.5
------
Released May 5, 2025

* New ``Plane.freeze()`` method
* Documentation updates
* Updates to GitHub Actions
* Add support for NumPy versions >= 2.0

v0.8.4
------
Released December 2, 2024

* New subarray() function
* Minor bugfixes and documentation updates

v0.8.3
------
Released June 15, 2024

* Minor bugfix and documentation updates

v0.8.2
------
Released April 20, 2024

* New propagate_dft() functionality for performing propagations
  in areas limited by a mask
* Plane directly sets its own internal attribues `#53`_
* Migrate to numpy random.default_rng() for random numbers
* Backend enhancements and bugfixes

.. _#53: https://github.com/andykee/lentil/issues/53

v0.8.1
------
Released February 27, 2024

* Deprecate :class:`~lentil.Image` ``shape`` parameter `#50`_
* :class:`~Plane` ``mask``, ``pixelscale``, ``diameter``, and ``ptype`` are 
  now immutable
* Fix bugs in creation and handling of :class:`~lentil.Plane` ``mask`` and 
  ``slice``
* Fix scipy deprecation warnings
* Update GitHub actions `#51`_ `#52`_
* Update docs `#49`_

.. _#49: https://github.com/andykee/lentil/issues/49
.. _#50: https://github.com/andykee/lentil/issues/50
.. _#51: https://github.com/andykee/lentil/issues/51
.. _#52: https://github.com/andykee/lentil/issues/52

v0.8.0
------
Released February 2, 2024


* New :func:`~lentil.propagate_dft` method for performing far-field
  propagations using the matrix triple product DFT
* New :func:`~lentil.propagate_fft` method for performing far-field 
  propagations using the FFT `#46`_
* Deprecate ``Wavefront.propagate_image()``
* Propagation backend updates and enhancements
* DFT implementation improvements
* Plane ``phase`` attribute renamed to ``opd``
* Allow ``amp`` as an alias for ``amplitude`` in Plane constructor
* New ``ptype`` object for managing plane and wavefront types
* Enforce ``ptype``-based rules in :func:`Plane.multiply`
* Streamline :class:`~lentil.Wavefront` constructor
* New :class:`DispersiveTilt` object
* Deprecate :class:`Grism` -- :class:`DispersiveTilt` should be used instead
* Utility shapes now accept an ``antialias`` argument
* New :func:`~lentil.spider` method for drawing spiders
* New :func:`~lentil.hex_segments` method for drawing segmented apertures made
  of rings of hexagonal segments
* Standardize shape naming
* Deprecate in-place operations `#43`_
* Many documentation updates and improvements

.. _#43: https://github.com/andykee/lentil/issues/43
.. _#46: https://github.com/andykee/lentil/issues/46

v0.7.0
------
Released March 7, 2022

* Fix complex amplitude sign flip introduced in `v0.6.0`_
* Remove unused parameter from ``Wavefront.insert()`` function
  signature `#42`_
* Scipy compatibility - Fix Scipy map_coordinates import `#40`_
* Python 3.9 compatibility - Ensure ``math.factorial()`` always
  receives an int

.. _#40: https://github.com/andykee/lentil/issues/40
.. _#42: https://github.com/andykee/lentil/issues/42

v0.6.0
------
Released January 21, 2022

* Entirely new approach to how diffraction propagations are performed:

  * New ``propagate_image()`` method for propagating between Pupil and
    Image planes

  * Deprecate core ``propagate()`` method

  * Include negative sign in complex phasor complex exponential

* Wavefront complex field data is now managed using a new internal Field
  class
* Standardize around (row, col) aka. ij indexing
* New methods for Plane resampling (``Plane.resample()``) and rescaling
  (``Plane.rescale()``)
* Collapse Plane segmask and mask functionality `#24`_
* Allow in-place operations on Wavefront `#38`_
* Relocate contents of ``zerenike``, ``wfe``, ``convolvable``, and ``util``
  modules to the core ``lentil`` namespace
* Allow floating point plane masks, which are automatically cast to bool
* Documentation updates
* Extend unit test coverage slightly
* Switch to GitHub Actions for unit testing and code coverage

.. _#24: https://github.com/andykee/lentil/issues/24
.. _#38: https://github.com/andykee/lentil/issues/38

v0.5.0
------
Released August 13, 2021

* Propagations with ``tilt='angle'`` have tilt projected out of each
  plane once before the entire propagation rather than at each monochromatic
  propagation
* Rework ``Plane.pixelscale`` to always store (r,c) pixelscale
* Fix bug in ``Plane.mask`` on the fly calculation that was overwriting
  ``Plane.amplitude`` with a binary mask
* No longer cache ``Plane.ptt_vector``
* Deprecate ``Plane.cache_propagate()`` and ``Plane.clear_cache_propagate()``.
  This functionality has been migrated to ``propagate._prepare_planes()``
  and ``propagate._cleanup_planes()``
* New ``Plane.rescale()`` method to rescale Plane pixelscale
* Update ``util.rescale()`` to choose a more conservative (better sampled)
  result when having to choose an integer output shape
* Define ``Wavefront.__slots__`` to increase attribute access speed and reduce
  memory footprint
* ``util.circle()`` `center` parameter is now called `shift`
* Deprecate ``cache.Cache`` in favor of a simple dictionary
* New function ``fourier.expc()`` to more quickly compute a complex exponential
* ``fourier.dft2()`` now accepts an offset parameter
* New function ``Plane.fit_tilt()`` to handle tilt fitting and removal of in the
  Plane's ``phase`` attribute. This is now called once
* New function ``Plane.slice()`` for computing avaliable slices from the plane
  attributes to speed up propagation performance
* New ``Detector()`` plane that returns intensity
* Update ``zernike.zernike_coordinates()`` to automatically compute shift that
  locates the origin at the mask centroid if no shift is provided.

v0.4.1
------
Released October 7, 2020

* Fix implementation error in Grism model dispersion calculations

v0.4.0
------
Released October 6, 2020

.. note::

  The Grism model updates are broken in this release. It has been yanked from
  PyPi. The issue is fixed in v0.4.1.

* Update Grism model to use correct definition of dispersion, accomodate
  trace and dispersion models with polynomial order > 1
* Establish coordinate system `#12`_
* Fix direction and orientation of Tilt `#12`_
* Allow spectral inputs to radiometry.path_emission

.. _#12: https://github.com/andykee/lentil/issues/12


v0.3.4
------
Released September 8, 2020

* Fix implementation error in Gaussian detector.shot_noise
* Add better exception handling for detector.shot_noise `#10`_
* No longer check Python version on import
* Update ``np.ediff1d`` usage to be compatible with Numpy 1.19

.. _#10: https://github.com/andykee/lentil/issues/10

v0.3.3
------
Released August 17, 2020

* Make FPN seed optional in ``detector.dark_current``

v0.3.2
------
Released July 20, 2020

* Update ``detector.adc`` to prevent negative values from being returned.

v0.3.1
------
Released July 16, 2020

* Imaging artifact classes have been rearchitected as functions.
* Legacy functionality from the ``detector.Windowable`` class has been resurrected into
  :func:`lentil.util.window`
* Deprecate ``util.col_major_to_util_major()``
* Lentil is now compatible with Python 3.6 and newer. `#9`_

.. _#9: https://github.com/andykee/lentil/issues/9

v0.3.0
------
Released July 8, 2020

* The Plane attribute caching approach has been entirely reworked, eliminating the need
  for end-users to explicitly decorate attributes defined in subclasses:

  * Users are now able to explicitly choose which attributes are cached when
    ``cache_propagate()`` is called by specifying them in ``Plane.cache_attrs``. The
    only accepted values right now are ``amplitude`` and ``phase``. Note that
    ``ptt_vector`` is always cached and is not allowed to be specified in
    ``cache_attrs``.

  * Rather than checking for and returning cached values at the attribute getter level,
    it is now done inside ``Plane.multiply()``. This change streamlines both the plane
    attribute getter code and the creation of planes with phase attributes that should
    be random with each access.

  * The ``cache_propagate`` decorator has been deprecated, and the documentation and
    tests have been updated to reflect the changes in functionality. `#7`_

* Fix bug in ``zernike_coordinates`` that was causing modes over off-centered masks to
  be incorrectly computed. `#8`_
* Change default behavior of ``zernike_basis`` to return a stack of matrices rather than
  a single vectorized matrix.

.. _#7: https://github.com/andykee/lentil/issues/7
.. _#8: https://github.com/andykee/lentil/issues/8

v0.2.0
------
Released June 29, 2020

* Collapse ``Detector`` and ``Image`` planes into single ``Image`` plane. The pupil to
  image plane propagation method is now chosen based on whether the ``Image`` plane has
  a defined ``pixelscale`` (propagate via matrix triple product DFT) or if
  ``pixelscale`` is None (propagate via FFT - eventually). ``Detector`` class has been
  deprecated. `#5`_
* Deprecate ``FPA`` and ``BayerFPA``. Some functionality has been retained but converted
  to atomic functions in the :ref:`detector<api.detector>` module. `#6`_
* Completely rework the contents of the :ref:`detector<api.detector>` module. All
  objects have been deprecated. Some functionality has been retained but converted to
  atomic functions instead. `#6`_
* Deprecate ``util.coordinates``
* Change the way ``Rotate`` angle is interpreted to behave more intuitively
* A number of small bugfixes and enhancements
* Updated documentation
* More unit tests

.. _#5: https://github.com/andykee/lentil/issues/5
.. _#6: https://github.com/andykee/lentil/issues/6

v0.1.1
------
Released June 21, 2020

* Update ``propagate`` to support :class:`~lentil.Tilt` planes `#1`_
* Streamline the innards of :func:`~lentil.propagate`
* Update :func:`lentil.wfe.power_spectrum` to return phases with a slightly more correct
  RMS
* Remove unused code
* Increase unit testing coverage
* Set up Travis CI, Coveralls

.. _#1: https://github.com/andykee/lentil/issues/1

v0.1.0
------
Released June 12, 2020

* Initial public release
