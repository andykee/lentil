.. _user_guide.diffraction:

.. |Wavefront| replace:: :class:`~lentil.Wavefront`
.. |Plane| replace:: :class:`~lentil.Plane`
.. |Pupil| replace:: :class:`~lentil.Pupil`
.. |Image| replace:: :class:`~lentil.Image`
.. |Detector| replace:: :class:`~lentil.Detector`

********************
Modeling diffraction
********************
Lentil uses Fourier transform-based algorithms to numerically model the propagation of an
electromagnetic field through an optical system. The electromagnetic field is represented
by a |Wavefront| object which stores the complex amplitude of the field at a discretely
sampled grid of points. The optical system is represented by a set of user-defined |Plane|
objects which interact with the wavefront as it propagates through the system.

.. note::

    This section of the User Guide assumes an undergraduate-level understanding of
    physical and Fourier optics. In-depth mathematical background and an extensive
    discussion of the validity of each diffraction approximation is available in [1]_.

How Wavefront works
===================


Numerical diffraction propagation
=================================
Lentil numerically models diffraction by propagating a |Wavefront| object through
any number of |Plane| objects representing an optical system. This propagation always
follows the same basic flow:

1. **Create a new wavefront** - A |Wavefront| represents a monochromatic, discretely
   sampled complex field that may propagate through space. By default, when a new
   |Wavefront| is constructed it represents an infinite plane wave (``1+0j``). Note
   that "infinite" in this case really just means that the wavefront field is
   `broadcastable <https://numpy.org/doc/stable/user/basics.broadcasting.html>`_ to
   any shape necessary.

    .. code-block:: pycon

        >>> w1 = lentil.Wavefront(wavelength=500e-9)
        >>> w1.field
        1+0j
        >>> w1.focal_length
        inf

2. **Propagate the wavefront through the first plane in the optical system** - the
   planes that describe an optical system typically change a propagating wavefront
   in some way. By multiplying a |Wavefront| and a |Plane| together, a new
   |Wavefront| is returned representing the state of the complex field after
   propagating through the plane:

    .. code-block:: pycon

        >>> pupil = lentil.Pupil(amplitude=lentil.circle((256, 256), 120),
        ...                      pixelscale=1/240, focal_length=10)
        >>> w2 = w1 * pupil

    Note the complex field of ``w2`` now clearly shows the effect of propagating through the
    circular aperture of ``pupil``:

    .. code-block:: pycon

        >>> plt.imshow(np.abs(w2.field))


    Additionally, because ``w2`` was propagated through a |Pupil| plane, it has inherited the
    pupil's focal length:

    .. code-block:: pycon

        >>> w2.focal_length
        10

    It is also possible to perform the multiplication in-place, reducing the memory footprint
    of the propagation:

    .. code-block:: pycon

        >>> w1 *= pupil

    .. note::

        Additional details on the plane-wavefront interaction are given in
        :ref:`user_guide.optical_systems.plane_wavefront`.

3. **Propagate the wavefront to the next plane in the optical system** - the |Wavefront|
   object provides a number of methods to propagate between planes. The appropriate method
   should be chosen based on the plane types the wavefront is propagating between.

   ======= ======= =========================================
   From    To      Method
   ======= ======= =========================================
   |Pupil| |Image| :func:`~lentil.Wavefront.propagate_image`
   |Image| |Pupil| :func:`~lentil.Wavefront.propagate_pupil`
   |Pupil| |Pupil| N/A
   |Image| |Image| N/A
   ======= ======= =========================================

   Propagations are defined by the following attributes:

   * :attr:`pixelscale` - the spatial sampling of the output plane
   * :attr:`npix` - the shape of the output plane
   * :attr:`npix_prop` - the shape of the propagation plane. See
     :ref:`user_guide.diffraction.npix` for additional details.
   * :attr:`oversample` - the number of times to oversample the output plane.
     See the section on :ref:`user_guide.diffraction.sampling` for more
     details.


   For example, to propagate a |Wavefront| from a |Pupil| to an |Image| plane:

   .. code-block:: pycon

        >>> w3 = w2.propagate_image(pixelscale=5e-6, npix=64, oversample=5)
        >>> plt.imshow(w3.intensity)

   .. note::

        When propagating between like planes (pupil to pupil or image to image),
        no additional propagation step must be taken.

4. **Repeat steps 2 and 3 until the propagation is complete** - if multiple planes
   are required to model the desired optical system, steps 2 and 3 should be
   repeated until the |Wavefront| has been propagated through all of the planes.

Polychromatic or broadband propagations
---------------------------------------
The steps outlined above propagate a single monochromatic |Wavefront| through an
optical system. The example below performs the same operation for a number of
different wavelengths and accumulates the resulting image plane intensity:

.. code-block:: python

    import numpy as np
    import lentil

    pupil = lentil.Pupil(amplitude=lentil.circle((256, 256), 120),
                         pixelscale=1/240, focal_length=10)

    wavelengths = [450e-9, 500e-9, 550e-9, 600e-9, 650e-9]
    img = np.zeros((320,320))

    for wl in wavelengths:
        w = lentil.Wavefront(wl)
        w *= pupil
        w.propagate_image(pixelscale=5e-6, npix=64, oversample=5)
        img += w.intensity

Note that the output ``img`` array must be sized to accommodate the oversampled
wavefront intensity given by ``npix`` * ``oversample``.

The ``propagate()`` method
--------------------------
In addition to the more direct/manual propagation approach described above, Lentil
provides a :func:`~lentil.propagate` method that will propagate a |Wavefront| (or
number of |Wavefront| objects in the case of a broadband propagation) through a
list of planes, selecting the appropriate propagator between planes based on
plane type.

.. image:: /_static/img/propagate.png
    :width: 800px
    :align: center

This method may struggle to do the right thing for more complicated
propagations, but it is very well tested and should be the default choice for
propagations including:

* Propagating from a |Pupil| to an |Image| or |Detector| plane
* Propagating from a |Pupil| to an |Image| or |Detector| plane including
  additional :class:`~lentil.Tilt` planes or other
  :ref:`user_guide.planes.transformations`
* Propagating from a |Pupil| to an |Image| or |Detector| plane including
  :ref:`user_guide.planes.dispersive`.

When using :func:`~lentil.propagate`, the propagation is defined by the following
arguments:

* :attr:`planes` - a list of |Plane| objects representing an unfolded optical system
* :attr:`wave` - the propagation wavelength. A single value results in a monochromatic
  propagation. If a list of wavelengths are provided, multiple monochromatic propagations
  will be run - once for each wavelength in ``wave``. used for the propagation.
* :attr:`weight` - the weight associated with each wavelength in :attr:`wave`. Note that
  weights can be either relative or absolute depending on the use case.

Additional arguments provide further customization of the propagation behavior:

* :attr:`npix` - the shape of the result.
* :attr:`npix_prop` - the shape of the propagation result. Note that if
  ``npix_prop`` is not specified, Lentil uses the value provided for ``npix``.
* :attr:`oversample` - the number of times to oversample the output plane. See the
  sectionon :ref:`user_guide.diffraction.sampling` for more details.
* :attr:`rebin` - specifies whether to return the output plane to its native sampling or
  return the oversampled result.
* :attr:`fit_tilt` - specifies the tilt handling strategy. See
  :ref:`user_guide.diffraction.tilt` for more details.
* :attr:`min_q` - specifies the minimum allowable propagation sampling term :math:`q`.
  Planes may be resampled as needed to satisfy ``min_q``. If ``min_q = 0`` (default),
  no resampling is done. See :ref:`user_guide.diffraction.sampling` for more details.
* :attr:`flatten` - specifies whether to collapse wavelength-specific output planes to a
  single array or return a 3D cube of results

.. _user_guide.diffraction.npix:

``npix`` vs ``npix_prop``
-------------------------
Lentil's propagation methods have two arguments for controlling the shape of
the propagation output: ``npix`` and ``npix_prop``.

``npix`` specifies the shape of the entire output plane while ``npix_prop``
specifies the shape of the propagation result. If ``npix_prop`` is not
specified, it defaults to ``npix``. The propagation result is placed in the
appropriate location in the (potentially larger) output plane when a |Wavefront|
:attr:`~lentil.Wavefront.field` or :attr:`~lentil.Wavefront.intensity`
attribute is accessed.

.. image:: /_static/img/propagate_npix_chip.png
    :width: 450px
    :align: center

It can be advantageous to specify ``npix_prop`` < ``npix`` for performance
reasons, although care must be taken to ensure needed data is not accidentally
left out:

.. image:: /_static/img/npix_prop.png
    :width: 500px
    :align: center

For most pupil to image plane propagations, setting ``npix_prop`` to 128 or 256
pixels provides an appropriate balance of performance and propagation plane size.

For image to pupil plane propagations, ``npix_prop`` must be sized to ensure
the pupil extent is adequately captured. Because the sampling constraints on
image to pupil plane propagations are typically looser, it is safest to let
``npix_prop`` default to the same value as ``npix``.

Discrete Fourier transform algorithms
-------------------------------------
Most diffraction modeling tools use the Fast Fourier Transform (FFT) to evaluate the
discrete Fourier transform (DFT) when propagating between planes. While the FFT provides
great computational and memory efficiency, high-fidelity optical simulations may require
working with exceptionally large zero-padded arrays to satisfy the sampling requirements
imposed by the FFT.

Lentil implements a more general form of the DFT sometimes called the matrix triple
product (MTP DFT) to perform the Fourier transform to propagate between planes. While the
MTP DFT is slower than the FFT for same sized arrays, the MTP DFT provides independent
control over the input and output plane sizing and sampling. This flexibility makes the
MTP DFT ideally suited for performing propagations to discretely sampled image planes
where it is often necessary to compute a finely sampled output over a relatively small
number of pixels.

The chirp Z-transform provides additional efficiency when transforming large arrays.
Lentil selects the most appropriate DFT method automatically based on the plane size and
sampling requirements.

.. _user_guide.diffraction.sampling:

Sampling considerations
=======================

.. image:: /_static/img/discrete_Q_sweep.png
    :width: 800px
    :align: center

.. image:: /_static/img/q_sweep.png
    :width: 800px
    :align: center


.. image:: /_static/img/propagate_fourier_period.png
    :width: 550px
    :align: center

.. _user_guide.diffraction.tilt:

Working with large tilts
========================
.. image:: /_static/img/propagate_tilt_phase.png
    :width: 450px
    :align: center


.. image:: /_static/img/propagate_tilt_phase_wrap.png
    :width: 650px
    :align: center


.. image:: /_static/img/propagate_tilt_angle.png
    :width: 600px
    :align: center


.. image:: /_static/img/propagate_tilt_angle_steps.png
    :width: 600px
    :align: center

Differences for segmented apertures
===================================




.. [1] Goodman, *Introduction to Fourier Optics*.:w
