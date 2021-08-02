.. _diffraction.overview:

.. |Wavefront| replace:: :class:`~lentil.wavefront.Wavefront`
.. |Plane| replace:: :class:`~lentil.Plane`

********
Overview
********

This section of the User Guide assumes an undergraduate-level understanding of
physical and Fourier optics. In-depth mathematical background and an extensive
discussion of the validity of each diffraction approximation is available in [1]_.

.. Theoretical background
.. ======================
.. 
.. Maxwell's equations describe how an electromagnetic field propagates through free space
.. and behaves when encountering an obstruction or aperture (diffraction). Practical
.. solutions to Maxwell's equations for common optics problems are made possible by making a
.. few key assumptions:
..
.. * `Scalar diffraction theory` assumes that the propagating electromagnetic field can be
..   treated treated as a scalar field (thus ignoring the vector nature of the field and any
..   polarization effects that may be present).
.. * The `paraxial approximation` assumes the electromagnetic field propagates through an
..   optical system in a direction geneally aligned with the optical axis.
..
.. These two assumptions form the basis for the `Fresnel` or `near-field` approximation for
.. modeling diffraciton. By assuming the electromagnetic field is observed at a sufficiently
.. large distance beyond the diffracting obscuration or if the optical system imparts a
.. quadratic phase term (by a focusing lens, for example), we are able to use the simpler
.. `Fraunhofer` or `far-field` approximation for modeling diffraction.

.. note::

    Lentil attempts to provide a consistent interface to both the Fraunhofer and Fresnel
    algorithms. The remainder of this section of the user guide describes the Fraunhofer
    approximation but identifies where any differences may exist. In-depth documentation
    specific to the Fresnel approximation is available in the :ref:`diffraction.fresnel`
    section of the user guide.

Propagation interface
=====================

Lentil uses Fourier transform-based algorithms to numerically model the propagation of an
electromagnetic field through an optical system. The electromagnetic field is represented
internally by a |Wavefront| object which stores the complex amplitude of the field at a
discretely sampled grid of points. The optical system is represented by a set of
user-defined |Plane| objects which interact with the wavefront as it propagates through
the system.

Lentil's :func:`~lentil.propagate` method implements a number of algorithms to
numerically model the propagation of a wavefront through an optical system using the
Fraunhofer diffraction approximation. The propagation is defined by the following
arguments:

* :attr:`planes` - a list of |Plane| objects representing an unfolded optical system
* :attr:`wave` - the propagation wavelength. A single value results in a monochromatic
  propagation. If a list of wavelengths are provided, multiple monochromatic propagations
  will be run - once for each wavelength in ``wave``. used for the propagation.
* :attr:`weight` - the weight associated with each wavelength in :attr:`wave`. Note that
  weights can be either relative or absolute depending on the use case.

Additional arguments provide further customization of the propagation behavior:

* :attr:`npix` - the shape of the result.
* :attr:`npix_chip` - the shape of the propagation result. Note that if
  ``npix_chip`` is not specified, Lentil uses the value provided for ``npix``.
* :attr:`oversample` - the number of times to oversample the output plane. See the
  sectionon :ref:`diffraction.sampling` for more details.
* :attr:`rebin` - specifies whether to return the output plane to its native sampling or
  return the oversampled result.
* :attr:`tilt` - specifies the tilt handling strategy. See :ref:`diffraction.tilt` for
  more details.
* :attr:`interp_phasor` - specifies whether |Plane| objects should be resampled to avoid
  Fourier domain wraparound contamination. See :ref:`diffraction.sampling` for more
  details.
* :attr:`flatten` - specifies whether to collapse wavelength-specific output planes to a
  single array or return a 3D cube of results

Available planes
----------------
pupil
image

subclasses

described in planes section



.. _diffraction.npix_vs_npix_chip:

``npix`` vs. ``npix_chip``
--------------------------
:func:`~lentil.propagate` has two arguments for defining the shape of the propagation
output: ``npix`` and ``npix_chip``.

Note that ``npix`` specifies the entire
  output shape while ``npix_chip`` specifies the shape of the propagation result. This
  distinction exists because of how Lentil handles tilt.


.. image:: /_static/img/propagate_npix.png
    :width: 450px
    :align: center

.. image:: /_static/img/propagate_npix_chip.png
    :width: 450px
    :align: center

Propagation algorithm
=====================
The general propagation algorithm propagates a monochromatic wavefront through a set of
optical planes in the following way:

1. Create a new wavefront with the desired wavelength attribute. At this point, the
   wavefront represents a plane wave.

2. Iterate through the supplied list of planes. For each plane:

   a. Multiply the wavefront's complex amplitude by the plane's complex amplitude
   b. Propagate the wavefront to the next plane according to the following rules:

      * Pupil to image: apply a Fourier transform to the wavefront's complex amplitude
      * Image to pupil: apply an inverse Fourier transform to the wavefront's complex
        amplitude
      * Between two planes of the same type or if the current plane is a transformation:
        do nothing

    c. If the final plane is an image plane, compute the intensity of the resulting
       wavefront's complex amplitude and apply any specified weight

3. Accumulate the results and perform and requesdted post-processing

Graphically, this looks like

.. image:: /_static/img/propagate.png
    :width: 800px
    :align: center

Discrete Fourier transform algorithms
-------------------------------------
Most diffraction modeling tools use the Fast Fourier Transform (FFT) to evaluate the
discrete Fourier transform (DFT) when propagating between planes. While the FFT provides
great computational and memory efficiency, high-fidelity optical simulations may require
working with exceptionally large zero-padded arrays to satisfy the sampling requirements
imposed by the FFT.

Lentil implements a more general form of the DFT sometimes called the matrix triple
product (MTP DFT) to perform the Fourier transform to propagate between planes. While the
MTP DFT is slower than an FFT of the same sized array, the MTP DFT provides independent
control over the input and output plane sizing and sampling. This flexibility makes the
MTP DFT ideally suited for performing propagations to discretely sampled image planes
where it is often necessary to compute a finely sampled output over a relatively small
number of pixels.

The chirp Z-transform provides additional efficiency when transforming large arrays.
Lentil selects the most appropriate DFT method automatically based on the plane size and
sampling requirements.


References
==========

.. [1] Goodman, *Introduction to Fourier Optics*.
