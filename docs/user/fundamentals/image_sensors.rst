.. _user.fundamentals.image_sensors:

*************
Image sensors
*************
Lentil does not provide an integrated detector modeling capability, but instead provides
a collection of functions in the :ref:`detector module<api.detector>` to help model
image sensors and represent some of the more commonly encountered noise sources.

Signal flow
===========
A large body of technical literature exists describing both how modern image sensors
work as well as the various sources of noise that impact their overall performance.
The details of the signal flow through an image sensor and decisions about which specific
noise sources should be modeled depend on the application, but either of the following
two references provide a good starting point:

* `EVMA Standard 1288 <https://www.emva.org/standards-technology/emva-1288/>`_
* Janesick, *Photon Transfer*

.. currentmodule:: lentil

Charge collection
=================
Imaging sensors convert light intensity to an electrical signal. The ratio of electrons
collected to the number of photons striking a sensor's photoreactive surface is called
the quantum efficiency (commonly represented by :math:`\eta` or abbreviated QE). Because
both light intensity and QE are spectrally-dependent but electron count is not, for
performance reasons it is advantageous to begin working in units of electrons as early 
as possible. 


This process is slightly more complicated for a sensor with a 
`Bayer filter <https://en.wikipedia.org/wiki/Bayer_filter>`_.

Lentil's
:func:`detector.collect_charge` function models this charge collection process for a
monochromatic sensor and the :func:`detector.collect_charge_bayer` models the same process
for a sensor with a `Bayer filter <https://en.wikipedia.org/wiki/Bayer_filter>`_.

.. image:: /_static/img/bayer.png
    :width: 550px
    :align: center

There are two different approaches to performing charge collection:

1. Provide a datacube of photon counts or fluxes separated spectrally (i.e. a
   nwave x nrows x ncols cube) and a 

Performance considerations
--------------------------
blah

Pixel effects
=============
An image sensor is constructed from an array of pixels that sample a continuous light
field to produce a digital image. Because Lentil models diffraction numerically by
propagating a finite set of points through an optical system, the discretely sampled
image plane intensity must be convolved with the pixel's aperture function to accurately
represent the intensity signal sensed by each pixel. Lentil's :func:`detector.pixel`
function implements this convolution. After convolving the image plane intensity with
the pixel MTF (the sinc function), the data should be resampled to native detector 
sampling using :func:`rescale`. The :func:`detector.pixelate` function combines the 
convolution and resampling operations into a single method.

.. code-block:: pycon

    >>> import matplotlib.pyplot as plt
    >>> import lentil
    >>> psf = ...  # PSF calculation details omitted
    >>> psf_pixelate = lentil.detector.pixelate(psf, oversample=2)
    >>> plt.subplot(121), plt.imshow(psf, origin='lower')
    >>> plt.subplot(122), plt.imshow(psf_pixelate, origin='lower')

.. plot:: _img/python/pixelate.py
    :scale: 50

Noise sources
=============


Analog to digital conversion
============================


Cosmic rays
===========
Cosmic rays are high energy particles that travel through space. When they strike
an image sensor they create a signal that appears as a saturated pixel or streak of
saturated pixels. Lentil's :func:`detector.cosmic_rays` function simulates random
cosmic ray strikes during an integration with hit rates and fluxes taken from [1]_.



Imaging artifacts
=================

.. note::

    To ensure accuracy and avoid introducing ailiasing artifacts, the input data should
    be at least Nyquist sampled when applying imaging artifacts via convolution.

Smear
-----
Smear is used to represent image motion with a relatively low temporal frequency relative
to integration time. The motion occurs in a slowly varying or fixed direction over one
integration time. Lentil's :func:`smear` method represents smear as a linear
directional blur over some distance (or number of pixels):

.. code:: pycon

    >>> import lentil
    >>> import matplotlib.pyplot as plt
    >>> psf = ...  # PSF calculation details omitted
    >>> psf_smear = lentil.smear(psf, distance=5e-5,
    ...                          pixelscale=5e-6,
    ...                          oversample=5)
    >>> plt.subplot(121), plt.imshow(psf, origin='lower')
    >>> plt.subplot(122), plt.imshow(psf_smear, origin='lower')

.. plot:: _img/python/smear.py
    :scale: 50

As an alternative to specifying physical distance and pixelscale, a number of pixels can also
be provided:

.. code:: pycon

    >>> import lentil
    >>> import matplotlib.pyplot as plt
    >>> psf = ...  # PSF calculation details omitted
    >>> psf_smear = lentil.smear(psf, distance=10,
    ...                          oversample=5)
    >>> plt.subplot(121), plt.imshow(psf, origin='lower')
    >>> plt.subplot(122), plt.imshow(psf_smear, origin='lower')

.. plot:: _img/python/smear.py
    :scale: 50

The default behavior is to choose a new random smear direction each time :func:`smear`
is called, but a static direction can optionally be specified as needed:

.. code:: pycon

    >>> import lentil
    >>> import matplotlib.pyplot as plt
    >>> psf = ...  # PSF calculation details omitted
    >>> psf_smear = lentil.smear(psf, distance=25,
    ...                          angle=30)
    >>> plt.subplot(121), plt.imshow(psf, origin='lower')
    >>> plt.subplot(122), plt.imshow(psf_smear, origin='lower')

.. plot:: _img/python/smear_directional.py
    :scale: 50

Jitter
------
Jitter is used to represent image motion with a relatively high temporal frequency relative
to integration time. Lentil's :func:`jitter` method represents jitter with a
Gaussian blur operation. Note this approach is only valid if the motion occurs randomly in all
directions during one integration time.

.. code:: pycon

    >>> import lentil
    >>> import matplotlib.pyplot as plt
    >>> psf = ...  # PSF calculation details omitted
    >>> psf_jitter = lentil.jitter(psf, scale=2, oversample=5)
    >>> plt.subplot(121), plt.imshow(psf)
    >>> plt.subplot(122), plt.imshow(psf_jitter)

.. plot:: _img/python/jitter.py
    :scale: 50

If the jitter being modeled is not sufficiently random during a typical integration time, a
timeseries should be used instead. This can have a major impact on propagation performance
but will provide the most accurate results.


.. [1] Offenberg, J.D. et. al. Multi-Read Out Data Simulator. (2000).
