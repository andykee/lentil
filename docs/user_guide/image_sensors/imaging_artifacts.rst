*****************
Imaging Artifacts
*****************

.. currentmodule:: lentil

.. note::

    To ensure accuracy and avoid introducing ailiasing artifacts, the input data should
    be at least Nyquist sampled.

Smear
=====
Smear is used to represent image motion with a relatively low temporal frequency relative
to integration time. The motion occurs in a slowly varying or fixed direction over one
integration time. Lentil's :func:`smear` method represents smear as a linear
directional blur over some distance (or number of pixels):

.. code:: pycon

    >>> import lentil
    >>> import matplotlib.pyplot as plt
    >>> psf = ...  # PSF calculation details omitted
    >>> psf_smear = lentil.convolvable.smear(psf, distance=5e-5,
    ...                                      pixelscale=5e-6,
    ...                                      oversample=3)
    >>> plt.subplot(121), plt.imshow(psf)
    >>> plt.subplot(122), plt.imshow(psf_smear)

.. image:: /_static/img/api/convolvable/smear_px.png
    :width: 500px

As an alternative to specifying physical distance and pixelscale, a number of pixels can also
be provided:

.. code:: pycon

    >>> import lentil
    >>> import matplotlib.pyplot as plt
    >>> psf = ...  # PSF calculation details omitted
    >>> psf_smear = lentil.convolvable.smear(psf, distance=10,
    ...                                      oversample=3)
    >>> plt.subplot(121), plt.imshow(psf)
    >>> plt.subplot(122), plt.imshow(psf_smear)

.. image:: /_static/img/api/convolvable/smear_px.png
    :width: 500px

The default behavior is to choose a new random smear direction each time :func:`smear`
is called, but a static direction can optionally be specified as needed:

.. code:: pycon

    >>> import lentil
    >>> import matplotlib.pyplot as plt
    >>> psf = ...  # PSF calculation details omitted
    >>> psf_smear = lentil.convolvable.smear(psf, distance=50,
    ...                                      angle=30)
    >>> plt.subplot(121), plt.imshow(psf)
    >>> plt.subplot(122), plt.imshow(psf_smear)

.. image:: /_static/img/api/convolvable/smear_m.png
    :width: 500px


Jitter
======
Jitter is used to represent image motion with a relatively high temporal frequency relative
to integration time. Lentil's :func:`jitter` method represents jitter with a
Gaussian blur operation. Note this approach is only valid if the motion occurs randomly in all
directions during one integration time.

.. code:: pycon

    >>> import lentil
    >>> import matplotlib.pyplot as plt
    >>> psf = ...  # PSF calculation details omitted
    >>> psf_jitter = lentil.convolvable.jitter(psf, scale=2)
    >>> plt.subplot(121), plt.imshow(psf)
    >>> plt.subplot(122), plt.imshow(psf_jitter)

.. image:: /_static/img/api/convolvable/jitter_px.png
    :width: 500px

If the jitter direction is not sufficiently random during a typical integration time, a timeseries
should be used instead. Note this will have a major impact on propagation performance but will
provide the most accurate results.

Pixel MTF
=========
A focal plane array samples a continuous light field to produce a digital image. Because
Lentil models diffraction numerically by propagating a finite set of points through an
optical system, the discretely sampled image plane intensity must be convolved with the
pixel's aperture function to accurately represent the intensity signal seen by each
pixel. Lentil's :func:`pixel` method implements this convolution.
After convolving the image plane intensity with the pixel MTF, the data should be
resampled to native detector sampling using :func:`rescale`. The
:func:`detector.pixelate` method combines the convolution and resampling operations
into a single method.
