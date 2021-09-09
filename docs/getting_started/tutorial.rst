.. _tutorial:

************************
Getting started tutorial
************************

In this short example, we'll walk through the steps to develop a very simple Lentil
model of an imaging system with a single pupil and propagate light to an image plane. 
The imaging system has a 1m diameter primary mirror, a secondary mirror obscuration 
of 0.33m centered over the primary, a focal length of 10m, and a focal plane with 
5um pixels.

.. image:: /_static/img/telescope.png
    :width: 650 px
    :align: center

First, we'll import Lentil and matplotlib:

.. code-block:: pycon

    >>> import matplotlib.pyplot as plt
    >>> import lentil

Now we can define the system amplitude and plot it:

.. code-block:: pycon

    >>> amplitude = lentil.circle(shape=(256, 256), radius=128) -
    ...             lentil.circle(shape=(256, 256), radius=128/3)
    >>> plt.imshow(amplitude)

.. image:: /_static/img/getting_started_amp.png
    :width: 350px

We'll create some Zernike-based wavefront error:

.. code-block:: pycon

    >>> coeffs = [0, 0, 0, 300e-9, 50e-9, -100e-9, 50e-9]
    >>> opd = lentil.zernike_compose(mask=amplitude, coeffs=coeffs)
    >>> plt.imshow(opd)

.. image:: /_static/img/getting_started_opd.png
    :width: 350px

Next we'll define the system's pupil plane. Note that the
:attr:`~lentil.Pupil.pixelscale` attribute represents the physical sampling of each
pixel in the pupil (in meters/pixel). Because our amplitude has a diameter of 256 pixels
and the system diameter was specified as 1m, the pixelscale is 1/256.

.. code-block:: pycon

    >>> pupil = lentil.Pupil(amplitude=amplitude, phase=opd, diameter=1, focal_length=10,
    ...                      pixelscale=1/256)

We will create a monochromatic :class:`~lentil.Wavefront` with wavelength of 650nm, 
propagate the wavefront through the pupil plane, and finally on to an image plane with
5um pixels. We'll also oversample the image plane by a factor of 10.

.. code-block:: pycon

    >>> w = lentil.Wavefront(wavelength=650e-9)
    >>> w = w * pupil
    >>> w = lentil.propagate_image(w, pixelscale=5e-6, npix=32, oversample=10)
    >>> plt.imshow(w.intensity)

.. image:: /_static/img/getting_started_psf_oversample.png
    :width: 350px

Finally, we will rescale the oversampled image to native sampling and include the
blurring effects of the pixel MTF:

.. code-block:: pycon

    >>> img = lentil.detector.pixellate(w.intensity, oversample=10)
    >>> plt.imshow(img)

.. image:: /_static/img/getting_started_psf_detector.png
    :width: 350px
