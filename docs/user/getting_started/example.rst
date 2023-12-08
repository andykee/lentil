.. _user.getting_started.example:

****************
A simple example
****************
This is an example that demonstrates the core functionality of Lentil. It
is mainly written for new users. More detailed examples and specific design 
patterns are available :ref:`here<examples>`.

First, we'll import Lentil:

.. code-block:: pycon

    >>> import lentil

We'll also import `Matplotlib <https://matplotlib.org>`_ to visualize results:

.. code-block:: pycon

    >>> import matplotlib.pyplot as plt

Creating planes
===============
Most simple diffraction simulations can be represented by a single far-field 
propagation from a pupil plane to an image plane. First, we'll create a 
pupil amplitude map and corresponding optical path difference (OPD) map which
represents the wavefront error of the system. We'll construct the OPD map from 
a combination of Zernike modes:

.. plot::
    :context: reset
    :include-source:
    :scale: 50

    >>> amp = lentil.circle(shape=(256,256), radius=120)
    >>> coef = [0, 0, 0, 300e-9, 50e-9, -100e-9, 50e-9]
    >>> opd = lentil.zernike_compose(mask=amp, coeffs=coef)
    >>> fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(5, 4))
    >>> ax1.imshow(amp)
    >>> ax1.set_title('Amplitude')
    >>> ax2.imshow(opd)
    >>> ax2.set_title('OPD')

Now we can use the amplitude and OPD maps to construct a |Pupil| plane with
a focal length of 20 meters and a diameter of 1 meter:

.. plot::
    :context: close-figs
    :include-source:

    >>> pupil = lentil.Pupil(amplitude=amp, opd=opd, pixelscale=1/240, 
    ...                      focal_length=20)

Note the diameter is implicitly defined via the 
:attr:`~lentil.Pupil.pixelscale` attribute:

.. image:: /_static/img/pixelscale.png
    :width: 500px
    :align: center

.. note::

    Lentil is "unitless" in the sense that it doesn't enforce a specific base
    unit. All calculations are well behaved for both metric and imperial units.
    It is important that units are consistent however, and this task is left to
    the user.

    That being said, it is recommended that all calculations be performed in
    terms of either meters, millimeters, or microns.

Diffraction
===========

Pupil to image plane propagation
--------------------------------
The simplest diffraction propagation is from a pupil to image plane. Here, we
construct a |Wavefront| with wavelength of 500 nm, again represented
in meters:

.. plot::
    :context:
    :include-source:

    >>> w0 = lentil.Wavefront(wavelength=500e-9)

Next, we'll propagate the wavefront through the pupil plane we defined above.
Lentil uses multiplication represent the interaction between a |Plane| and
|Wavefront|:

.. plot::
    :context:
    :include-source:

    >>> w1 = w0 * pupil

Finally, we'll propagate the wavefront to a discreetly sampled image plane
using :func:`~lentil.propagate_dft`. In this case, we'll sample
the result on a grid with spacing of 5e-6 meters and perform the propagation 
2 times oversampled:

.. plot::
    :context:
    :include-source:

    >>> w2 = lentil.propagate_dft(w1, shape=(64,64), pixelscale=5e-6, oversample=2)

The resulting intensity (point spread function) can now be observed:

.. plot::
    :context:
    :include-source:
    :scale: 50

    >>> plt.imshow(w2.intensity)

Finally, we will rescale the oversampled image to native sampling and include the
blurring effects due to the discrete pixel sampling of the image plane:

.. plot::
    :context: close-figs
    :include-source:
    :scale: 50

    >>> img = lentil.detector.pixelate(w2.intensity, oversample=2)
    >>> plt.imshow(img)

.. Focal planes
.. ============


.. Radiometry
.. ==========


