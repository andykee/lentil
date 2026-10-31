First Steps
===========

Lentil provides objects for defining optical planes and combining them to 
model diffraction and other imaging effects. It can be used to create models 
with varying levels of complexity - from simple prototypes to extremely 
high-fidelity models of ground or space telescopes.

Lentil uses one or more :ref:`Plane<user.planes>` objects to numericlly model 
an optical system. A :ref:`Wavefront<user.wavefront>` object representing a 
discretely sampled monochromatic plane wave is then propagated from plane 
to plane using one of the available :ref:`numerical diffraction propagation 
routines<user.diffraction>`. The resulting complex field can be visualized, 
analyzed, or further propagated through an optical model. Finally, Lentil 
includes a number of additional tools for :ref:`representing common imaging 
artifacts<user.artifacts>`, :ref:`modeling focal plane arrays
<user.image_sensors>` including many common noise sources, and 
:ref:`working with radiometric quantities<user.radiometry>`.

Below, we'll walk through the steps to develop a very simple Lentil model of 
an imaging system with a single pupil and image plane. The imaging system has 
a 1 m diameter primary mirror, a secondary mirror obscuration of 0.33 m 
centered over the primary, a focal length of 20 m, and a focal plane with 
5 um pixels.

Once installed, Lentil can be imported  with:

.. code-block:: python

    import lentil

This example also uses `matplotlib <https://matplotlib.org>`_ for plotting

.. code-block:: python

    import matplotlib.pyplot as plt

Now we can define the system amplitude and plot it:

.. plot::
    :context: reset
    :include-source:

    circ = lentil.circle(shape=(256, 256), radius=120)
    obsc = lentil.circle(shape=(256, 256), radius=120/3)
    amplitude = circ - obsc

.. plot::
    :context: close-figs
    :scale: 50

    fig, ax = plt.subplots()
    ax.imshow(amplitude, cmap='gray', origin='lower',
              extent=[-128, 127, -128, 127])

We'll create some wavefront error using a small number of Zernike modes:

.. plot::
    :context:
    :include-source:
    :nofigs:

    coeffs = [0, 0, 0, 300e-9, 50e-9, -100e-9, 50e-9]
    opd = lentil.zernike_compose(mask=amplitude, coeffs=coeffs)
    
.. plot::
    :context: close-figs
    :scale: 50

    fig, ax = plt.subplots(figsize=(3, 2.5))
    im = ax.imshow(opd, cmap=opd_cmap, origin='lower',
                   extent=[-128, 127, -128, 127])
    fig.colorbar(im, ax=ax, shrink=0.98)

Next we'll define the system's pupil plane. Note that the ``pixelscale`` 
attribute represents the physical sampling of each pixel in the pupil (in 
meters/pixel). Because our amplitude has a diameter of 256 pixels and the 
system diameter was specified as 1m, the pixelscale is 1/256.

.. plot::
    :context:
    :include-source:
    :nofigs:

    pupil = lentil.Pupil(amplitude=amplitude, opd=opd, diameter=1,
                         focal_length=20, pixelscale=1/256)

Here we create a ``Wavefront`` and propagate it through the pupil plane we
created above by multiplying the two objects together:

.. plot::
    :context:
    :include-source:
    :nofigs:

    wf = lentil.Wavefront(wavelength=650e-9)
    wf = wf * pupil

Now we can propagate the wavefront to the image plane using a discrete 
Fourier transform and plot the resulting image plane intensity:

.. plot::
    :context:
    :include-source:
    :nofigs:

    wf = lentil.propagate_dft(wf, shape=(64,64), pixelscale=5e-6, oversample=2)

.. plot::
    :context: close-figs
    :scale: 50

    fig, ax = plt.subplots()
    ax.imshow(wf.intensity, cmap=psf_cmap, origin='lower',
              extent=[-0.32, 0.315, -0.32, 0.315])
    ax.set_xlabel('mm')
    ax.set_ylabel('mm')

With a little extra code, we can perform a broadband propagation including
many wavelengths in the simulation:

.. plot::
    :context:
    :include-source:
    :nofigs:

    img = 0
    for wl in np.arange(450e-9, 650e-9, 5e-9):
        wf = lentil.Wavefront(wl)
        wf = wf * pupil
        wf = lentil.propagate_dft(wf, shape=(64,64), pixelscale=5e-6, oversample=2)
        img += wf.intensity

.. plot::
    :context: close-figs
    :scale: 50

    fig, ax = plt.subplots()
    ax.imshow(img, cmap=psf_cmap, origin='lower',
              extent=[-0.32, 0.315, -0.32, 0.315])
    ax.set_xlabel('mm')
    ax.set_ylabel('mm')

Finally, let's convolve the oversampled image with the MTF of a square pixel
and resample the image to native detector sampling:

.. plot::
    :context:
    :include-source:
    :nofigs:

    img_det = lentil.detector.pixelate(img, oversample=2)

.. plot::
    :context: close-figs
    :scale: 50

    fig, ax = plt.subplots()
    ax.imshow(img_det, cmap=psf_cmap, origin='lower',
              extent=[-0.32, 0.315, -0.32, 0.315])
    ax.set_xlabel('mm')
    ax.set_ylabel('mm')
