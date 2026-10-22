Imaging Artifacts
=================

.. note::

    To ensure accuracy and avoid introducing ailiasing artifacts, the input
    data should be at least Nyquist sampled when applying imaging artifacts
    via convolution.

.. plot::
    :context: reset

    mask = lentil.circle((256, 256), 120)
    opd = lentil.zernike_compose(mask, [0, 0, 0, -300e-9, 50e-9, -100e-9, 50e-9])
    pupil = lentil.Pupil(amplitude=mask, opd=opd, focal_length=10, pixelscale=1 / 240)
    w = lentil.Wavefront(650e-9)
    w *= pupil
    w = lentil.propagate_dft(w, pixelscale=5e-6, shape=64, oversample=5)
    psf = w.intensity

Jitter
------
Jitter is used to represent image motion with a high temporal frequency
relative to integration time (i.e. motion occurring within an integration 
time). Lentil's :func:`jitter` method represents jitter with a Gaussian blur
operation. Note this approach is only valid if the motion occurs randomly in
all directions during one integration time.

.. plot::
    :context: close-figs
    :include-source:
    :scale: 50

    psf_jitter = lentil.jitter(psf, scale=2, oversample=5)

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(4,2))
    ax[0].imshow(psf, cmap='turbo')
    ax[0].set_title('Input image'); ax[0].axis('off')
    ax[1].imshow(psf_jitter, cmap='turbo')
    ax[1].set_title('Jitter image'); ax[1].axis('off')

If the jitter being modeled is not sufficiently random during a typical integration time, a
timeseries should be used instead. This can have a major impact on propagation performance
but will provide the most accurate results.

Smear
-----
Smear is used to represent image motion with a relatively low temporal frequency relative
to integration time. The motion occurs in a slowly varying or fixed direction over one
integration time. Lentil's :func:`smear` method represents smear as a linear
directional blur over some physical distance or number of pixels:

.. plot::
    :context: close-figs
    :include-source:
    :scale: 50

    psf_smear = lentil.smear(psf, distance=5e-5, pixelscale=5e-6,
                             oversample=5)

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(4,2))
    ax[0].imshow(psf, cmap='turbo')
    ax[0].set_title('Input image'); ax[0].axis('off')
    ax[1].imshow(psf_smear, cmap='turbo')
    ax[1].set_title('Smeared image'); ax[1].axis('off')


As an alternative to specifying physical distance and pixelscale, a number of pixels can also
be provided:

.. plot::
    :context: close-figs
    :include-source:
    :scale: 50

    psf_smear = lentil.smear(psf, distance=10.2, oversample=5)

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(4,2))
    ax[0].imshow(psf, cmap='turbo')
    ax[0].set_title('Input image'); ax[0].axis('off')
    ax[1].imshow(psf_smear, cmap='turbo')
    ax[1].set_title('Smeared image'); ax[1].axis('off')

The default behavior is to choose a new random smear direction each time :func:`smear`
is called, but a static direction can optionally be specified as needed:

.. plot::
    :context: close-figs
    :include-source:
    :scale: 50

    psf_smear = lentil.smear(psf, distance=25, oversample=5,
                             angle=30)

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(4,2))
    ax[0].imshow(psf, cmap='turbo')
    ax[0].set_title('Input image'); ax[0].axis('off')
    ax[1].imshow(psf_smear, cmap='turbo')
    ax[1].set_title('Smeared image'); ax[1].axis('off')
