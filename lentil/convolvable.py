import numpy as np


def jitter(img, scale, pixelscale=1, oversample=1):
    """
    Apply image jitter via convolution.

    Parameters
    ----------
    img : array_like
        Input image
    scale : float
        1-sigma jitter motion

        If ``pixelscale = 1``, ``scale`` should give jitter in terms of
        fractional pixels. If ``pixelscale`` is the physical dimension of
        one pixel, ``scale`` should give jitter in terms of the same units.
    pixelscale : float, optional
        Pixel size. If ``pixelscale = 1`` (default), the ``scale`` parameter
        provided when Jitter is called should give jitter in terms of fractional
        pixels. If ``pixelscale`` is the physical dimension of one pixel, ``the
        ``scale`` parameter provided when Jitter is called should give jitter in
        terms of the same units.
    oversample : int, optional
        Oversampling factor of img. Default is 1.

    Returns
    -------
    out : ndarray
        Image with jitter applied.

    Examples
    --------
    Apply 2 pixels jitter to a 5x oversmapled PSF:

    .. plot::
        :context: reset

        mask = lentil.circle((256, 256), 120)
        opd = lentil.zernike_compose(mask, [0, 0, 0, -300e-9, 50e-9, -100e-9, 50e-9])
        pupil = lentil.Pupil(amplitude=mask, opd=opd, focal_length=10, pixelscale=1 / 240)
        w = lentil.Wavefront(650e-9)
        w *= pupil
        w = lentil.propagate_dft(w, pixelscale=5e-6, shape=64, oversample=5)
        psf = w.intensity

    .. plot::
        :context:
        :include-source:
        :nofigs:

        psf_jitter = lentil.jitter(psf, scale=2, oversample=5)

    .. plot::
        :context:
        :scale: 50

        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(4,2))
        ax[0].imshow(psf, cmap=psf_cmap)
        ax[0].set_title('Input image'); ax[0].axis('off')
        ax[1].imshow(psf_jitter, cmap=psf_cmap)
        ax[1].set_title('Jitter image'); ax[1].axis('off')

    Apply 10 um jitter to a 5x oversampled PSF. Note that because we are
    specifying jitter in terms of linear distance on the focal plane, we must
    provide the detector pixelscale when applying the Jitter. We also provide
    the oversampling factor, ensuring the convolution kernel is appropriately
    sized:

    .. plot::
        :context: close-figs
        :include-source:
        :nofigs:

        psf_jitter = lentil.jitter(psf, scale=10e-6, pixelscale=5e-6,
                                   oversample=5)

    .. plot::
        :context:
        :scale: 50

        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(4,2))
        ax[0].imshow(psf, cmap=psf_cmap)
        ax[0].set_title('Input image'); ax[0].axis('off')
        ax[1].imshow(psf_jitter, cmap=psf_cmap)
        ax[1].set_title('Jitter image'); ax[1].axis('off')

    References
    ----------
    [1] https://en.wikipedia.org/wiki/Convolution_theorem

    """
    img = np.asarray(img)
    x = np.fft.fftfreq(img.shape[1])
    y = np.fft.fftfreq(img.shape[0])
    xx, yy = np.meshgrid(x, y)
    rho = np.sqrt(xx ** 2 + yy ** 2)
    kernel = np.exp(-2 * (np.pi * (scale / pixelscale) * oversample * rho) ** 2)

    out = np.abs(np.fft.ifft2(np.fft.fft2(img)*kernel))
    return out * np.sum(img) / np.sum(out)  # rescale to preserve input weight


def smear(img, distance, angle=None, pixelscale=1, oversample=1):
    """
    Apply image smear via convolution.

    Parameters
    ----------
    img : array_like
        Input image
    distance : float
        Linear smear distance

        If ``pixelscale = 1``, ``distance`` should give smear in terms of
        fractional pixels. If ``pixelscale`` is the physical dimension of
        one pixel, ``distance`` should give smear in terms of the same
        units.
    angle : float, optional
        Smear direction in degrees measured clockwise from the x-axis. If
        ``None`` (default), a new direction is randomly chosen every time
        :class:`Smear` is called.
    pixelscale : float, optional
        Pixel size. If ``pixelscale = 1`` (default), the ``scale`` parameter
        provided when Smear is called should give jitter in terms of fractional
        pixels. If ``pixelscale`` is the physical dimension of one pixel, ``the
        ``scale`` parameter provided when Smear is called should give jitter in
        terms of the same units.
    oversample : int, optional
        Oversampling factor of img. Default is 1.

    Returns
    -------
    out : ndarray
        Image with smear applied.

    Examples
    --------
    Apply 25 pixels smear at 30 degrees to a 5x oversampled PSF:

    .. plot::
        :context: reset

        mask = lentil.circle((256, 256), 120)
        opd = lentil.zernike_compose(mask, [0, 0, 0, -300e-9, 50e-9, -100e-9, 50e-9])
        pupil = lentil.Pupil(amplitude=mask, opd=opd, focal_length=10, pixelscale=1 / 240)
        w = lentil.Wavefront(650e-9)
        w *= pupil
        w = lentil.propagate_dft(w, pixelscale=5e-6, shape=64, oversample=5)
        psf = w.intensity    

    .. plot::
        :context: close-figs
        :include-source:
        :nofigs:

        psf_smear = lentil.smear(psf, distance=25, angle=30,
                                 oversample=5)

    .. plot::
        :context:
        :scale: 50

        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(4,2))
        ax[0].imshow(psf, cmap=psf_cmap)
        ax[0].set_title('Input image'); ax[0].axis('off')
        ax[1].imshow(psf_smear, cmap=psf_cmap)
        ax[1].set_title('Smeared image'); ax[1].axis('off')

    Apply 50 um smear at a random angle to a 5x oversampled PSF. Note that
    because we are specifying smear in terms of linear distance on the focal
    plane, we must also provide the detector pixelscale:

    .. plot::
        :context: close-figs
        :include-source:
        :nofigs:

        psf_smear = lentil.smear(psf, distance=50e-6,
                                 pixelscale=5e-6, oversample=5)

    .. plot::
        :context:
        :scale: 50

        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(4,2))
        ax[0].imshow(psf, cmap=psf_cmap)
        ax[0].set_title('Input image'); ax[0].axis('off')
        ax[1].imshow(psf_smear, cmap=psf_cmap)
        ax[1].set_title('Smeared image'); ax[1].axis('off')

    References
    ----------
    [1] https://en.wikipedia.org/wiki/Convolution_theorem

    """
    img = np.asarray(img)
    x = np.fft.fftfreq(img.shape[1])
    y = np.fft.fftfreq(img.shape[0])

    xx, yy = np.meshgrid(x, y)

    if angle is None:
        angle = np.random.uniform(0, 2 * np.pi)
    else:
        angle = np.radians(angle)

    yy_rot = np.sin(angle) * yy + np.cos(angle) * xx

    kernel = np.sinc(yy_rot * (distance / pixelscale) * oversample)

    out = np.abs(np.fft.ifft2(np.fft.fft2(img)*kernel))
    return out * np.sum(img) / np.sum(out)  # rescale to preserve input weight

