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

    .. code:: pycon

        >>> import lentil
        >>> import matplotlib.pyplot as plt
        >>> psf = ...  # PSF calculation details omitted
        >>> psf_jitter = lentil.jitter(psf, scale=2, oversample=5)
        >>> plt.subplot(121), plt.imshow(psf)
        >>> plt.subplot(122), plt.imshow(psf_jitter)

    ..plot:: _img/python/jitter.py
        :scale: 50

    Apply 10 um jitter to a 5x oversampled PSF. Note that because we are
    specifying jitter in terms of linear distance on the focal plane, we must
    provide the detector pixelscale when applying the Jitter. We also provide
    the oversampling factor, ensuring the convolution kernel is appropriately
    sized:

    .. code:: pycon

        >>> import lentil
        >>> import matplotlib.pyplot as plt
        >>> psf = ...  # PSF calculation details omitted
        >>> psf_jitter = lentil.jitter(psf,
        ...                            scale=10e-6,
        ...                            pixelscale=5e-6,
        ...                            oversample=5)
        >>> plt.subplot(121), plt.imshow(psf)
        >>> plt.subplot(122), plt.imshow(psf_jitter)

    .. plot:: _img/python/jitter.py
        :scale: 50

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

    .. code:: pycon

        >>> import matplotlib.pyplot as plt
        >>> import lentil
        >>> psf = ...  # PSF calculation details omitted
        >>> psf_smear = lentil.smear(psf, distance=25,
        ...                          angle=30)
        >>> plt.subplot(121), plt.imshow(psf)
        >>> plt.subplot(122), plt.imshow(psf_smear)

    .. plot:: _img/python/smear_directional.py
        :scale: 50

    Apply 50 um smear at a random angle to a 5x oversampled PSF. Note that
    because we are specifying smear in terms of linear distance on the focal
    plane, we must also provide the detector pixelscale:

    .. code:: pycon

        >>> import lentil
        >>> import matplotlib.pyplot as plt
        >>> psf = ...  # PSF calculation details omitted
        >>> psf_smear = lentil.smear(psf, distance=5e-5,
        ...                          pixelscale=5e-6, oversample=5)
        >>> plt.subplot(121), plt.imshow(psf)
        >>> plt.subplot(122), plt.imshow(psf_smear)

    .. plot:: _img/python/smear.py
        :scale: 50

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

