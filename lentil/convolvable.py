import numpy as np

__all__ = ['Pixel', 'Jitter', 'Smear']


class Convolvable:
    """Base class for representing spatial artifacts that can be applied via
    multiplication with a kernel in the frequency domain [1].

    See Also
    --------
    * :class:`~lentil.convolvable.Jitter`
    * :class:`~lentil.convolvable.Smear`
    * :class:`~lentil.convolvable.Pixel`

    References
    ----------
    [1] https://en.wikipedia.org/wiki/Convolution_theorem

    """

    def __call__(self, *args, **kwargs):
        """Perform convolution

        This method defines the callable interface. Subclasses should provide a
        concrete implementation that computes an appropriate kernel and returns
        the result of :func:`convolve`.

        """
        raise NotImplementedError

    @staticmethod
    def convolve(f, G):
        r"""Convolve an array with a frequency domain kernel.

        Given an array :math:`f` in the spatial domain and a kernel :math:`G` in
        the frequency domain, the convolution is computed as

        .. math::

            f * g = \mathcal{F}^{-1} \left\{ \mathcal{F}\left\{f\right\}\cdot G \right\}


        Parameters
        ----------
        f : ndarray
            Spatial domain array

        G : ndarray
            Frequency domain kernel

        Returns
        -------
        out : ndarray
            f * g

        """
        return np.fft.ifft2(np.fft.fft2(f)*G)


class Pixel(Convolvable):
    """Callable object for applying the aperture effects of a square pixel on a
    discretely sampled irradiance distribution.

    Example
    -------
    Apply pixel MTF to a 3x oversampled PSF:

    .. code:: pycon

        >>> import lentil
        >>> import matplotlib.pyplot as plt
        >>> psf = ...  # PSF calculation details omitted
        >>> pixel_mtf = lentil.convolvable.Pixel()
        >>> psf_mtf = pixel_mtf(psf, oversample=3)
        >>> psf_detector = lentil.util.rescale(psf_mtf, 1/3)

    Note that both the pixel MTF and detector resampling operations preserve
    radiometry:

    .. code:: pycon

        >>> print(np.sum(psf), np.sum(psf_mtf), np.sum(psf_detector))
        50398.80556524441 50398.80556524441 50398.80556524441

    """

    def __call__(self, img, oversample=1):
        """Convolve ``img`` with a square pixel transfer function.

        Parameters
        ----------
        img : array_like
            Input image

        oversample : int, optional
            Oversampling factor of img. Default is 1.

        Returns
        -------
        out : ndarray
            Image with pixel sampling effects applied.

        """
        img = np.asarray(img)
        x = np.fft.fftfreq(img.shape[1])
        y = np.fft.fftfreq(img.shape[0])
        kernel = self.kernel(x, y, oversample)

        return np.abs(self.convolve(img, kernel))

    def kernel(self, x, y, oversample):
        mtf_x = np.sinc(x*oversample)
        mtf_y = np.sinc(y*oversample)
        return np.dot(mtf_x[:, np.newaxis], mtf_y[np.newaxis, :])


class Jitter(Convolvable):
    """Callable object for applying image jitter via convolution.

    Parameters
    ----------
    pixelscale : float, optional
        Pixel size. If ``pixelscale = 1`` (default), the ``scale`` parameter
        provided when Jitter is called should give jitter in terms of fractional
        pixels. If ``pixelscale`` is the physical dimension of one pixel, ``the
        ``scale`` parameter provided when Jitter is called should give jitter in
        terms of the same units.

    Examples
    --------
    Apply 2 pixels jitter to a natively-sampled PSF:

    .. code:: pycon

        >>> import lentil
        >>> import matplotlib.pyplot as plt
        >>> psf = ...  # PSF calculation details omitted
        >>> jitter = lentil.convolvable.Jitter()
        >>> psf_jitter = jitter(psf, scale=2)
        >>> plt.subplot(121), plt.imshow(psf)
        >>> plt.subplot(122), plt.imshow(psf_jitter)

    .. image:: ../_static/img/api/convolvable/jitter_px.png
        :width: 500px

    Apply 20 um jitter to a 3x oversampled PSF. Note that because we are
    specifying jitter in terms of linear distance on the focal plane, we must
    provide the detector pixelscale when creating the Jitter object. We also
    provide the oversampling factor, ensuring the convolution kernel is
    appropriately sized:

    .. code:: pycon

        >>> import lentil
        >>> import matplotlib.pyplot as plt
        >>> psf = ...  # PSF calculation details omitted
        >>> jitter = lentil.convolvable.Jitter(pixelscale=5e-6)
        >>> psf_jitter = jitter(psf, scale=20e-6, oversample=3)
        >>> plt.subplot(121), plt.imshow(psf)
        >>> plt.subplot(122), plt.imshow(psf_jitter)

    .. image:: ../_static/img/api/convolvable/jitter_m.png
        :width: 500px

    """

    def __init__(self, pixelscale=1):
        self.pixelscale = pixelscale

    def __call__(self, img, scale, oversample=1):
        """Convolve ``img`` with a Gaussian jitter transfer function.

        Parameters
        ----------
        img : array_like
            Input image

        scale : float
            1-sigma jitter motion

            If ``pixelscale = 1``, ``scale`` should give jitter in terms of
            fractional pixels. If ``pixelscale`` is the physical dimension of
            one pixel, ``scale`` should give jitter in terms of the same units.

        oversample : int, optional
            Oversampling factor of img. Default is 1.

        Returns
        -------
        out : ndarray
            Image with jitter applied.

        """
        img = np.asarray(img)
        x = np.fft.fftfreq(img.shape[1])
        y = np.fft.fftfreq(img.shape[0])
        kernel = self.kernel(x, y, scale, oversample)

        return np.abs(self.convolve(img, kernel))

    def kernel(self, x, y, scale, oversample):
        xx, yy = np.meshgrid(x, y)
        rho = np.sqrt(xx**2 + yy**2)
        return np.exp(-2*(np.pi*(scale/self.pixelscale)*oversample*rho)**2)


class Smear(Convolvable):
    """Callable object for applying image smear via convolution.

    Parameters
    ----------
    pixelscale : float, optional
        Pixel size. If ``pixelscale = 1`` (default), the ``scale`` parameter
        provided when Smear is called should give jitter in terms of fractional
        pixels. If ``pixelscale`` is the physical dimension of one pixel, ``the
        ``scale`` parameter provided when Smear is called should give jitter in
        terms of the same units.

    Examples
    --------
    Apply 100 um smear at 60 degrees to a natively-sampled PSF. Note that
    because we are specifying smear in terms of linear distance on the focal
    plane, we must provide the detector pixelscale when creating the Smear
    object:

    .. code:: pycon

        >>> import lentil
        >>> import matplotlib.pyplot as plt
        >>> psf = ...  # PSF calculation details omitted
        >>> smear = lentil.convolvable.Smear(pixelscale=5e-6)
        >>> psf_smear = smear(psf, distance=100e-6, angle=60)
        >>> plt.subplot(121), plt.imshow(psf)
        >>> plt.subplot(122), plt.imshow(psf_smear)

    .. image:: ../_static/img/api/convolvable/smear_m.png
        :width: 500px

    Apply 10 pixels smear at a random angle to a 3x oversampled PSF:

    .. code:: pycon

        >>> import lentil
        >>> import matplotlib.pyplot as plt
        >>> psf = ...  # PSF calculation details omitted
        >>> smear = lentil.convolvable.Smear()
        >>> psf_smear = smear(psf, distance=10, oversample=3)
        >>> plt.subplot(121), plt.imshow(psf)
        >>> plt.subplot(122), plt.imshow(psf_smear)

    .. image:: ../_static/img/api/convolvable/smear_px.png
        :width: 500px

    """

    def __init__(self, pixelscale=1):
        self.pixelscale = pixelscale

    def __call__(self, img, distance, angle=None, oversample=1):
        """Convolve ``img`` with a rotated sinc smear transfer function.

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

        oversample : int, optional
            Oversampling factor of img. Default is 1.

        Returns
        -------
        out : ndarray
            Image with smear applied.

        """
        img = np.asarray(img)
        x = np.fft.fftfreq(img.shape[1])
        y = np.fft.fftfreq(img.shape[0])
        kernel = self.kernel(x, y, distance, angle, oversample)

        return np.abs(self.convolve(img, kernel))

    def kernel(self, x, y, distance, angle, oversample):
        xx, yy = np.meshgrid(x, y)

        if angle is None:
            angle = np.random.uniform(0, 2*np.pi)
        else:
            angle = np.radians(angle)

        yy_rot = -np.sin(angle)*xx + np.cos(angle)*yy

        return np.sinc(yy_rot*(distance/self.pixelscale)*oversample)
