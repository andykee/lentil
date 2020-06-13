import warnings

import numpy as np
from scipy import signal

from lentil import util

__all__ = ['Gain', 'PolynomialGain', 'ShotNoise', 'GaussianShotNoise',
           'DarkCurrent', 'ReadNoise', 'Rule07DarkCurrent', 'FPN', 'RandomFPN',
           'NormalFPN', 'UniformFPN', 'LognormalFPN', 'PRNU', 'DarkCurrentFPN',
           'PixelOffsetFPN', 'ColumnOffsetFPN', 'ChargeDiffusion', 'Defect',
           'PixelMask', 'CosmicRay']


class Windowable:
    """Base class for representing data that can have a specific sub-window
    extracted.

    Parameters
    ----------
    data : float or array_like
        Windowable data.

    """

    def __init__(self, data):
        self.data = np.asarray(data)

    def window(self, shape=None, window=None):
        """Return an appropriately sized, potentially windowed :attr:`data`
        array.

        Parameters
        ----------
        shape : array_like or None, optional
            Output shape given as (nrows, ncols). If ``None`` (default), an
            uncropped array is returned.

        window : array_like or None, optional
            Indices of :attr:`data` array to return given as (r_start, r_end,
            c_start, c_end). This definition follows standard numpy indexing.

        Returns
        -------
        data : ndarray
            Trimmed and windowed :attr:`data` array

        Notes
        -----
        * If ``data`` is a single value (data.size == 1), self.data is returned
          regardless of what ``shape`` and ``window`` are.
        * If ``shape`` is given but ``window`` is ``None``, the returned ndarray
          is trimmed about the center of the array using
          :func:`lentil.util.pad`.
        * If ``window`` is given but ``shape`` is ``None``, the returned ndarray
          is extracted from :attr:`data` according to the indices in ``window``
        * If both ``shape`` and ``window`` are given, the returned ndarray is
          extracted from :attr:`data` according to the indices in ``window`` and
          the following expressions must be true:

        .. code::

            shape[0] = (window[1] - window[0]) = (r_end - r_start)
            shape[1] = (window[3] - window[2]) = (c_end - c_start)

        """
        if self.data.size == 1:
            return self.data

        if shape is None and window is None:
            # return the entire array
            return self.data

        elif window is not None:
            if shape is not None:
                # ensure consistency
                assert (window[1] - window[0]) == shape[0]
                assert (window[3] - window[2]) == shape[1]
            # return the windowed array.
            # Note that numpy implicitly handles a third dimension if one is
            # present
            return self.data[window[0]:window[1], window[2]:window[3]]

        else:  # shape is None
            # return the padded array
            # The same comment about numpy implicitly handling a third dimension
            # (see above) holds true here, although in this case the logic
            # resides in pad
            return util.pad(self.data, shape)


class Gain(Windowable):
    """Linear gain model.

    Parameters
    ----------
    gain : float or array_like
        Conversion gain in DN/e-. If specified as a single value, it is applied
        uniformly. If specified as an array, it is assumed to provide a
        pixel-by-pixel gain map.

    saturation_capacity : int or None
        Electron count resulting in pixel saturation. If None, pixels will not
        saturate. This is obviously nonphysical, but can be useful for testing
        or debugging.

    dtype : data-type or None, optional
        Output data-type. If None (default), no data-type casting is performed.

    Note
    ----
    The saturation capacity should not be confused with the full-well capacity.
    Saturation capacity is typically smaller than the full well capacity because
    the signal is clipped before the physical saturation of the pixel is
    reached.
    """

    def __init__(self, gain, saturation_capacity, dtype=None):
        super().__init__(gain)
        self.saturation_capacity = saturation_capacity
        self.dtype = dtype

    def __call__(self, img, window=None, warn_saturate=False):
        """Apply gain model to an array

        Parameters
        ---------
        img : ndarray
            Array of electron counts

        window : array_like or None, optional
            Indices of gain to return given as (r_start, r_end, c_start, c_end).
            This definition follows standard Numpy indexing.

        warn_saturate : bool, optional
            Raise a warning when pixels saturate. Default is False.

        Returns
        -------
        img : ndarray
            Array of DN

        """
        # enforce the saturation limit
        if self.saturation_capacity:

            if warn_saturate:
                if self.check_saturate(img):
                    warnings.warn('Frame has saturated pixels.')

            # Apply the saturation limit
            img[img > self.saturation_capacity] = self.saturation_capacity

        # convert to DN
        return np.floor(img * self.gain(img.shape, window), dtype=self.dtype)

    def gain(self, shape=None, window=None):
        """Return a trimmed and/or windowed representation of the gain model."""
        return self.window(shape, window)

    def check_saturate(self, img):
        if np.any(img > self.saturation_capacity):
            return True
        else:
            return False


class PolynomialGain(Gain):
    """Gain model defined by polynomial coefficients.

    Parameters
    ----------
    gain : array_like
        Conversion gain in DN/e-. Can be specified in two ways:

            * As a one-dimensional array of polynomial coefficients applied
              globally to each pixel
            * As a three-dimensional array of pixel-by-pixel gain where the
              first dimension gives polynomial coefficients of each pixel

    saturation_capacity : int or None
        Electron count resulting in pixel saturation. If None, pixels will not
        saturate. This is obviously nonphysical, but can be useful for testing
        or debugging.

    dtype : data-type or None, optional
        Output data-type. If None (default), no data-type casting is performed.

    Note
    ----
    The saturation capacity should not be confused with the full-well capacity.
    Saturation capacity is typically smaller than the full well capacity because
    the signal is clipped before the physical saturation of the pixel is
    reached.

    """

    def __init__(self, gain, saturation_capacity, dtype=None):

        gain = np.asarray(gain)
        if not np.any([gain.ndim == 1, gain.ndim == 3]):
            raise ValueError

        super().__init__(gain, saturation_capacity, dtype)

    def __call__(self, img, window=None, warn_saturate=False):
        """Apply gain model to an array

        Parameters
        ---------
        img : ndarray
            Array of electron counts

        window : array_like or None, optional
            Indices of gain to return given as (r_start, r_end, c_start, c_end).
            This definition follows standard Numpy indexing.

        warn_saturate : bool, optional
            Raise a warning when pixels saturate. Default is False.

        Returns
        -------
        img : ndarray
            Array of DN

        """
        # enforce the saturation limit
        if self.saturation_capacity:

            if warn_saturate:
                if self.check_saturate(img):
                    warnings.warn('Frame has saturated pixels.')

            # Apply the saturation limit
            img[img > self.saturation_capacity] = self.saturation_capacity

        # Determine the polynomial order
        if self.data.ndim == 1:
            model_order = self.data.shape[0]
        else:
            model_order = self.data.shape[2]

        # Prepare a cube of electron counts to apply the polynomial gain to
        img_cube = np.repeat(img[:, :, np.newaxis], model_order-1, axis=2)
        for order in np.arange(model_order-1, 0, -1):
            d = model_order - order - 1
            img_cube[..., d] = img_cube[..., d]**order

        # Get the appropriately sized and windowed gain model
        gain_model = self.gain(img.shape, window)

        # Apply the gain model and convert to DN
        if self.data.ndim == 1:
            gain, offset = gain_model[0:-1], gain_model[-1]
            return np.einsum('ijk,k->ij', img_cube, gain) + offset
        else:
            # self.data.ndim == 3
            gain, offset = gain_model[..., 0:-1], gain_model[..., -1]
            return np.einsum('ijk,ijk->ij', img_cube, gain) + offset

    def gain(self, shape=None, window=None):
        """Return a trimmed and/or windowed representation of the gain model."""
        if self.data.ndim == 1:
            # we have a global nonlinear gain where the shape and window don't
            # matter
            return self.data
        else:
            # self.data.ndim == 3
            return self.window(shape, window)


class RandomNoise:
    """Base class for all radiometry of random noise with an optional
    user-definable seed.

    If a seed is provided, an instance of ``numpy.random.RandomState`` is
    created and all random numbers are generated from this object in a
    repeatable (but still random) way. If no seed is provided, random numbers
    are generated using the appropriate ``numpy.random`` distribution.

    Parameters
    ----------
    seed : None, int, or array_like, optional
        Random seed used to initialize ``numpy.random.RandomState``. If
        ``None``, no random number generator object is created.

    Attributes
    ----------
    rng : numpy.random.RandomState or None
        Random number generator object. If ``seed`` is ``None``, ``rng`` will
        also be ``None``
    """

    def __init__(self, seed=None):
        self.seed = seed

        if seed:
            self.rng = np.random.RandomState(seed)
        else:
            self.rng = None


class ShotNoise(RandomNoise):
    r"""Shot noise modeled by a Poisson process.

    Parameters
    ----------
    seed : None, int, or array_like, optional
        Seed for random number generator. If ``None`` (default), the Numpy
        default is used.

    Note
    ----
    For sufficiently large values of :math:`\lambda`,  (say :math:`\lambda >
    1000`), the Normal(:math:`\mu = \lambda`, :math:`\sigma^2 = \lambda`)
    distribution is an excellent approximation to the Poisson(:math:`\lambda`)
    distribution and is about 25% faster to compute.

    See Also
    --------
    :class:`~lentil.detector.GaussianShotNoise` Shot noise modeled by a normal
    distribution

    """

    def __init__(self, seed=None):
        super().__init__(seed)

    def __call__(self, img):
        """Apply shot noise model to an array

        Parameters
        ---------
        img : ndarray
            Array of electron counts

        Returns
        -------
        img : ndarray
            Array of noisy electron counts

        """
        if self.rng:
            return self.rng.poisson(img)
        else:
            return np.random.poisson(img)


class GaussianShotNoise(RandomNoise):
    r"""Shot noise modeled by a normal distribution.

    Parameters
    ----------
    seed : None, int, or array_like, optional
        Seed for random number generator. If ``None`` (default), the Numpy
        default is used.

    Note
    ----
    For sufficiently large values of :math:`\lambda`,  (say :math:`\lambda >
    1000`), the Normal(:math:`\mu = \lambda`, :math:`\sigma^2 = \lambda`)
    distribution is an excellent approximation to the Poisson(:math:`\lambda`)
    distribution and is about 25% faster to compute.

    See Also
    --------
    :class:`~lentil.detector.ShotNoise` Shot noise modeled by a Poisson process

    """

    def __init__(self, seed=None):
        super().__init__(seed)

    def __call__(self, img):
        """Apply shot noise model to an array

        Parameters
        ---------
        img : ndarray
            Array of electron counts

        Returns
        -------
        img : ndarray
            Array of noisy electron counts

        """
        if self.rng:
            return self.rng.normal(loc=img, scale=img)
        else:
            return np.random.normal(loc=img, scale=img)


class ReadNoise(RandomNoise):
    """Read noise model.

    Parameters
    ----------
    electrons : int
        Read noise per frame

    seed : None, int, or array_like, optional
        Seed for random number generator. If ``None`` (default), the Numpy
        default is used.

    """
    def __init__(self, electrons, seed=None):
        super().__init__(seed)
        self.electrons = electrons

    def __call__(self, frame):
        if self.rng:
            noise = self.rng.normal(loc=0.0, scale=self.electrons, size=frame.shape)
        else:
            noise = np.random.normal(loc=0.0, scale=self.electrons, size=frame.shape)
        return frame + noise


class DarkSignal(Windowable):
    """Base class for representing a dark signal.

    Parameters
    ----------
    value : float, or array_like
        Dark signal value. If specified as a single value, it is applied
        uniformly. If specified as an array, it is assumed to provide a
        pixel-by-pixel map.

    """
    def __init__(self, value):
        super().__init__(value)

    def __call__(self, shape, window=None):
        """Create dark signal array

        Parameters
        ---------
        shape : array_like
            Dimensions of dark signal array to create

        window : array_like or None, optional
            Indices of dark signal to return given as (r_start, r_end, c_start,
            c_end). This definition follows standard Numpy indexing.

        Returns
        -------
        dark_signal : ndarray
            Dark signal array

        """
        # Return ones * the dark() value in the event that a single value is
        # provided (and thus returned by dark)
        dark_signal = self.dark(shape, window) * np.ones(shape)

        return dark_signal

    def dark(self, shape, window):
        return self.window(shape, window)


class DarkCurrent(DarkSignal):
    """Generic dark current model.

    Parameters
    ----------
    rate : float or array_like
        Dark current rate in electrons/sec/px. If specified as a single value,
        it is applied uniformly. If specified as an array, it is assumed to
        provide a pixel-by-pixel map.

    See Also
    --------
    :class:`~lentil.detector.Rule07DarkCurrent` Rule 07 dark current model for
    IR detectors

    """

    def __init__(self, rate):
        super().__init__(rate)

    @property
    def rate(self):
        return self.data


class Rule07DarkCurrent(DarkCurrent):
    """Model for representing dark current in HgCdTe infrared detectors using
    Rule 07.

    Parameters
    ----------
    temperature : float
        Focal plane temperature in K

    cutoff_wavelength : float
        Focal plane cutoff wavelength in m

    pixelscale : float
        Size of one pixel in m

    References
    ----------
    * Tennant, W.E. et al. MBE HgCdTe Technology: A Very General Solution to IR
      Detection, Described by "Rule 07", a Very Convenient Heuristic. Journal of
      Electronic Materials (2008).

    """

    def __init__(self, temperature, cutoff_wavelength, pixelscale):

        self.temperature = temperature
        self.cutoff_wavelength = cutoff_wavelength

        J0 = 8367.00001853855  # A/cm^2
        C = -1.16239134096245
        k = 1.3802e-23  # J/K - Boltzmanns' constant
        q = 1.6021e-19  # C - electron charge

        lambda_threshold = 4.63513642316149
        lambda_scale = 0.200847413564122
        P = 0.544071281108481

        # compute lambda_e
        lambda_cutoff = cutoff_wavelength * 1e6  # meters to microns
        if lambda_cutoff >= lambda_threshold:
            lambda_e = lambda_cutoff
        else:
            lambda_e = lambda_cutoff/(1-((lambda_scale/lambda_cutoff) -
                                         (lambda_scale/lambda_threshold))**P)

        # apply Rule 07 to compute the dark flux in A/cm^2
        J = J0*np.exp(C*(1.24*q/(k*lambda_e*temperature)))

        # convert J in A/cm^2 to J in e-/px/s
        px_area = (pixelscale*1e2)**2  # pixel area in cm^2
        super().__init__(rate=1/q * px_area * J)

    def __call__(self, shape, *args, **kwargs):
        """Create dark signal array

        Parameters
        ---------
        shape : array_like
            Dimensions of dark signal array to create

        Returns
        -------
        dark_signal : ndarray
            Dark signal array

        """
        return self.rate * np.ones(shape)


class FPN(Windowable):
    """Base class for all radiometry of fixed-pattern noise .

    Parameters
    ----------
    data : array_like
        FPN array.

    """
    def __init__(self, data):
        super().__init__(data)

    def __call__(self, frame, window=None):
        """Apply fixed pattern noise to an array.

            The output array :math:`A'` is related to the input array :math:`A`
            by

            .. math::

                A' = A*F

            where :math:`F` is the :attr:`data` FPN array.

        Parameters
        ----------
        frame : ndarray

        window : array_like or None, optional
            Indices of :attr:`data` FPN array to return given as (r_start,
            r_end, c_start, c_end). This definition follows standard Numpy
            indexing.

        Returns
        -------
        frame : ndarray
            Array with  FPN applied

        """
        return frame * self.fpn(frame.shape, window)

    def fpn(self, shape=None, window=None):
        return self.window(shape, window)


class RandomFPN(FPN):
    """Base class for randomly generated fixed-pattern noise with an optional
    user-definable seed.

    Parameters
    ----------
    seed : None, int, or array_like, optional
        Seed for random number generator. If ``None`` (default), the Numpy
        default is used.

    Attributes
    ----------
    rng : numpy.random.RandomState or None
        Random number generator object

    data : ndarray
        FPN array.

    """
    def __init__(self, seed=None):
        self.seed = seed
        if seed:
            self.rng = np.random.RandomState(seed)
        else:
            self.rng = None

        self.data = None


class NormalFPN(RandomFPN):
    """Base class for fixed pattern noise represented by a normal distribution.

    Parameters
    ----------
    loc : float

    scale : float

    size : float

    """
    def __init__(self, loc, scale, size, seed=None):
        super().__init__(seed)
        if self.rng:
            self.data = self.rng.normal(loc, scale, size)
        else:
            self.data = np.random.normal(loc, scale, size)


class LognormalFPN(RandomFPN):
    """Base class for fixed pattern noise represented by a log-normal
    distribution.

    Parameters
    ----------
    mean : float

    sigma : float

    size : float

    """
    def __init__(self, mean, sigma, size, seed=None):
        super().__init__(seed)
        if self.rng:
            self.data = self.rng.lognormal(mean, sigma, size)
        else:
            self.data = np.random.lognormal(mean, sigma, size)


class UniformFPN(RandomFPN):
    """Base class for fixed pattern noise represented by a uniform distribution.

    Parameters
    ----------
    low : float

    high : float

    size : float

    """
    def __init__(self, low, high, size, seed=None):
        super().__init__(seed)
        if self.rng:
            self.data = self.rng.uniform(low, high, size)
        else:
            self.data = np.random.uniform(low, high, size)


class PRNU(NormalFPN):
    """Model for representing photoresponse non-uniformity.

    """
    def __init__(self, prnu_factor, size, seed=None):
        self.prnu_factor = prnu_factor
        super().__init__(loc=1.0, scale=prnu_factor, size=size, seed=seed)


class DarkCurrentFPN(LognormalFPN):
    """Dark current FPN model, represented by a log-normal distribution.

    Parameters
    ----------
    sigma : float
        Dark current FPN factor, usually between 0.1 and 0.4 for CCD and CMOS
        detectors [2].

    size : None, tuple, or array_like, optional
        Precomputed FPN matrix size. Typically this will be equal to the
        detector size. See :class:`GainFPN` note for details on the two
        different modes for when ``size`` is specified or not.

    seed : None, int, or array_like, optional
        Random seed used to initialize ``numpy.random.RandomState``. If
        ``None``, then ``RandomState`` will try to read data from /dev/urandom
        (or the Windows analogue) if available or seed from the clock otherwise.

    References
    ----------
    [1] `Log-normal distribution - Wikipedia <https://en.wikipedia.org/wiki/Log-normal_distribution>`_

    [2] Janesick, J. R., Photon Transfer, Vol. PM170, SPIE (2007)

    """

    def __init__(self, sigma, size, seed=None):
        super().__init__(mean=1.0, sigma=sigma, size=size, seed=seed)
        self.sigma = sigma


class PixelOffsetFPN(NormalFPN):
    r"""Pixel offset FPN model.

    Pixel offset FPN can be represented by an autoregressive random process with
    (:math:`\mu = 0.0, \sigma = ` :attr:`sigma`) [1], [2].

    Parameters
    ----------
    sigma : float
        1-sigma pixel offset FPN value in number of electrons

    ar_parameter : float
        Autoregressive parameter which characterizes the dependency of each
        pixel on its neighbors. Valid values are between 0 and 0.5 [1]. Small
        values correspond to less correlation.

    size : None or array_like, optional
        Precomputed FPN matrix size. Typically this will be equal to the
        detector size. See :class:`GainFPN` note for details on the two
        different modes for when ``size`` is specified or not.

    seed : None, int, or array_like, optional
        Random seed used to initialize ``numpy.random.RandomState``. If
        ``None``, then ``RandomState`` will try to read data from /dev/urandom
        (or the Windows analogue) if available or seed from the clock otherwise.

    References
    ----------
    [1] El Gamal, A. et. al., Modeling and Estimation of FPN Components in CMOS Image Sensors, SPIE Proc. 3301 (1998)

    [2] `Autoregressive model - Wikipedia <https://en.wikipedia.org/wiki/Autoregressive_model>`_

    """

    def __init__(self, sigma, ar_parameter, size, seed=None):
        assert (ar_parameter >= 0) and (ar_parameter <= 0.5)

        super().__init__(loc=1.0, scale=sigma, size=size, seed=seed)
        self.sigma = sigma
        self.ar_parameter = ar_parameter

        # apply the autoregressive filter
        self.data = signal.lfilter([1], [1, self.ar_parameter], self.data, axis=0)


class ColumnOffsetFPN(NormalFPN):
    r"""Column offset FPN model.

    Column offset FPN can be represented by an autoregressive random process
    with (:math:`\mu = 0.0, \sigma = ` :attr:`sigma`) [1], [2].

    Parameters
    ----------
    sigma : float
        1-sigma column offset FPN value in number of electrons

    ar_parameter : float
        Autoregressive parameter which characterizes the dependency of each
        column on its neighbors. Valid values are between 0 and 0.5 [1]. Small
        values correspond to less correlation.

    size : None or array_like, optional
        Precomputed FPN matrix size. Typically this will be equal to the
        detector size. See :class:`GainFPN` note for details on the two
        different modes for when ``size`` is specified or not.

    seed : None, int, or array_like, optional
        Random seed used to initialize ``numpy.random.RandomState``. If
        ``None``, then ``RandomState`` will try to read data from /dev/urandom
        (or the Windows analogue) if available or seed from the clock otherwise.

    References
    ----------
    [1] El Gamal, A. et. al., Modeling and Estimation of FPN Components in CMOS
        Image Sensors, SPIE Proc. 3301 (1998)

    [2] `Autoregressive model - Wikipedia <https://en.wikipedia.org/wiki/Autoregressive_model>`_

    """

    def __init__(self, sigma, ar_parameter, size, seed=None):
        assert (ar_parameter >= 0) and (ar_parameter <= 0.5)

        super().__init__(loc=1.0, scale=sigma, size=size[1], seed=seed)
        self.sigma = sigma
        self.ar_parameter = ar_parameter

        # apply the autoregressive filter
        self.data = signal.lfilter([1], [1, self.ar_parameter], self.data)
        self.data = np.tile(self.data, (size[0], 1))


class ChargeDiffusion:
    """Charge diffusion model represented by a gaussian blur."""
    def __init__(self, sigma):
        self.sigma = sigma
        self.kernel = util.gaussian2d(3, self.sigma)

    def __call__(self, img):
        return signal.convolve2d(img, self.kernel, mode='same')


class Defect(Windowable):
    """Base class for representing detector defects."""
    pass


class PixelMask(Defect):
    """Class for representing defective, nonlinear, or dead pixels by providing
    a mask.
    """
    pass


class CosmicRay:
    """Cosmic ray generator for simulating cosmic ray hits on a detector.

    The default values for cosmic ray rates and electron fluxes are taken from
    [1]. It is commonly accepted at about 90% of cosmic rays are protons and the
    majority of the remaining 10% are alpha particles [1], [2].

    Parameters
    ----------
    rate : float, optional
        Cosmic ray rate expressed as number of events per square meter of
        detector area. Default is 40000 hits/m^2.

    proton_flux : float, optional
        Number of electrons liberated per meter of travel for a proton. Default
        is 1e9 e-/m.

    alpha_flux : float, optional
        Number of electrons liberated per meter of travel for an alpha particle.
        By definition, alpha particles have four times the energy of protons.
        Default is 4e9 e-/m.

    Example
    -------
    Simulate the cosmic ray hits for a 300 second exposure over a 256 x 256 detector
    patch with 5 um x 5 um x 3 um pixels:

    .. code:: pycon

        >>> import matplotlib.pyplot as plt
        >>> import lentil
        >>> cosmic_generator = lentil.detector.CosmicRay()
        >>> cosmic_frame = cosmic_generator((256,256), (5e-6, 5e-6, 3e-6), 300)
        >>> plt.imshow(cosmic_frame)

    .. image:: ../_static/img/api/detector/cosmic_ray.png
        :width: 300px

    References
    ----------
    [1] Offenberg, J.D. et. al.  Multi-Read Out Data Simulator. (2000).

    [2] `Cosmic ray - Wikipedia <https://en.wikipedia.org/wiki/Cosmic_ray>`_
    """

    def __init__(self, rate=4e4, proton_flux=1e9, alpha_flux=4e9):
        self.rate = rate
        self.proton_flux = proton_flux
        self.alpha_flux = alpha_flux

    def __call__(self, shape, pixelscale, ts):
        """Simulate cosmic ray hits on a detector.

        Parameters
        ----------
        shape : array_like
            Frame size defined as (rows, cols)

        pixelscale : array_like
            Pixel dimensions as (y,x,z) where z is the pixel depth

        ts : float
            Integration time in seconds

        Returns
        -------
        cosmic_e : ndarray
            Array representing the number of electrons in a detector image due
            to cosmic ray hits


        """

        # compute the number of rays that strike the detector during the
        # integration time
        nrays = self.nrays(shape, pixelscale, ts)

        img = np.zeros(shape)
        for ray in np.arange(0, nrays):
            img += self.cosmic_ray(shape, pixelscale)
        return img

    def nrays(self, shape, pixelscale, ts):
        area = shape[0]*pixelscale[0]*shape[1]*pixelscale[1]

        nrays = area * self.rate * ts

        # even if we don't have any events as defined by the rate, there is
        # still a chance we will see an event
        if nrays < 1:
            if np.random.uniform() <= nrays:
                nrays = 1
            else:
                nrays = 0

        return int(nrays)

    def cosmic_ray(self, shape, pixelscale):

        # randomly choose which type of particle we have
        if np.random.uniform() > 0.9:
            # alpha particle
            electron_flux = self.alpha_flux
        else:
            # proton
            electron_flux = self.proton_flux

        img = np.zeros(shape)

        # random position
        r = np.random.rand() * (shape[0]-1)
        c = np.random.rand() * (shape[1]-1)
        z = 0
        position = np.array([r, c, z])

        # give the cosmic ray a random direction (theta) and angle of incidence
        # (phi) and compute a unit vector representing the ray's direction of
        # travel. we can make things slightly easier on ourselves by limiting
        # phi to the interval 0 < phi < pi. this is assumed to be valid because
        # a cosmic ray striking the detector from behind will have essentially
        # the same effect as one striking it from the front. phi is negative to
        # force the z-direction of v to be negative (traveling into the detector
        # from above)
        theta = np.random.uniform() * 2 * np.pi
        phi = np.random.uniform() * 1*np.pi
        direction = np.array([np.cos(theta)*np.cos(phi), np.sin(theta)*np.cos(phi),
                              -np.sin(phi)])

        # scale the direction vector so that we can use ray_box_exit() assuming
        # a cube instead of a potentially rectangular pixel size
        pixelscale = np.asarray(pixelscale)
        scale = pixelscale/np.max(pixelscale)
        direction /= scale
        direction /= np.linalg.norm(direction)

        extent = (0, shape[0]-1, 0, shape[1]-1, 0, -1)

        # propagate the ray
        ray = self.propagate(position, direction, extent)

        for i in np.arange(0, ray.shape[0]-1):

            # compute the distance the ray traveled in pixel space
            dr = (ray[i+1][0]-ray[i][0])*pixelscale[0]
            dc = (ray[i+1][1]-ray[i][1])*pixelscale[1]
            dz = (ray[i+1][2]-ray[i][2])*pixelscale[2]
            dist = np.sqrt(dr**2+dc**2+dz**2)

            row = int(np.floor(ray[i+1][0]))
            col = int(np.floor(ray[i+1][1]))
            img[row, col] += electron_flux*dist

        return img

    @staticmethod
    def propagate(position, direction, extent):
        position = position[np.newaxis, ...]
        direction = direction[np.newaxis, ...]

        intersections = cubeplane_ray_intersection(position, direction, extent)
        raw_ray_data = process_cube_intersections(intersections, position, direction)

        return raw_ray_data[3]


def cubeplane_ray_intersection(xyz, dxdydz, extent):
    """Intersects rays with a plane of cubes and returns data about
    intersections

    Suppose the walls of the cube are on integer boundaries in the space of
    extent[0] <= x <= extent[1], extent[2] <= y <= extent[3], extent[4] <= z <=
    extent[5]

    Parameters
    ----------
    xyz : ndarray, ndim=2
        xyz[:, 0] is the x coordinate of the starting point of the ray,
        xyz[:, 1] is the y coordinate of the starting point of the ray.
        xyz[:, 2] is the z coordinate of the starting point.

    dxdydz : double ndim=2
        Unit vector representing the direction of travel of the ray.

    extent : tuple of int
        Integer address for the detector edge on the -x, +x, -y,
        +y, -z and +z edges. Should be representable as an int16.

    Returns
    -------
    inter : structured ndarray
        This is an array of structures with the dtype ('raylength', 'f4'),
        ('indx', 'u8'), ('facedirection', 'u1') where indx is the number of the
        row for the ray in xyz/dxdydz, raylength is the distance along the ray
        of the intersection and facedirection is an indicator of the ray
        direction. If face_direction = 0, 2, 4, the ray passed through the
        constant X, Y or Z planes, respectively, going from positive to
        negative. If face_direction = 1, 3, 5, the ray was going from negative
        to positive X, Y or Z through a surface.

    Notes
    -----
    If you don't like the output of this function, feed it into
    process_cube_intersections().

    This code has been tested when the rays start outside the cube volume. Would
    recommend not starting rays inside cube or within 2e-7 of a cube boundary.

    Very very very slightly adapted from code kindly developed by Dustin Moore.

    """

    xyz = xyz.astype(np.float32)
    dxdydz = dxdydz.astype(np.float32)

    indx = np.array(range(xyz.shape[0]), dtype=np.uint64)

    extent_shaped = np.array(extent, dtype=np.float32).reshape((1, 3, 2))

    # CUBE INTERSECTION
    # distance along ray before it penetrates outer CCD grid on xyz, -/+ side
    max_Ls = np.true_divide(extent_shaped - xyz[:, :, None], dxdydz[:, :, None])
    max_Ls[np.isinf(max_Ls)] = np.nan

    # Find length of ray to point where it enters cube array
    min_L = np.nanmax(np.min(max_Ls, axis=2), axis=1)

    # Find total length of ray to point as it exits cube array
    max_L = np.nanmin(np.max(max_Ls, axis=2), axis=1)
    del max_Ls

    # Find rays that never intersect cube or only intersect if the ray traced
    # backward
    no_intersection_mask = np.logical_or(max_L < min_L, max_L <= 0.0)

    # ... and chuck them from the calculation
    if np.any(no_intersection_mask):
        keep_mask = np.invert(no_intersection_mask)
        xyz = xyz[keep_mask, :]
        dxdydz = dxdydz[keep_mask, :]
        indx = indx[keep_mask]
        min_L = min_L[keep_mask]
        max_L = max_L[keep_mask]
        del keep_mask
    del no_intersection_mask

    # CUBE EDGE INTERSECTIONS
    # Shortest travel along each axis before ray exits cube array bounding
    # extent
    min_L[min_L < 0.0] = 0.0  # Rays can not go backward
    nearest = xyz + (min_L[:, None] - 1e-5) * dxdydz
    del min_L
    # Furthest travel along each axis before ray exits cube array bounding
    # extent
    furthest = xyz + (max_L[:, None] + 1e-5) * dxdydz
    del max_L

    d_pos = dxdydz >= 0.0
    least = np.ceil(np.where(d_pos, nearest, furthest)).astype(np.int16)
    great = np.floor(np.where(d_pos, furthest, nearest)).astype(np.int16) + 1

    d_pos = d_pos.astype(np.uint8)
    d_pos[:, 1] += 2
    d_pos[:, 2] += 4

    # No iterate when intersection impossible
    np.putmask(great, dxdydz == 0.0, -32768)

    intersections = []  # O(1) appends for the win.

    for xy, dxd, lea, grea, id, dp in zip(xyz, dxdydz, least, great, indx, d_pos):
        # algorithm-limiting performance here
        this_ints = [((q - xy[0]) / dxd[0], id, dp[0]) for q in range(lea[0], grea[0])]

        this_ints += [((q - xy[1]) / dxd[1], id, dp[1]) for q in range(lea[1], grea[1])]

        this_ints += [((q - xy[2]) / dxd[2], id, dp[2]) for q in range(lea[2], grea[2])]

        intersections += sorted(this_ints)

    # We repack the list into a numpy object for ease of conveyence back to the
    # mother process and later numpy computations
    fancy_dtype = np.dtype([('raylength', 'f4'), ('indx', 'u8'),
                            ('facedirection', 'u1')], align=True)
    inter = np.array(intersections, dtype=fancy_dtype)

    return inter


def process_cube_intersections(inter, xyz, dxdydz, verbose=False):
    """Intersects rays with a plane of cubes and returns data about
    intersections

    Suppose the walls of the cube are on integer boundaries in the space of
    extent[0] <= x <= extent[1], extent[2] <= y <= extent[3], extent[4] <= z
    <= extent[5]

    Parameters
    ----------
    inter : list of tuple
        Output of cubeplane_ray_intersection() or _cubeplane_ray_intersection()

    xyz : ndarray, ndim=2
        xyz[:, 0] is the x coordinate of the starting point of the ray,
        xyz[:, 1] is the y coordinate of the starting point of the ray.
        xyz[:, 2] is the z coordinate of the starting point.

    dxdydz : double ndim=2
        Unit vector representing the direction of travel of the ray.

    verbose: bool, optional
        If True, prints every ray intersection in human readable format. Default
        is False.

    Returns
    -------
    indx : ndarray of np.uint64
        Index of the ray for each ray intersection

    raylength : ndarray of np.float32
        Length along each ray that the intersection happened

    draylength : ndarray of np.float32
        Length along ray travsersed since ray was born or last intersection,
        which ever came most recently

    interpoint : ndarray of np.float32
        One row for each ray intersection, values are x,y,z of intersection
        point

    last_coord : ndarray of np.int16
        One row for each ray intersection, values are the x,y,z coordinate of
        the cube that the ray is leaving. (May not be inside extent, you should
        check if you care.)

    next_coord : ndarray of np.int16
        One row for each ray intersection, values are the x,y,z coordinate of
        the cube that the ray is entering. (May not be inside extent, you should
        check if you care.)

    Notes
    -----
    This code has been tested when the rays start outside the cube
    volume. Would recommend not starting rays inside cube or within 2e-7 of a
    cube boundary.

    Originally developed by Dustin Moore.

    """

    _xyz = xyz[inter['indx'], :].astype(np.float32)
    _dxdydz = dxdydz[inter['indx'], :].astype(np.float32)

    draylength = np.ediff1d(inter['raylength'], to_begin=(0.0,))
    np.putmask(draylength, np.ediff1d(inter['indx'], to_begin=(1,)) != 0,
               inter['raylength'])

    inter_point = _xyz + inter['raylength'][:, None] * _dxdydz

    next_pixel_bias = np.array([[-1e-2, 0.0, 0.0], [1e-2, 0.0, 0.0],
                                [0.0, -1e-2, 0.0], [0.0, 1e-2, 0.0],
                                [0.0, 0.0, -1e-2],  [0.0, 0.0, 1e-2]])

    w = np.choose(inter['facedirection'][:, None], next_pixel_bias)
    next_coord = np.empty(_xyz.shape, dtype=np.int16)
    np.floor(inter_point + w, out=next_coord, casting='unsafe')

    last_coord = np.empty(_xyz.shape, dtype=np.int16)
    np.floor(inter_point - w, out=last_coord, casting='unsafe')

    if verbose:
        last_indx = -1
        for ind, rl, ip, lc, nc, dl in zip(inter['indx'], inter['raylength'],
                                           inter_point, last_coord, next_coord,
                                           draylength):
            if last_indx != ind:
                print('Ray ' + str(ind) + ' is born at ' + str(xyz[ind, :]) +
                      ' and heads off in direction ' + str(dxdydz[ind, :]))
            print('  then leaves cube ' + str(lc) + ' at ' + str(ip) +
                  ' for cube ' + str(nc) + ' after having gone  ' +
                  str(dl) + ' and traveling ' + str(rl) + ' from inception')
            last_indx = ind

    return (inter['indx'], inter['raylength'], draylength, inter_point,
            last_coord, next_coord)
