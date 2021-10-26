import warnings

import numpy as np
import scipy.signal
import scipy.ndimage

import lentil

def collect_charge(img, wave, qe, waveunit='nm'):
    """Convert photon count (or flux) to electron count (or flux) by
    applying the detector's wavelength-dependent quantum efficiency.

    Parameters
    ----------
    img : array_like
        The photons presented to the sensor. Should have shape (nwave,
        nrows, ncols)
    wave : array_like
        Wavelengths corresponding to each slice in ``count``. The length of
        ``wave`` must be equal to the number of samples ``nwave`` in ``count``.
    qe : {:class:`lentil.radiometry.Spectrum` array_like, or scalar}
        Quantum efficiency used to convert detected photons to electrons.
    waveunit : str, optional
        Units of ``wave``. Defaults is ``nm``

    Returns
    -------
    img : ndarray
        Electron count or flux

    Notes
    -----
    The units of ``count`` don't really matter, as long as the user is aware
    that this method converts photons per whatever to electrons per
    whatever. Whatever is nothing for counts and seconds for flux.

    """

    img = np.asarray(img)
    if img.ndim == 2:
        img = img[np.newaxis, ...]

    qe = qe_asarray(qe, wave, waveunit)

    return np.einsum('ijk,i->jk', img, qe)


def collect_charge_bayer(img, wave, qe_red, qe_green, qe_blue, bayer_pattern, oversample=1, waveunit='nm'):
    """Convert photon count (or flux) to electron count (or flux) by
    applying the detector's wavelength-dependent quantum efficiency.
    Additional processing is performed to apply separate QE curves and
    location masks for the separate red, green, and blue channels.

    Parameters
    ----------
    img : array_like
        The photons presented to the sensor. Should have shape (nwave,
        nrows, ncols)
    wave : array_like
        Wavelengths corresponding to each slice in ``count``. The length of
        ``wave`` must be equal to the first dimension in ``count``.
    qe_red : {:class:`lentil.radiometry.Spectrum`, array_like, scalar}
        Red channel quantum efficiency
    qe_green : {:class:`lentil.radiometry.Spectrum`, array_like, scalar}
        Green channel quantum efficiency
    qe_blue : {:class:`lentil.radiometry.Spectrum`, array_like, scalar}
        Blue channel quantum efficiency
    bayer_pattern : array_like
        Layout of the detector's Bayer pattern. For example, [['B','G'],['G','R']]
        describes a BGGR pattern.
    oversample : int, optional
        Oversampling factor present in ``count``. Default is 1.
    waveunit : str, optional
        Units of ``wave``. Defaults is ``nm``

    Returns
    -------
    img : ndarray
        Electron count or flux

    Notes
    -----
    The units of ``count`` don't really matter, as long as the user is aware
    that this method converts photons per whatever to electrons per
    whatever. Whatever is nothing for counts and seconds for flux.

    """

    img = np.asarray(img)
    if img.ndim == 2:
        img = img[np.newaxis, ...]

    qe_red = qe_asarray(qe_red, wave, waveunit)
    qe_green = qe_asarray(qe_green, wave, waveunit)
    qe_blue = qe_asarray(qe_blue, wave, waveunit)

    bayer_pattern = np.char.upper(bayer_pattern)
    red_kernel = np.where(bayer_pattern == 'R', 1, 0)
    green_kernel = np.where(bayer_pattern == 'G', 1, 0)
    blue_kernel = np.where(bayer_pattern == 'B', 1, 0)

    nrow = img.shape[1] // oversample
    ncol = img.shape[2] // oversample

    # build up the bayer image. we do this one channel at a time, and then
    # accumulate the individual channels into one final frame
    red_mosaic = np.tile(red_kernel, (nrow // red_kernel.shape[0], ncol // red_kernel.shape[1]))
    red_mosaic = scipy.ndimage.zoom(red_mosaic, oversample, order=0, mode='wrap')
    red_e = np.einsum('ijk,i->jk', img, qe_red) * red_mosaic

    green_mosaic = np.tile(green_kernel, (nrow // green_kernel.shape[0], ncol // green_kernel.shape[1]))
    green_mosaic = scipy.ndimage.zoom(green_mosaic, oversample, order=0, mode='wrap')
    green_e = np.einsum('ijk,i->jk', img, qe_green) * green_mosaic

    blue_mosaic = np.tile(blue_kernel, (nrow // blue_kernel.shape[0], ncol // blue_kernel.shape[1]))
    blue_mosaic = scipy.ndimage.zoom(blue_mosaic, oversample, order=0, mode='wrap')
    blue_e = np.einsum('ijk,i->jk', img, qe_blue) * blue_mosaic

    return red_e + green_e + blue_e


def qe_asarray(qe, wave, waveunit):
    # Ensure qe is well-formed
    wave = np.asarray(wave)
    if not isinstance(qe, lentil.radiometry.Spectrum):
        qe = np.asarray(qe)
        if qe.shape == ():
            qe = qe*np.ones(wave.size)
        else:
            assert qe.size == wave.size
    else:
        qe = qe.sample(wave, waveunit=waveunit)

    return qe


def pixel(img, oversample=1):
    """Apply the aperture effects of a square pixel on a discretely sampled
    image.

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

    Example
    -------
    Apply pixel MTF to a 3x oversampled PSF:

    .. code:: pycon

        >>> import lentil
        >>> import matplotlib.pyplot as plt
        >>> psf = ...  # PSF calculation details omitted
        >>> psf_mtf = lentil.pixel(psf, oversample=3)
        >>> psf_detector = lentil.rescale(psf_mtf, 1/3)

    Note that both the pixel MTF and detector resampling operations preserve
    radiometry:

    .. code:: pycon

        >>> print(np.sum(psf), np.sum(psf_mtf), np.sum(psf_detector))
        50398.80556524441 50398.80556524441 50398.80556524441

    See Also
    --------
    :func:`lentil.detector.pixelate`

    References
    ----------
    [1] https://en.wikipedia.org/wiki/Convolution_theorem

    """

    img = np.asarray(img)
    x = np.fft.fftfreq(img.shape[1])
    y = np.fft.fftfreq(img.shape[0])

    mtf_x = np.sinc(x*oversample)
    mtf_y = np.sinc(y*oversample)
    kernel = np.dot(mtf_x[:, np.newaxis], mtf_y[np.newaxis, :])

    return np.abs(np.fft.ifft2(np.fft.fft2(img)*kernel))


def pixelate(img, oversample):
    """Convolve an image with the pixel MTF before rescaling to native sampling

    Parameters
    ----------
    img : array_like
        Input image
    oversample : int
        Number of times `img` is oversampled

    Returns
    -------
    img : ndarray
        Rescaled image with pixel MTF applied

    Note
    ----
    :func:`pixelate` should only be used if ``oversample`` > 2

    See Also
    --------
    :func:`lentil.detector.pixel`

    """

    img = lentil.detector.pixel(img, oversample)
    return lentil.rescale(img, 1/oversample, order=3, mode='nearest', unitary=True)


def adc(img, gain, saturation_capacity=None, warn_saturate=False, dtype=None):
    """Analog to digital conversion

    Parameters
    ----------
    img : ndarray
        Array of electron counts
    gain : saclar or array_like
        Conversion gain in DN/e-. Can be specified in multiple ways:

            * As a scalar term applied globally to each pixel

            * As a one-dimensional array of polynomial coefficients applied
              globally to each pixel

            * As a two-dimensional array of pixel-by-pixel scalar gain applied
              individually to each pixel

            * As a three-dimensional array of pixel-by-pixel gain where the
              first dimension gives polynomial coefficients of each pixel
    saturation_capacity : int or None
        Electron count resulting in pixel saturation. If None, pixels will not
        saturate. This is obviously nonphysical, but can be useful for testing
        or debugging.
    warn_saturate : bool, optional
        Raise a warning when pixels saturate. Default is False.
    dtype : data-type or None, optional
        Output data-type. If None (default), no data-type casting is performed.

    Returns
    -------
    img : ndarray
        Array of DN

    Note
    ----
    The saturation capacity should not be confused with the full-well capacity.
    Saturation capacity is typically smaller than the full well capacity because
    the signal is clipped before the physical saturation of the pixel is
    reached.

    """

    img = np.asarray(img)

    # Enforce saturation capacity
    if saturation_capacity:
        if warn_saturate:
            if np.any(img > saturation_capacity):
                warnings.warn('Frame has saturated pixels.')

        # Apply the saturation limit
        img[img > saturation_capacity] = saturation_capacity

    # Determine the polynomial order
    gain = np.asarray(gain)
    if gain.ndim in [0, 2]:
        if gain.ndim == 0:
            gain = gain[..., np.newaxis]
        model_order = 1
    elif gain.ndim in [1, 3]:
        model_order = gain.shape[0]
    else:
        raise ValueError

    # Prepare a cube of electron counts to apply the polynomial gain to
    img_cube = np.repeat(img[np.newaxis, :, :], model_order, axis=0)
    for order in np.arange(model_order, 1, -1):
        d = model_order - order
        img_cube[d] = img_cube[d]**order

    # Apply the gain model and convert to DN
    if gain.ndim == 1:
        img = np.einsum('ijk,i->jk', img_cube, gain)
    elif gain.ndim == 2:
        img = np.einsum('ijk,jk->jk', img_cube, gain)
    else:
        # gain.ndim == 3
        img = np.einsum('ijk,ijk->jk', img_cube, gain)

    img = np.floor(img)
    img[img < 0] = 0

    if dtype is not None:
        img = img.astype(dtype)

    return img


def shot_noise(img, method='poisson', seed=None):
    r"""Apply shot noise to an image

    Parameters
    ----------
    img : array_like
        Array of counts. All values must be >= 0.
    method : 'poisson' or 'gaussian'
        Noise method.
    seed : None, int, or array_like, optional
        Random seed used to initialize ``numpy.random.RandomState``. If
        ``None``, then ``RandomState`` will try to read data from /dev/urandom
        (or the Windows analogue) if available or seed from the clock otherwise.

    Returns
    -------
    img : ndarray
        Array of noisy counts

    Notes
    -----
    The output of the Poisson distribution is limited to the range of the C
    int64 type. A ValueError is raised when img contains values within 10 sigma
    of the maximum representable value (lam > 9.223372006484771e+18).

    For sufficiently large values of :math:`\lambda`,  (say :math:`\lambda >
    1000`), the Normal(:math:`\mu = \lambda`, :math:`\sigma^2 = \lambda`)
    distribution is an excellent approximation to the Poisson(:math:`\lambda`)
    distribution and is about 25% faster to compute.

    """
    assert method in {'poisson', 'gaussian'}

    rng = np.random.RandomState(seed)

    if method == 'poisson':
        try:
            img = rng.poisson(img)
        except ValueError as e:
            if np.min(img) < 0:
                raise ValueError('Counts must be positive')
            elif np.max(img) > 9.223372006484771e+18:
                raise ValueError('Counts exceed max representable value')
            else:
                raise e
    else:
        # REF: https://stackoverflow.com/a/33701974
        with np.errstate(divide='raise'):
            try:
                img = np.asarray(rng.normal(loc=img, scale=np.sqrt(img)), dtype=int)
            except FloatingPointError:
                raise ValueError('Counts must be positive')

    return np.floor(img)


def read_noise(img, electrons, seed=None):
    """Apply read noise to a frame

    Parameters
    ----------
    img : array_like
        Array of electrons
    electrons : int
        Read noise per frame
    seed : None, int, or array_like, optional
        Random seed used to initialize ``numpy.random.RandomState``. If
        ``None``, then ``RandomState`` will try to read data from /dev/urandom
        (or the Windows analogue) if available or seed from the clock otherwise.

    Returns
    -------
    img : ndarray
        Input image with read noise applied

    """
    img = np.asarray(img)
    rng = np.random.RandomState(seed)
    noise = rng.normal(loc=0.0, scale=electrons, size=img.shape)
    return img + noise


def charge_diffusion(img, sigma, oversample=1):
    """Apply charge diffusion represented by a Gaussian blur

    Parameters
    ----------
    img : array_like
        Frame to apply charge diffusion to
    sigma : float
        Diffusion sigma
    oversample : int
        Oversample factor in img

    Returns
    -------
    img : ndarray
        Input frame blurred by charge diffusion

    """
    kernel = lentil.util.gaussian2d(3*oversample, sigma)
    return scipy.signal.convolve2d(img, kernel, mode='same')


def dark_current(rate, shape=1, fpn_factor=0, seed=None):
    """Create dark current frame

    If applied, dark current fixed pattern noise (FPN) is modeled by a
    log-normal distribution [1], [2] with mean = 1.0 and sigma = fpn_factor
    where rate = rate * FPN.

    Parameters
    ----------
    rate : int
        Dark current rate in electrons/second/pixel
    shape : array_like or 1
        Output shape. If not specified, a scalar is returned.
    fpn_factor : float
        Dark current FPN factor. Should be between 0.1 and 0.4 for CCD and CMOS
        sensors [1].
    seed : None, int, or array_like, optional
        Random seed used to initialize ``numpy.random.RandomState``. If
        ``None``, then ``RandomState`` will try to read data from /dev/urandom
        (or the Windows analogue) if available or seed from the clock otherwise.

        .. warning::
            Be aware that if fpn_factor is nonzero and a seed is not provided,
            the fixed pattern noise will be different each time dark_current is
            called.

    Returns
    -------
    dark : ndarray
        Dark current frame

    References
    ----------
    [1] `Log-normal distribution - Wikipedia <https://en.wikipedia.org/wiki/Log-normal_distribution>`_

    [2] Janesick, J. R., Photon Transfer, Vol. PM170, SPIE (2007)

    """
    if fpn_factor > 0:
        rng = np.random.RandomState(seed)
        fpn = rng.lognormal(mean=1.0, sigma=fpn_factor, size=shape)
    else:
        fpn = 1

    dark = np.floor(rate*np.ones(shape)*fpn)
    return dark


def rule07_dark_current(temperature, cutoff_wavelength, pixelscale, shape=1,
                        fpn_factor=0, seed=None):
    """Create dark current frame for HgCdTe infrared detectors using Rule 07 [1].

    If applied, dark current fixed pattern noise (FPN) is modeled by a
    log-normal distribution [2], [3] with mean = 1.0 and sigma = fpn_factor
    where rate = rate * FPN.

    Parameters
    ----------
    temperature : float
        Focal plane temperature in K
    cutoff_wavelength : float
        Focal plane cutoff wavelength in m
    pixelscale : float
        Size of one pixel in m
    shape : array_like or 1
        Output shape. If not specified, a scalar is returned.
    fpn_factor : float
        Dark current FPN factor. Should be between 0.1 and 0.4 for CCD and CMOS
        sensors [2]. When fpn_factor is nonzero, seed must be provided. When
        fpn_factor is 0 (default), dark current FPN is not applied.
    seed : None, int, or array_like, optional
        Random seed used to initialize ``numpy.random.RandomState``. If
        ``None``, then ``RandomState`` will try to read data from /dev/urandom
        (or the Windows analogue) if available or seed from the clock otherwise.

    Returns
    -------
    dark_current  : ndarray
        Dark current

    References
    ----------
    [1] Tennant, W.E. et al. MBE HgCdTe Technology: A Very General Solution to IR
    Detection, Described by "Rule 07", a Very Convenient Heuristic. Journal of
    Electronic Materials (2008).

    [2] `Log-normal distribution - Wikipedia <https://en.wikipedia.org/wiki/Log-normal_distribution>`_

    [3] Janesick, J. R., Photon Transfer, Vol. PM170, SPIE (2007)

    """
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
    rate = 1/q * px_area * J

    return dark_current(rate, shape, fpn_factor, seed)


def cosmic_rays(shape, pixelscale, ts, rate=4e4, proton_flux=1e9, alpha_flux=4e9):
    """Cosmic ray generator for simulating cosmic ray hits on a detector.

    The default values for cosmic ray rates and electron fluxes are taken from
    [1]. It is commonly accepted at about 90% of cosmic rays are protons and the
    majority of the remaining 10% are alpha particles [1], [2].

    Parameters
    ----------
    shape : array_like
        Frame size defined as (rows, cols)
    pixelscale : array_like
        Pixel dimensions as (y,x,z) where z is the pixel depth
    ts : float
        Integration time in seconds
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

    Returns
    -------
    cosmic_e : ndarray
        Array representing the number of electrons in a detector image due
        to cosmic ray hits

    Example
    -------
    Simulate the cosmic ray hits for a 300 second exposure over a 256 x 256 detector
    patch with 5 um x 5 um x 3 um pixels:

    .. code:: pycon

        >>> import matplotlib.pyplot as plt
        >>> import lentil
        >>> cosmic_frame = lentil.detector.cosmic_rays((256,256), (5e-6, 5e-6, 3e-6), 300)
        >>> plt.imshow(cosmic_frame)

    .. image:: /_static/img/api/detector/cosmic_ray.png
        :width: 300px

    References
    ----------
    [1] Offenberg, J.D. et. al.  Multi-Read Out Data Simulator. (2000).

    [2] `Cosmic ray - Wikipedia <https://en.wikipedia.org/wiki/Cosmic_ray>`_

    """
    # compute the number of rays that strike the detector during the
    # integration time
    nrays = _nrays(shape, pixelscale, ts, rate)

    img = np.zeros(shape)
    for ray in np.arange(0, nrays):
        img += _cosmic_ray(shape, pixelscale, alpha_flux, proton_flux)
    return img


def _nrays(shape, pixelscale, ts, rate):
    area = shape[0]*pixelscale[0]*shape[1]*pixelscale[1]

    nrays = area * rate * ts

    # even if we don't have any events as defined by the rate, there is
    # still a chance we will see an event
    if nrays < 1:
        if np.random.uniform() <= nrays:
            nrays = 1
        else:
            nrays = 0

    return int(nrays)


def _cosmic_ray(shape, pixelscale, alpha_flux, proton_flux):

    # randomly choose which type of particle we have
    if np.random.uniform() > 0.9:
        # alpha particle
        electron_flux = alpha_flux
    else:
        # proton
        electron_flux = proton_flux

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
    ray = _propagate_ray(position, direction, extent)

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


def _propagate_ray(position, direction, extent):
    position = position[np.newaxis, ...]
    direction = direction[np.newaxis, ...]

    intersections = _cubeplane_ray_intersection(position, direction, extent)
    raw_ray_data = _process_cube_intersections(intersections, position, direction)

    return raw_ray_data[3]


def _cubeplane_ray_intersection(xyz, dxdydz, extent):
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
    _process_cube_intersections().

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


def _process_cube_intersections(inter, xyz, dxdydz):
    """Intersects rays with a plane of cubes and returns data about
    intersections

    Suppose the walls of the cube are on integer boundaries in the space of
    extent[0] <= x <= extent[1], extent[2] <= y <= extent[3], extent[4] <= z
    <= extent[5]

    Parameters
    ----------
    inter : list of tuple
        Output of _cubeplane_ray_intersection() or _cubeplane_ray_intersection()
    xyz : ndarray, ndim=2
        xyz[:, 0] is the x coordinate of the starting point of the ray,
        xyz[:, 1] is the y coordinate of the starting point of the ray.
        xyz[:, 2] is the z coordinate of the starting point.
    dxdydz : double ndim=2
        Unit vector representing the direction of travel of the ray.

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
    np.putmask(draylength, np.ediff1d(inter['indx'], to_begin=np.array([1], dtype=np.uint64)) != 0,
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

    return (inter['indx'], inter['raylength'], draylength, inter_point,
            last_coord, next_coord)
