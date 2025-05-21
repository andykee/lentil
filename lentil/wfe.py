import numpy as np

import lentil

__all__ = ['power_spectrum', 'translation_defocus']


def power_spectrum(mask, pixelscale, rms, half_power_freq, exp, seed=None):
    """Wavefront error defined by a Power Spectral Density (PSD) function.

    Parameters
    ----------
    mask : array_like
        Binary [0,1] mask defining pupil extent
    pixelscale : float
        Physical size of each pixel in the resulting opd
    rms : float
        RMS value of the PSD error
    half_power_freq : float
        Half-power frequency in number of cycles per unit
    exp : float
        Exponent of the inverse-power law
    seed : int, optional
        Random seed used to initialize the pseudo-random number generator. If
        seed is `None` (default), the seed will be randomly generated from
        ``/dev/urandom`` if available or the system clock.

    Returns
    -------
    wfe : ndarray
        Masked wavefront error with requested PSD

    References
    ----------
    * Sidick (2009) Power Spectral Density Specification and Analysis of
      Large Optical Surfaces

    """

    mask = np.asarray(mask)
    rng = np.random.default_rng(seed)

    # Define a frequency grid
    n, m = mask.shape
    yy, xx = np.mgrid[0:m, 0:n]
    yy = (yy - (np.floor(m / 2) + 1)) / m
    xx = (xx - (np.floor(n / 2) + 1)) / n
    dr = np.sqrt(xx * xx + yy * yy)

    # Scale the half-power point in units of cycle/px
    half_power_freq = half_power_freq * pixelscale / np.sqrt(m ** 2 + n ** 2)

    # Generate the PSD function and noise filter
    psd = 1 / (1 + (dr / half_power_freq) ** exp)
    psd[dr == 0] = 0
    psd = psd / np.sum(psd)

    H = np.fft.fftshift(np.sqrt(psd))

    n, m = mask.shape

    # Generate noise, filter it to realize the requested PSD, and enforce
    # the pupil mask
    noise = rng.normal(size=[n, m])
    opd = np.real(np.fft.ifft2(np.fft.fft2(noise) * H)) * np.sqrt(m * n)

    opd *= mask
    opd = opd * np.sqrt(np.count_nonzero(opd)/np.sum(np.abs(opd)**2)) * rms

    return opd


def translation_defocus(mask, f_number, translation):
    """Defocus error characterized by a axial translation.

    Parameters
    ----------
    mask : array_like
        Binary [0,1] mask defining the pupil extent to compute the wavefront
        over.
    f_number : float
        System or beam f-number.
    translation : float
        Axial translation from best focus given in the same units as the
        resulting OPD.

    Returns
    -------
    wfe : ndarray
        Focus error due to axial translation

    Warnings
    --------
    Defocus is computed over a circular pupal defined by the maximum
    horizontal or vertical extent of the supplied mask before applying the
    mask. As a result, the RMS and peak-to-valley wavefront measurements
    of the returned OPD will be different than what was requested for non-
    circular masks.

    """
    mask = np.asarray(mask)

    rmin, rmax, cmin, cmax = lentil.boundary(mask)
    extent = max(rmax-rmin, cmax-cmin)

    circlemask = lentil.circle(mask.shape, np.ceil(extent/2), antialias=False)

    coeff = translation/(8*f_number**2)  # p-v defocus

    z4 = lentil.zernike(circlemask, 4)
    z4 = z4/(z4.max() - z4.min())  # normalize p-v = 1

    return z4*mask*coeff


# def kolmogorov_wfe():
#     """Wavefront error representing a turbulent phase screen."""
#     pass
