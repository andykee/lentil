import numpy as np

__all__ = ['normalize_power', 'iterable_amplitude', 'iterable_mask', 'iterable_phase',
           'iterable_segmask']


def normalize_power(array, power=1):
    r"""Normalizie the power in an array.

    The total power in an array is given by

    .. math::

        P = \sum{\left|\mbox{array}\right|^2}

    A normalization coefficient is computed as

    .. math::

        c = \sqrt{\frac{p}{\sum{\left|\mbox{array}\right|^2}}}

    The array returned by a will be scaled by the normalization coefficient so
    that its power is equal to :math:`p`.

    Parameters
    ----------
    array : array_like
        Array to be normalized

    power : float, optional
        Desired power in normalized array. Default is 1.

    Returns
    -------
    array : ndarray
        Normalized array

    """
    array = np.asarray(array)
    return array * np.sqrt(power/np.sum(np.abs(array)**2))


def iterable_amplitude(iterable):
    """Construct a common amplitude array from an iterable of planes

    Parameters
    ----------
    iterable : list_like
        List of planes

    Returns
    -------
    amplitude : ndarray

    """
    amplitude = np.array(1)
    for plane in iterable:
        amplitude = amplitude * plane.amplitude
    return amplitude


def iterable_phase(iterable):
    """Construct a common phase array from an iterable of planes

    Parameters
    ----------
    iterable: list_like
        List of planes

    Returns
    -------
    phase: ndarray

    """

    phase = np.array(0)
    for plane in iterable:
        phase = phase + plane.phase
    return phase


def iterable_mask(iterable):
    """Construct a common mask array from an iterable of planes

    Parameters
    ----------
    iterable: list_like
        List of planes

    Returns
    -------
    mask: ndarray

    """
    mask = np.array(1)
    for plane in iterable:
        mask = mask * plane.mask
    return mask


def iterable_segmask(iterable):
    """Construct a common segmask array from an iterable of planes

    Parameters
    ----------
    iterable: list_like
        List of planes

    Returns
    -------
    segmask: ndarray

    """
    if all(plane.segmask is None for plane in iterable):
        return None
    else:
        segmask = np.array(1)
        for plane in iterable:
            if plane.segmask is not None:
                segmask = segmask * plane.segmask
        return segmask
