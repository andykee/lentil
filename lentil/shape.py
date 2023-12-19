import numpy as np

import lentil

def circle(shape, radius, shift=(0, 0)):
    """Compute a circle with anti-aliasing.

    Parameters
    ----------
    shape : array_like
        Size of output in pixels (nrows, ncols)

    radius : float
        Radius of circle in pixels

    shift : (2,) array_like, optional
        How far to shift center in float (rows, cols). Default is (0, 0).

    Returns
    -------
    circle : ndarray

    """
    rr, cc = lentil.helper.mesh(shape)
    r = np.sqrt(np.square(rr - shift[0]) + np.square(cc - shift[1]))
    return np.clip(radius + 0.5 - r, 0.0, 1.0)


def circlemask(shape, radius, shift=(0, 0)):
    """Compute a circular mask.

    Parameters
    ----------
    shape : array_like
        Size of output in pixels (nrows, ncols)

    radius : float
        Radius of circle in pixels

    shift : array_like, optional
        How far to shift center in float (rows, cols). Default is (0, 0).

    Returns
    -------
    mask : ndarray

    """

    mask = lentil.circle(shape, radius, shift)
    mask[mask > 0] = 1
    return mask


def hexagon(shape, radius, rotate=False):
    """Compute a hexagon mask.

    Parameters
    ----------
    shape : array_like
        Size of output in pixels (nrows, ncols)

    radius : int
        Radius of outscribing circle (which also equals the side length) in
        pixels.

    rotate : bool
        Rotate mask so that flat sides are aligned with the Y direction instead
        of the default orientation which is aligned with the X direction.

    Returns
    -------
    mask : ndarray

    """

    inner_radius = radius * np.sqrt(3)/2
    side_length = radius/2

    rr, cc = lentil.helper.mesh(shape)

    rect = np.where((np.abs(cc) <= side_length) & (np.abs(rr) <= inner_radius))
    left_tri = np.where((cc <= -side_length) & (cc >= -radius) & (np.abs(rr) <= (cc + radius)*np.sqrt(3)))
    right_tri = np.where((cc >= side_length) & (cc <= radius) & (np.abs(rr) <= (radius - cc)*np.sqrt(3)))

    mask = np.zeros(shape)
    mask[rect] = 1
    mask[left_tri] = 1
    mask[right_tri] = 1

    if rotate:
        return mask.transpose()
    else:
        return mask


def slit(shape, width, length=None):
    """Compute a slit mask.

    Parameters
    ----------
    shape : array_like
        Size of output in pixels (nrows, ncols)
    width : float
        Slit width in pixels
    length : float, optional
        Slit length in pixels. If not specified, the slit spans
        the entire column shape (default).

    Returns
    -------
    mask : ndarray

    """

    rr, cc = lentil.helper.mesh(shape)
    slit = np.ones(shape)

    length = shape[1] if length is None else length
    width_clip = np.clip(0.5 + (width/2) - np.abs(rr), 0, 1)
    length_clip = np.clip(0.5 + (length/2) - np.abs(cc), 0, 1)

    slit = np.minimum(np.minimum(slit, width_clip), length_clip)

    return slit