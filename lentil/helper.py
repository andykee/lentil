# Lentil internal helper functions

import numpy as np

import lentil

def mesh(shape, shift=(0, 0)):
    """Generate a standard mesh."""

    nr = shape[0]
    nc = shape[1]

    r = np.arange(nr) - np.floor(nr/2.0) - shift[0]
    c = np.arange(nc) - np.floor(nc/2.0) - shift[1]

    return np.meshgrid(r, c, indexing='ij')


def gaussian2d(size, sigma):
    """2D Gaussian kernel."""
    x, y = np.meshgrid(np.linspace(-1, 1, size), np.linspace(-1, 1, size))

    G = np.exp(-((x**2/(2*sigma**2)) + (y**2/(2*sigma**2))))
    return G/np.sum(G)


def sanitize_shape(shape, default=()):
    if shape is None:
        shape = default
    shape = np.asarray(shape)
    if shape.shape == ():
        shape = np.append(shape, shape)
    return shape


def sanitize_bandpass(vec, default=()):
    if vec is None:
        vec = default
    vec = np.asarray(vec)
    if vec.shape == ():
        vec = vec[np.newaxis, ...]
    return vec


def dft_alpha(dx, du, wave, z, oversample):
    return (dx*du)/(wave*z*oversample)


def boundary_slice(x, threshold=0, pad=(0, 0)):
    """Find bounding row and column indices of data within an array and
    return the results as slice objects.

    Parameters
    ----------
    x : array_like
        Input array

    threshold : float, optional
        Masking threshold to apply before boundary finding. Only values
        in x that are larger than threshold are considered in the boundary
        finding operation. Default is 0.

    pad : int or tuple of ints
        Additional number of pixels to pad the boundary finding result by.
        Default is (0,0).

    Returns
    -------
    row_slice, col_slice : tuple of slices
        Boundary slices

    """
    pad = np.asarray(pad)
    if pad.shape == ():
        pad = np.append(pad, pad)

    rmin, rmax, cmin, cmax = lentil.boundary(x, threshold)

    rmin = np.max((rmin-pad[0], 0))
    rmax = np.min((rmax+pad[0]+1, x.shape[0]))
    cmin = np.max((cmin-pad[1], 0))
    cmax = np.min((cmax+pad[1]+1, x.shape[1]))

    return np.s_[rmin:rmax, cmin:cmax]


def slice_offset(slice, shape):
    """Compute the offset of the center of a 2D slice relative to the center of a larger
    array.

    It is assumed the center of the larger containing array with shape (m,n) is at:

        r = m - np.floor(m/2)
        c = n - np.floor(n/2)

    Parameters
    ----------
    slice : (2,) array_like
        2D slice. Entries must be ``slice`` objects or ``Ellipsis``.
    shape : (2,) array_like
        Shape of containing array the slice is taken from.

    Returns
    -------
    offset : tuple or None
        Offset ordered according to ``indexing``.

    See Also
    --------
    :func:`boundary_slice`

    """

    # There are a couple of different ways numpy arrays can be sliced with Ellipsis:
    # np.s_[...]      -> Ellipsis
    # np.s_[..., :]   -> (Ellipsis, slice(None, None, None))
    # np.s_[..., 2]   -> (Ellipsis, 2)
    # np.s_[..., 2:4] -> (Ellipsis, slice(2, 4, None))
    # np.s_[..., ...] -> (Ellipsis, Ellipsis) ! Note that trying to slice a numpy array
    #                    like this will return IndexError: an index can only have a single
    #                    ellipsis ('...')

    if slice == Ellipsis:
        offset = (0, 0)
    elif Ellipsis in slice:
        # The only case we know enough to deal with is (Ellipsis, slice(None, None, None))
        if slice(None, None, None) in slice:
            offset = (0, 0)
        else:
            raise ValueError(f"Can't compute offset from slice {slice}")
    else:
        slice_shape = np.array((slice[0].stop-slice[0].start-1, slice[1].stop-slice[1].start-1))
        slice_center = slice_shape//2

        shape = np.asarray(shape)
        center = shape//2

        slice_offset = np.array((slice[0].start+slice_center[0], slice[1].start+slice_center[1])) - center

        if np.all(slice_offset == 0):
            offset = (0, 0)
        else:
            offset = tuple(slice_offset)

    return offset
