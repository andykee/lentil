# Lentil internal helper functions

import numpy as np


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
