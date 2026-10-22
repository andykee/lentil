import warnings

import numpy as np
from scipy.ndimage import map_coordinates
import lentil


def centroid(x, where=None, mode='absolute', indexing='xy'):
    """Compute array centroid location.

    Parameters
    ----------
    x : array_like
        Input array.
    where: array_like of bool, optional
        Elements to include in the centroid calculation. If None (default),
        all finite and non-NaN values are used.
    mode : {'absolute', 'center'}, optional
        Specifies how the centroid is computed. If 'absolute' (default),
        the absolute centroid within the input is returned. If 'center', the
        centroid relative to the center of the input is returned.
    indexing : {'xy', 'ij'}, optional
        Cartesian ('xy', default) or matrix ('ij', default) indexing of mesh.

    Returns
    -------
    centroid : tuple
        (r, c) or (x, y) centroid location

    Examples
    --------
    .. plot::
        :include-source:
        :context: reset
        :scale: 50

        >>> circ = lentil.circle(shape=(256, 256), radius=25, shift=(-80, 50))
        >>> plt.imshow(circ, cmap='gray')

    .. code:: pycon

        >>> lentil.centroid(circ)
        (48.0, 178.0)
        >>> lentil.centroid(circ, mode='center')
        (-80.0, 50.0)
        >>> lentil.centroid(circ, mode='center', indexing='xy')
        (50.0, -80.0)

    """
    x = np.asarray(x)

    if mode not in ('absolute', 'center'):
        raise ValueError(f'Unknown kind {mode}')

    if indexing not in ('ij', 'xy'):
        raise ValueError("Valid values for indexing are 'xy' and 'ij'.")

    if where is None:
        where = np.isfinite(x)
    else:
        where = np.asarray(where, dtype=bool)

    if np.isnan(x[where]).any():
        warnings.warn('Unmasked NaN in input', RuntimeWarning,
                      stacklevel=2)

    anorm = x[where]/np.sum(x[where])

    nr, nc = x.shape
    rr, cc = np.mgrid[0:nr, 0:nc]

    r = np.dot(rr[where].ravel(), anorm.ravel())
    c = np.dot(cc[where].ravel(), anorm.ravel())

    if mode == 'center':
        rc, cc = np.floor(np.asarray(x.shape)/2)
        r = r - rc
        c = c - cc

    if indexing == 'xy':
        if mode == 'center':
            # Cartesian offset from center: +x right, +y up
            r, c = c, -r
        else:
            # absolute pixel coordinates (col, row); +y down
            r, c = c, r

    return r, c


def pad(x, shape, fill=0):
    """Pad or crop an array about its center.

    Note that pad accepts both two and three dimensional arrays.

    Parameters
    ----------
    x : array_like
        Array to be padded.
    shape : array_like of ints
        Shape of output array in ``(nrows, ncols)``.
    fill : scalar
        Fill value used when pad operation increases array size.

    Returns
    -------
    padded_array : ndarray
        Zero-padded array with shape ``(nrows, ncols)``. If ``x`` has a
        third dimension, the return shape will be ``(depth, nrows, ncols)``.

    See Also
    --------
    :func:`crop`

    Examples
    --------
    .. plot::
        :include-source:
        :context: reset
        :scale: 50

        >>> circ = lentil.circle(shape=(128, 128), radius=64)
        >>> circ_pad = lentil.pad(circ, shape=(200, 200))
        >>> fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(5, 2))
        >>> ax[0].imshow(circ, cmap='gray')
        >>> ax[0].set_title('Original array')
        >>> ax[1].imshow(circ_pad, cmap='gray')
        >>> ax[1].set_title('Padded array')

    .. plot::
        :include-source:
        :context: reset
        :scale: 50

        >>> circ = lentil.circle(shape=(128, 128), radius=64)
        >>> circ_pad = lentil.pad(circ, shape=(110, 110))
        >>> fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(5, 2))
        >>> ax[0].imshow(circ, cmap='gray')
        >>> ax[0].set_title('Original array')
        >>> ax[1].imshow(circ_pad, cmap='gray')
        >>> ax[1].set_title('Padded array')

    """
    x = np.atleast_2d(x)
    shape = np.broadcast_to(np.asarray(shape, dtype=int), (2,))

    squeeze = x.ndim == 2
    if squeeze:
        out_shape = (1, int(shape[0]), int(shape[1]))
        x = x[np.newaxis, :]
    else:  # a.ndim == 3
        out_shape = (x.shape[0], int(shape[0]), int(shape[1]))

    # The row and col indices here are 1 and 2 respectively since we've
    # forced both x and out to have a 3rd dimension
    rmin = min(x.shape[1] // 2, out_shape[1] // 2)
    rmax = min(x.shape[1] - x.shape[1] // 2, out_shape[1] - out_shape[1] // 2)
    cmin = min(x.shape[2] // 2, out_shape[2] // 2)
    cmax = min(x.shape[2] - x.shape[2] // 2, out_shape[2] - out_shape[2] // 2)

    # Extract the overlapping region from the source array
    a_slc_r = slice(x.shape[1] // 2 - rmin, x.shape[1] // 2 + rmax)
    a_slc_c = slice(x.shape[2] // 2 - cmin, x.shape[2] // 2 + cmax)
    cropped = x[:, a_slc_r, a_slc_c]

    # Pad the cropped region to the output shape
    pad_r_before = out_shape[1] // 2 - rmin
    pad_r_after = out_shape[1] - out_shape[1] // 2 - rmax
    pad_c_before = out_shape[2] // 2 - cmin
    pad_c_after = out_shape[2] - out_shape[2] // 2 - cmax

    out = np.pad(cropped,
                 ((0, 0), (pad_r_before, pad_r_after), (pad_c_before, pad_c_after)),
                 mode='constant', constant_values=fill)

    # only remove the depth axis added above
    return out[0] if squeeze else out


def subarray(a, shape, shift=(0,0)):
    """Extract a contiguous subarray from a larger array.

    The subarray is extracted about the center of the source array unless
    a shift is specified.

    Parameters
    ----------
    a : array_like
        Source array
    shape : array_like of ints
        Shape of subarray array in ``(nrows, ncols)``.
    shift : array_like of ints
        Relative shift of the center of the subarray in ``(row, col)``.
    
    Returns
    -------
    out : ndarray
        Subarray extracted from the source array.
    """

    a = np.asarray(a)
    shape = np.asarray(shape)

    rmin = a.shape[0]//2 - shape[0]//2 + shift[0]
    cmin = a.shape[1]//2 - shape[1]//2 + shift[1]
    rmax = rmin + shape[0]
    cmax = cmin + shape[1]

    if any((rmin<0, cmin<0, rmax>a.shape[0], cmax>a.shape[1])):
        raise ValueError('window lies outside of array')

    return a[rmin:rmax, cmin:cmax]


def window(img, shape=None, slice=None):
    """Extract an appropriately sized, potentially windowed array

    Parameters
    ----------
    img : array_like
        Data to window

    shape : array_like or None, optional
        Output shape given as (nrows, ncols). If ``None`` (default), an
        uncropped array is returned.

    slice : array_like or None, optional
        Indices of ``img`` array to return given as (r_start, r_end, c_start,
        c_end). This definition follows standard numpy indexing.

    Returns
    -------
    data : ndarray
        Trimmed and windowed ``img`` array

    Notes
    -----
    * If ``img`` is a single value (img.size == 1), self.data is returned
      regardless of what ``shape`` and ``slice`` are.
    * If ``shape`` is given but ``slice`` is ``None``, the returned ndarray
      is trimmed about the center of the array using :func:`lentil.pad`.
    * If ``slice`` is given but ``shape`` is ``None``, the returned ndarray
      is extracted from ``img`` according to the indices in ``slice``
    * If both ``shape`` and ``slice`` are given, the returned ndarray is
      extracted from ``img`` according to the indices in ``slice`` and the
      following expressions must also be true:

    .. code::

        shape[0] = (slice[1] - slice[0]) = (r_end - r_start)
        shape[1] = (slice[3] - slice[2]) = (c_end - c_start)

    """

    img = np.asarray(img)
    if img.size == 1:
        return img

    if shape is None and slice is None:
        return img
    elif slice is not None:
        if shape is not None:
            # ensure size consistency
            assert(slice[1] - slice[0]) == shape[0]
            assert(slice[3] - slice[2]) == shape[1]

        # return the requested view. Note that numpy will implicitly handle a
        # third dimension if one is present
        return img[slice[0]:slice[1], slice[2]:slice[3]]

    else:
        # return the padded array. Note that pad will handle a third dimension
        # if one exists
        return lentil.pad(img, shape)


def boundary(x, threshold=0):
    """Find bounding row and column indices of data within an array.

    Parameters
    ----------
    x : array_like
        Input array

    threshold : float, optional
        Masking threshold to apply before boundary finding. Only values
        in x that are larger than threshold are considered in the boundary
        finding operation. Default is 0.

    Returns
    -------
    rmin, rmax, cmin, cmax : ints
        Boundary indices

    """
    x = np.asarray(x)
    x = (x > threshold)

    rows = np.any(x, axis=1)
    cols = np.any(x, axis=0)

    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]

    return rmin, rmax, cmin, cmax


def sample(x, oversample, offset=0, axes=(-2, -1), indexing='ij'):
    """Extract samples from an oversampled array.

    The sampled axes are tiled into contiguous blocks of ``oversample``
    samples, with one sample being taken from each block. No averaging or
    filtering is performed.

    Parameters
    ----------
    x : array_like
        Array to sample
    oversample : int or tuple of ints
        Decimation factor. Providing two values defines asymmetric
        oversampling with ordering interpreted according to ``indexing``.
    offset : int or tuple of ints, optional
        Position of the sample with each block, measured from
        ``oversample//2``. The sub-sampled position reachable by ``offset`` is
        quantized to ``1/oversample``. Default is 0. If two values are
        provided, their ordering is interpreted according to ``indexing``.
    axes : (2,) sequence of ints, optional
        Axes over which to sample. If not provided, the last two axes are
        used.
    indexing : {'ij', 'xy'}, optional
        Matrix (‘ij’, default) or cartesian (‘xy’) indexing of ``x`` grid.

    Returns
    -------
    x : ndarray
        Each sampled axis has length ``x.shape[axis]//oversample``. If axis
        length is not a multiple of ``oversample``, samples taken beyond the
        last whole block are ignored.

    """
    x = np.asarray(x)

    if indexing not in ('ij', 'xy'):
        raise ValueError("Valid values for indexing are 'xy' and 'ij'.")

    oversample = tuple(int(v) for v in np.broadcast_to(oversample, (2,)))
    offset = tuple(int(v) for v in np.broadcast_to(offset, (2,)))
    axes = tuple(int(a) % x.ndim for a in axes)

    if indexing == 'xy':
        offset = (-offset[1], offset[0])
        oversample = (oversample[1], oversample[0])

    idx = [slice(None)] * x.ndim

    for ax, m, o in zip(axes, oversample, offset):

        lo, hi = -(m // 2), m - 1 - m // 2
        if not lo <= o <= hi:
            raise ValueError(f'offset {o} lies outside one pixel for oversample '
                             f'{m} on axis {ax}; require {lo} <= offset <= {hi}')

        n = x.shape[ax] // m
        start = m // 2 + o
        stop = start + m * n

        idx[ax] = slice(start, stop, m)

    return x[tuple(idx)]


def rebin(img, factor):
    """Rebin an image by an integer factor.

    Parameters
    ----------
    img : array_like
        Array or cube of arrays to rebin. If a cube is provided, the first dimension
        should index the image slices.

    factor : int
        Rebinning factor

    Returns
    -------
    img : ndarray
        Rebinned image

    See Also
    --------
    :func:`rescale`

    """
    img = np.asarray(img)

    if np.iscomplexobj(img):
        raise ValueError('rebin is not defined for complex data')

    if img.ndim == 3:
        rebinned_shape = (img.shape[0], img.shape[1]//factor, img.shape[2]//factor)
        img_rebinned = np.zeros(rebinned_shape, dtype=img.dtype)
        for i in range(img.shape[0]):
            img_rebinned[i] = img[i].reshape(rebinned_shape[1], factor,
                                             rebinned_shape[2], factor).sum(-1).sum(1)
    else:
        img_rebinned = img.reshape(img.shape[0]//factor, factor, img.shape[1]//factor,
                                   factor).sum(-1).sum(1)

    return img_rebinned


def rescale(img, scale, shape=None, mask=None, order=3, mode='nearest',
            unitary=True):
    """Rescale an image by interpolation.

    Parameters
    ----------
    img : array_like
        Image to rescale

    scale : float
        Scaling factor. Scale factors less than 1 will shrink the image. Scale
        factors greater than 1 will grow the image.

    shape : array_like or int, optional
        Output shape. If None (default), the output shape will be the input img
        shape multiplied by the scale factor.

    mask : array_like, optional
        Binary mask applied after rescaling. If None (default), a mask is
        created from the nonzero portions of img. To skip masking operation,
        set ``mask = np.ones_like(img)``

    order : int, optional
        Order of spline interpolation used for rescaling operation. Default is
        3. Order must be in the range 0-5.

    mode : {'constant', 'nearest', 'reflect', 'wrap'}, optional
        Points outside the boundaries of the input are filled according to the
        given mode. Default is 'constant'.

    unitary : bool, optional
        Normalization flag. If True (default), a normalization is performed on
        the output such that the rescaling operation is unitary and image power
        (if complex) or intensity (if real) is conserved.

    Returns
    -------
    img : ndarray
        Rescaled image.

    Notes
    -----
    The post-rescale masking operation should have no real effect on the
    resulting image but is included to eliminate interpolation artifacts that
    sometimes appear in large clusters of zeros in rescaled images.

    See Also
    --------
    :func:`rebin`

    """

    img = np.asarray(img)

    if mask is None:
        # take the real portion to ensure that even if img is complex, mask will
        # be real
        mask = np.zeros_like(img).real
        mask[img != 0] = 1

    if shape is None:
        shape = np.ceil((img.shape[0]*scale, img.shape[1]*scale)).astype(int)
    else:
        if np.isscalar(shape):
            shape = np.ceil((shape*scale, shape*scale)).astype(int)
        else:
            shape = np.ceil((shape[0]*scale, shape[1]*scale)).astype(int)

    x = (np.arange(shape[1], dtype=np.float64) - shape[1]/2.)/scale + img.shape[1]/2.
    y = (np.arange(shape[0], dtype=np.float64) - shape[0]/2.)/scale + img.shape[0]/2.

    xx, yy = np.meshgrid(x, y)

    mask = map_coordinates(mask, [yy, xx], order=1, mode='nearest')
    mask[mask < np.finfo(mask.dtype).eps] = 0

    if np.iscomplexobj(img):
        out = np.zeros(shape, dtype=np.complex128)
        out.real = map_coordinates(img.real, [yy, xx], order=order, mode=mode)
        out.imag = map_coordinates(img.imag, [yy, xx], order=order, mode=mode)
    else:
        out = map_coordinates(img, [yy, xx], order=order, mode=mode)

    if unitary:
        out *= np.sum(img)/np.sum(out)

    out *= mask

    return out


def pixelscale_nyquist(wave, f_number):
    """Compute the output plane sampling which is Nyquist sampled for
    intensity.

    Parameters
    ----------
    wave : float
        Wavelength in meters

    f_number : float
        Optical system F/#

    Returns
    -------
    float
        Sampling in meters

    """
    return f_number * wave / 2


def min_sampling(wave, z, du, shape, min_q):
    num = np.min(wave) * z
    return num/(min_q * du[0] * shape[0]), num/(min_q * du[1] * shape[1])


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


def sanitize_shape(shape):
    shape = np.asarray(shape)
    if shape.size == 0:
        shape = ()
    else:
        if shape.shape == ():
            shape = np.append(shape, shape)
    return tuple(shape)


def sanitize_bandpass(vec):
    vec = np.asarray(vec)
    if vec.shape == ():
        vec = vec[np.newaxis, ...]
    return vec
