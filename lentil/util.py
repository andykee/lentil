import numpy as np
from scipy.ndimage.interpolation import map_coordinates
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
    x, y = lentil.util.mesh(shape)
    r = np.sqrt(np.square(x - shift[1]) + np.square(y - shift[0]))
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

    y, x = lentil.util.mesh(shape)

    rect = np.where((np.abs(x) <= side_length) & (np.abs(y) <= inner_radius))
    left_tri = np.where((x <= -side_length) & (x >= -radius) & (np.abs(y) <= (x + radius)*np.sqrt(3)))
    right_tri = np.where((x >= side_length) & (x <= radius) & (np.abs(y) <= (radius - x)*np.sqrt(3)))

    mask = np.zeros(shape)
    mask[rect] = 1
    mask[left_tri] = 1
    mask[right_tri] = 1

    if rotate:
        return mask.transpose()
    else:
        return mask


def slit(shape, width):
    """Compute a slit mask.

    Parameters
    ----------
    shape : array_like
        Size of output in pixels (nrows, ncols)

    width : int
        Slit width in pixels

    Returns
    -------
    mask : ndarray

    """

    y, x = lentil.util.mesh(shape)

    mask = np.zeros(shape)
    mask[np.abs(x) <= width/2] = 1

    return mask


def centroid(img):
    """Compute image centroid location.

    Parameters
    ----------
    img : array_like
        Input array.

    Returns
    -------
    tuple
        ``(x,y)`` centroid location.
    """

    img = np.asarray(img)
    img = img/np.sum(img)
    nr, nc = img.shape
    yy, xx = np.mgrid[0:nr, 0:nc]

    x = np.dot(xx.ravel(), img.ravel())
    y = np.dot(yy.ravel(), img.ravel())

    return x, y


def pad(array, shape):
    """Zero-pad an array.

    Note that pad works for both two and three dimensional arrays.

    Parameters
    ----------
    array : array_like
        Array to be padded.

    shape : tuple of ints
        Shape of output array in ``(nrows, ncols)``.

    Returns
    -------
    padded : ndarray
        Zero-padded array with shape ``(nrows, ncols)``. If ``array`` has a
        third dimension, the return shape will be ``(nrows, ncols, depth)``.
    """

    array = np.asarray(array)

    offset = 0
    if array.ndim == 3:
        offset = 1


    dr = shape[0] - array.shape[0+offset]
    dc = shape[1] - array.shape[1+offset]

    if dr <= 0:
        rmin0 = (array.shape[0+offset] - shape[0])//2
        rmax0 = rmin0 + shape[0]
        rmin1 = 0
        rmax1 = shape[0]
    else:
        rmin0 = 0
        rmax0 = array.shape[0+offset]
        rmin1 = (shape[0] - array.shape[0+offset])//2
        rmax1 = rmin1 + array.shape[0+offset]

    if dc <= 0:
        cmin0 = (array.shape[1+offset] - shape[1])//2
        cmax0 = cmin0 + shape[1]
        cmin1 = 0
        cmax1 = shape[1]
    else:
        cmin0 = 0
        cmax0 = array.shape[1]
        cmin1 = (shape[1] - array.shape[1+offset])//2
        cmax1 = cmin1 + array.shape[1+offset]

    if array.ndim < 3:
        padded = np.zeros((shape[0], shape[1]), dtype=array.dtype)
        padded[rmin1:rmax1, cmin1:cmax1] = array[rmin0:rmax0, cmin0:cmax0]
    else:
        padded = np.zeros((array.shape[0], shape[0], shape[1]), dtype=array.dtype)
        padded[:, rmin1:rmax1, cmin1:cmax1] = array[:, rmin0:rmax0, cmin0:cmax0]

    return padded


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
      is trimmed about the center of the array using :func:`lentil.util.pad`.
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

    Note
    ----
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


def slice_offset(slice, shape, indexing='xy'):
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

    indexing : {'xy', 'ij'}, optional
        Offset ordering. Default is 'xy'.

    Returns
    -------
    offset : tuple or None
        Offset ordered according to ``indexing``. Note that if the computed offset
        is (0, 0), None is returned instead.

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
        offset = None
    elif Ellipsis in slice:
        # The only case we know enough to deal with is (Ellipsis, slice(None, None, None))
        if slice(None, None, None) in slice:
            offset = None
        else:
            raise ValueError(f"Can't compute offset from slice {slice}")
    else:
        slice_shape = np.array((slice[0].stop-slice[0].start-1, slice[1].stop-slice[1].start-1))
        slice_center = slice_shape - np.floor(slice_shape/2)

        shape = np.asarray(shape)
        center = shape - np.floor(shape/2)

        slice_offset = np.array((slice[0].start+slice_center[0], slice[1].start+slice_center[1])) - center

        if np.all(slice_offset == 0):
            offset = None
        else:
            if indexing == 'xy':
                offset = tuple(slice_offset[::-1])
            elif indexing == 'ij':
                offset = tuple(slice_offset)
            else:
                raise ValueError(f"Unknown indexing {indexing}. indexing must be 'ij' or 'xy'.")

    return offset


def mesh(shape, shift=(0, 0)):
    """Generate a standard mesh."""

    nr = shape[0]
    nc = shape[1]

    x = np.arange(shape[1]) - np.floor(shape[1]/2.0) - shift[1]
    y = np.arange(shape[0]) - np.floor(shape[0]/2.0) - shift[0]

    return np.meshgrid(x, y)


def gaussian2d(size, sigma):
    """2D Gaussian kernel."""
    x, y = np.meshgrid(np.linspace(-1, 1, size), np.linspace(-1, 1, size))

    G = np.exp(-((x**2/(2*sigma**2)) + (y**2/(2*sigma**2))))
    return G/np.sum(G)


def make_index(mat):
    """Create a sparse coordinate list (COO) index dictionary.

    Parameters
    ----------
    mat : array_like
        Full matrix to ve vectorized

    Returns
    -------
    dict
        Index dictionary with the following attributes:

        * ``row`` List of row indices which contain nonzero data
        * ``col`` List of column indices which contain nonzero data
        * ``shape`` Tuple of dimensions of :attr:`~lentil.wfe.make_index.mat`

    See Also
    --------
    * :func:`~lentil.sparse.v2m` Convert sparse vectorized data to a full matrix
    * :func:`~lentil.sparse.m2v` Convert a full matrix to sparse vectorized data

    """
    mat = np.asarray(mat)
    row, col = mat.nonzero()
    shape = mat.shape
    return {'row': row, 'col': col, 'shape': shape}


def make_mask(index):
    """Create a mask from a sparse oordinate list (COO) index dictionary."""
    return v2m(np.ones(index['row'].shape), index)


def v2m(vec, index):
    """Convert a sparse vector to a full matrix.

    Parameters
    ----------
    vec : array_like
        Sparse vector to be reformed as a full matrix

    index : dict
        Corresponding index dictionary

    Returns
    -------
    ndarray
        Full matrix

    See Also
    --------
    * :func:`~lentil.sparse.m2v` Convert a full matrix to sparse vectorized data
    * :func:`~lentil.sparse.make_index` Create a sparse coordinate list (COO)
      index dictionary

    """
    mat = np.zeros(index['shape'])
    nnz = len(index['row'])
    for i in range(nnz):
        mat[index['row'][i], index['col'][i]] = vec[i]
    return mat


def m2v(mat, index):
    """Convert a full matrix to a sparse vector.

    Parameters
    ----------
    mat : array_like
        Full matrix to be reformed as a sparse vector

    index : dict
        Corresponding index dictionary

    Returns
    -------
    ndarray
        Sparse vector

    See Also
    --------
    * :func:`~lentil.sparse.v2m` Convert sparse vectorized data to a full matrix
    * :func:`~lentil.sparse.make_index` Create a sparse coordinate list (COO)
      index dictionary

    """
    raveled_indices = np.ravel_multi_index((index['row'], index['col']),
                                           index['shape'], order='C')
    vec = mat.ravel()
    return vec[raveled_indices]


def omega(f_number):
    return np.pi*np.sin(np.arctan(1/(2*f_number)))**2


def expc(x):
    """Computes np.exp(1.0j * x) for real valued x... FAST!"""
    x = np.asarray(x)
    if x.dtype == np.float32:
        out = np.empty(x.shape, dtype=np.complex64)
    elif x.dtype == np.float64:
        out = np.empty(x.shape, dtype=np.complex128)

    out.real = np.cos(x)
    out.imag = np.sin(x)
    return out
