import numpy as np
from scipy.ndimage.interpolation import map_coordinates


def circle(shape, radius, center=(0, 0)):
    """Compute a circle with anti-aliasing.

    Parameters
    ----------
    shape : array_like
        Size of output in pixels (nrows, ncols)

    radius : float
        Radius of circle in pixels

    center : (2,) array_like, optional
        How far to shift center in float (rows, cols). Default is (0, 0).

    Returns
    -------
    circle : ndarray

    """
    x, y = mesh(shape)
    r = np.sqrt(np.square(x - center[1]) + np.square(y - center[0]))
    return np.clip(radius + 0.5 - r, 0.0, 1.0)


def circlemask(shape, radius, center=(0, 0)):
    """Compute a circular mask.

    Parameters
    ----------
    shape : array_like
        Size of output in pixels (nrows, ncols)

    radius : float
        Radius of circle in pixels

    center : array_like, optional
        How far to shift center in float (rows, cols). Default is (0, 0).

    Returns
    -------
    mask : ndarray

    """

    mask = circle(shape, radius, center)
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

    y, x = mesh(shape)

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

    y, x = mesh(shape)

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

    output_shape = (shape[0], shape[1])
    if array.ndim == 3:
        output_shape = np.append(output_shape, array.shape[2])

    padded = np.zeros(output_shape, dtype=array.dtype)

    dr = shape[0] - array.shape[0]
    dc = shape[1] - array.shape[1]

    if dr <= 0:
        rmin0 = (array.shape[0] - shape[0])//2
        rmax0 = rmin0 + shape[0]
        rmin1 = 0
        rmax1 = shape[0]
    else:
        rmin0 = 0
        rmax0 = array.shape[0]
        rmin1 = (shape[0] - array.shape[0])//2
        rmax1 = rmin1 + array.shape[0]

    if dc <= 0:
        cmin0 = (array.shape[1] - shape[1])//2
        cmax0 = cmin0 + shape[1]
        cmin1 = 0
        cmax1 = shape[1]
    else:
        cmin0 = 0
        cmax0 = array.shape[1]
        cmin1 = (shape[1] - array.shape[1])//2
        cmax1 = cmin1 + array.shape[1]

    padded[rmin1:rmax1, cmin1:cmax1] = array[rmin0:rmax0, cmin0:cmax0]

    return padded


def boundary(img, threshold=0):
    """Find bounding row and column indices of data within an array."""
    img = np.asarray(img)
    img = (img > threshold)

    rows = np.any(img, axis=1)
    cols = np.any(img, axis=0)

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


def gaussian2d(size, sigma):
    """2D Gaussian kernel."""
    x, y = np.meshgrid(np.linspace(-1, 1, size), np.linspace(-1, 1, size))

    G = np.exp(-((x**2/(2*sigma**2)) + (y**2/(2*sigma**2))))
    return G/np.sum(G)


def mesh(shape, shift=(0, 0)):
    """Generate a standard mesh."""

    nr = shape[0]
    nc = shape[1]

    x = np.arange(shape[1]) - np.floor(shape[1]/2.0) - shift[1]
    y = np.arange(shape[0]) - np.floor(shape[0]/2.0) - shift[0]

    return np.meshgrid(x, y)


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
        mask = np.zeros_like(img)
        mask[img != 0] = 1

    if shape is None:
        shape = (int(img.shape[0]*scale), int(img.shape[1]*scale))
    else:
        if np.isscalar(shape):
            shape = (int(shape*scale), int(shape*scale))
        else:
            shape = (int(shape[0]*scale), int(shape[1]*scale))

    x = (np.arange(shape[1], dtype=np.float64) - shape[1]/2.)/scale + img.shape[1]/2.
    y = (np.arange(shape[0], dtype=np.float64) - shape[0]/2.)/scale + img.shape[0]/2.

    xx, yy = np.meshgrid(x, y)

    mask = map_coordinates(mask, [yy, xx], order=1, mode='nearest')
    mask[mask < np.finfo(mask.dtype).eps] = 0

    if np.iscomplexobj(img):
        out = np.zeros(shape, dtype=np.complex128)
        out.real = map_coordinates(img.real, [yy, xx], order=order, mode=mode)
        out.imag = map_coordinates(img.imag, [yy, xx], order=order, mode=mode)
        if unitary:
            out /= scale
    else:
        out = map_coordinates(img, [yy, xx], order=order, mode=mode)
        if unitary:
            out /= scale**2

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


def col_major_to_row_major(col_major_sens, col_major_index):
    """Convert a column-major sensitivity matrix and corresponding index
    dictionary (MATLAB-style) to a row-major sensitivity matrix and
    corresponding index disctionary (Python-style).

    Parameters
    ----------
    col_major_sens : array_like
        Sensitivity matrix ordered as column-major

    col_major_index : dict
        Corresponding index dictionary

    Returns
    -------
    ndarray
        Sensitivity matrix reordered as row-major
    dict
        Corresponding index dictionary

    Warning
    -------
    Keep in mind that when importing MATLAB-created index files to create
    :attr:`~lentil.wfe.col_major_to_row_major.col_major_index`, MATLAB used
    1-based indexing and Python uses 0-based indexing. As a result, both the
    row and col lists in
    :attr:`~lentil.wfe.col_major_to_row_major.col_major_index` must have one
    subtracted from them to ensure consistency. See example below.

    Examples
    --------
    .. code:: python3

        import numpy as np
        import lentil

        # import the column-major sensitivity matrix
        dwdx = np.genfromtxt('dwdx.csv', delimiter=',')

        # index_csv contains the row and column indices as its first two
        # columns. we'll have to manually make note of shape and directly
        # specify it
        index_csv = np.genfromtxt('index.csv', delimiter=',')

        index = {}
        index['row'] = np.asarray(index_csv[:,0], dtype=int)
        index['col'] = np.asarray(index_csv[:,1], dtype=int)
        index['shape'] = (255,255)

        # convert from 1-based (MATLAB) indexing to 0-based (Python) indexing
        index['row'] -= 1
        index['col'] -= 1

        # do the conversion
        dwdx, index = lentil.util.col_major_to_row_major(dwdx, index)

    See Also
    --------
    * :func:`~lentil.sparse.v2m` Convert sparse vectorized data to a full matrix
    * :func:`~lentil.sparse.m2v` Convert a full matrix to sparse vectorized data

    """

    col_major_sens = np.asarray(col_major_sens)

    # handle 1-dim vectors
    if col_major_sens.ndim == 1:
        col_major_sens = col_major_sens[..., np.newaxis]

    # make row_major_index
    mat = np.zeros(col_major_index['shape'])
    nnz = len(col_major_index['row'])
    for i in range(nnz):
        mat[col_major_index['row'][i], col_major_index['col'][i]] = 1
    row_major_index = make_index(mat)

    row_major_sens = np.zeros(col_major_sens.shape)

    for i in range(col_major_sens.shape[1]):
        m = v2m(col_major_sens[:, i], col_major_index)
        row_major_sens[:, i] = m2v(m, row_major_index)

    return row_major_sens, row_major_index


def omega(f_number):
    return np.pi*np.sin(np.arctan(1/(2*f_number)))**2
