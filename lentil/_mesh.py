import numpy as np


def mesh(shape, shift=(0, 0), angle=0, indexing='xy'):
    """Generate a standard mesh.

    Parameters
    ----------
    shape : array_like
        Size of output in pixels (nrows, ncols)
    shift : (2,) array_like, optional
        Coordinate translation from center of containing array.
        Default is (0, 0). Shift is applied according to the specified
        ``indexing`` scheme.
    angle : float, optional
        Mesh rotation in degrees. Default is 0. When ``indexing = 'xy'``,
        the rotation is measured from the +x axis toward the +y axis, which
        appears counterclockwise when displayed with ``origin='lower'``.
        When ``indexing = 'ij'``, the rotation is measured from the +i axis
        toward the +j axis. Note that rotations occur about the shifted center if 
        ``shift`` is nonzero.
    indexing : {'xy', 'ij'}, optional
        Cartesian ('xy', default) or matrix ('ij') indexing of mesh.

    Returns
    -------
    list of ndarrays

    """
    shape = np.broadcast_to(np.asarray(shape), (2,))
    angle = np.radians(angle)
    if indexing == 'xy':
        x1, x2 = _meshxy(shape, shift, angle)
    elif indexing == 'ij':
        x1, x2 = _meshij(shape, shift, angle)
    else:
        raise ValueError(
            "Valid values for indexing are 'xy' and 'ij'.")
    return x1, x2


def _meshxy(shape, shift, angle):
    # return x and y coordinate matrices where x increases as columns increase
    # and y increases as rows increase. angle specifies a CCW rotation from
    # the +x axis
    nx = shape[1]
    ny = shape[0]
    xx, yy = np.meshgrid(np.arange(nx) - np.floor(nx/2.0) - shift[0],
                         np.arange(ny) - np.floor(ny/2.0) - shift[1],
                         indexing='xy')
    x = xx * np.cos(angle) + yy * np.sin(angle)
    y = xx * -np.sin(angle) + yy * np.cos(angle)
    return x, y


def _meshij(shape, shift, angle):
    nr = shape[0]
    nc = shape[1]
    rr, cc = np.meshgrid(np.arange(nr) - np.floor(nr/2.0) - shift[0],
                         np.arange(nc) - np.floor(nc/2.0) - shift[1],
                         indexing='ij')
    r = rr * np.cos(angle) + cc * np.sin(angle)
    c = rr * -np.sin(angle) + cc * np.cos(angle)
    return r, c
