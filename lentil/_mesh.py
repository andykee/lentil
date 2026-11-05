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
        Mesh rotation in degrees. When ``indexing = 'xy'`` the rotation is
        measured from the +x axis toward the +y axis. When ``indexing = 'ij'``
        the rotation is measured from the +i axis toward the +j axis. Note
        that rotations occur about the shifted center if ``shift`` is nonzero.
        Default is 0.
    indexing : {'xy', 'ij'}, optional
        Cartesian ('xy', default) or matrix ('ij') indexing of mesh.

    Returns
    -------
    list of ndarrays

    """
    shape = np.broadcast_to(np.asarray(shape), (2,))
    if indexing == 'xy':
        shape = np.flip(shape)
    angle = np.radians(angle)

    x1, x2 = np.meshgrid(np.arange(shape[0]) - np.floor(shape[0]/2) - shift[0],
                         np.arange(shape[1]) - np.floor(shape[1]/2) - shift[1],
                         indexing=indexing)

    X1 = x1 * np.cos(angle) + x2 * np.sin(angle)
    X2 = x1 * -np.sin(angle) + x2 * np.cos(angle)
    return X1, X2