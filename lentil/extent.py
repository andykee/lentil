# Functions for working with extent (rmin, rmax, cmin, cmax) data

import numpy as np

def array_extent(shape, shift, parent_shape=None):
    """Compute the extent of a shifted array.

    Parameters
    ----------
    shape : (2,) array_like
        Array shape
    shift : (2,) array_like
        Array shift in (r, c)
    parent_shape : (2,) array like or None, optional
        Enclosing parent shape. If None, the returned extent is relative
        to the origin (0,0). If provided, the returned extent is relative
        to the upper left corner of the parent shape.

    Notes
    -----
    To use the values returned by ``extent()`` in a slice, 
    ``rmax`` and ``cmax`` should be increased by 1.
    """

    if len(shape) < 2:
        shape = (1, 1)

    rmin = int(-(shape[0]//2) + shift[0])
    cmin = int(-(shape[1]//2) + shift[1])
    rmax = int(rmin + shape[0] - 1)
    cmax = int(cmin + shape[1] - 1)

    if parent_shape is not None:
        parent_center = np.asarray(parent_shape)//2
        rmin += parent_center[0]
        rmax += parent_center[0]
        cmin += parent_center[1]
        cmax += parent_center[1]

    return rmin, rmax, cmin, cmax


def array_center(extent):
    """Compute the center of an extent

    Parameters
    ----------
    a : (4,) array_like
        Array extent (rmin, rmax, cmin, cmax)

    Returns
    -------
    tuple
    """

    rmin, rmax, cmin, cmax = extent
    nrow = rmax - rmin + 1
    ncol = cmax - cmin + 1
    return rmin + nrow//2, cmin + ncol//2


def intersect(a, b):
    """Return True if two extents intersect, otherwise False

    Parameters
    ----------
    a, b : (4,) array like
        Array extents (rmin, rmax, cmin, cmax)
    
    Returns
    -------
    bool
    """

    armin, armax, acmin, acmax = a
    brmin, brmax, bcmin, bcmax = b
    return armin <= brmax and armax >= brmin and acmin <= bcmax and acmax >= bcmin


def intersection_extent(a, b):
    """Compute the extent of two overlapping extents
    
    Parameters
    ----------
    a, b : (4,) array_like
        Array extents (rmin, rmax, cmin, cmax)

    Returns
    -------
    tuple
    """

    armin, armax, acmin, acmax = a
    brmin, brmax, bcmin, bcmax = b

    rmin, rmax = max(armin, brmin), min(armax, brmax)
    cmin, cmax = max(acmin, bcmin), min(acmax, bcmax)

    return rmin, rmax, cmin, cmax


def intersection_shape(a, b):
    """Compute the shape of two overlapping extents. If there is no
    overlap, an empty tuple is returned.

    Parameters
    ----------
    a, b : (4,) array like
        Array extents (rmin, rmax, cmin, cmax)
    
    Returns
    -------
    tuple
    """

    rmin, rmax, cmin, cmax = intersection_extent(a, b)
    nr, nc = rmax - rmin + 1, cmax - cmin + 1

    if nr <= 0 or nc <= 0:
        shape = ()
    else:
        shape = (nr, nc)
    
    return shape


def intersection_slices(a, b):
    """Compute slices of overlapping areas between two overlapping extents

    Parameters
    ----------
    a, b : (4,) array like
        Array extents (rmin, rmax, cmin, cmax)
    
    Returns
    -------
    tuples of slices
    """

    rmin, rmax, cmin, cmax = intersection_extent(a, b)

    armin, armax, acmin, acmax = a
    brmin, brmax, bcmin, bcmax = b

    arow = slice(rmin-armin, rmax-armin+1)
    acol = slice(cmin-acmin, cmax-acmin+1)
    brow = slice(rmin-brmin, rmax-brmin+1)
    bcol = slice(cmin-bcmin, cmax-bcmin+1)

    return (arow, acol), (brow, bcol)


def intersection_shift(a, b):
    """Compute the shift between two overlapping extents

    Parameters
    ----------
    a, b : (4,) array like
        Array extents (rmin, rmax, cmin, cmax)
    
    Returns
    -------
    tuple
    """
    
    rmin, rmax, cmin, cmax = intersection_extent(a, b)
    nrow = rmax - rmin + 1
    ncol = cmax - cmin + 1
    return rmin + nrow//2, cmin + ncol//2