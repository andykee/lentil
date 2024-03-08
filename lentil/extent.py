# Functions for working with extent (rmon, rmax, cmin, cmax) data

import numpy as np

def array_extent(shape, shift, parent_shape=None):
    """Compute the extent of a shifted array.

    Parameters
    ----------
    shape : (2,) array_like

    shift : (2,) array_like

    parent_shape : (2,) array like or None, optional
    
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


def intersect(a, b):
    """Return True if two extents intersect, otherwise False

    Parameters
    ----------
    a, b : (4,) array like
        Two array extents (rmin, rmax, cmin, cmax)
    
    Returns
    -------
    bool
    """

    armin, armax, acmin, acmax = a
    brmin, brmax, bcmin, bcmax = b
    return armin <= brmax and armax >= brmin and acmin <= bcmax and acmax >= bcmin


def intersection_extent(a, b):
    # bounding  array indices to be multiplied
    armin, armax, acmin, acmax = a
    brmin, brmax, bcmin, bcmax = b

    rmin, rmax = max(armin, brmin), min(armax, brmax)
    cmin, cmax = max(acmin, bcmin), min(acmax, bcmax)

    return rmin, rmax, cmin, cmax


def intersection_shape(a, b):
    """Compute the shape 
    """

    rmin, rmax, cmin, cmax = intersection_extent(a, b)
    nr, nc = rmax - rmin + 1, cmax - cmin + 1

    if nr < 0 or nc < 0:
        shape = ()
    else:
        shape = (nr, nc)
    
    return shape


def intersection_slices(a, b):
    rmin, rmax, cmin, cmax = intersection_extent(a, b)

    armin, armax, acmin, acmax = a
    brmin, brmax, bcmin, bcmax = b

    arow = slice(rmin-armin, rmax-armin+1)
    acol = slice(cmin-acmin, cmax-acmin+1)
    brow = slice(rmin-brmin, rmax-brmin+1)
    bcol = slice(cmin-bcmin, cmax-bcmin+1)

    return (arow, acol), (brow, bcol)


def intersection_shift(a, b):
    rmin, rmax, cmin, cmax = intersection_extent(a, b)
    nrow = rmax - rmin + 1
    ncol = cmax - cmin + 1
    return rmin + nrow//2, cmin + ncol//2