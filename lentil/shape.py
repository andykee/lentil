import numpy as np

import lentil


def circle(shape, radius, shift=(0, 0), antialias=True):
    """Draw a circle

    Parameters
    ----------
    shape : array_like
        Size of output in pixels (nrows, ncols)
    radius : float
        Radius of circle in pixels
    shift : (2,) array_like, optional
        How far to shift center in float (rows, cols). Default is (0, 0).
    antialias : bool, optional
        If True (default), the shape edges are antialiased.

    Returns
    -------
    ndarray

    """
    rr, cc = lentil.helper.mesh(shape)
    r = np.sqrt(np.square(rr - shift[0]) + np.square(cc - shift[1]))
    mask = np.clip(radius + 0.5 - r, 0.0, 1.0)
    if not antialias:
        mask[mask > 0] = 1
    return mask


def hexagon(shape, radius, shift=(0, 0), rotate=False, antialias=True):
    """Draw a hexagon

    Parameters
    ----------
    shape : array_like
        Size of output in pixels (nrows, ncols)
    radius : float
        Radius of outscribing circle (which also equals the side length) in
        pixels.
    shift : tuple of floats, optional
        How far to shift center in (rows, cols). Default is (0, 0).
    rotate : bool, optional
        Rotate mask so that flat sides are aligned with the Y direction instead
        of the default orientation which is aligned with the X direction.
    antialias : bool, optional
        If True (default), the shape edges are antialiased.
        
    Returns
    -------
    ndarray

    """
    shape = np.broadcast_to(shape, (2,))
    r, c = lentil.helper.mesh(shape, shift)
    mask = np.ones(shape)

    inner_radius = radius * np.sqrt(3)/2

    for n in range(6):
    
        theta = n * np.pi/3 if rotate else n * np.pi/3 + np.pi/6
        rho = r * np.sin(theta) + c * np.cos(theta)
    
        if antialias:
            slc = np.clip(inner_radius + 0.5 - rho, 0.0, 1.0)
        else:
            slc = np.ones(shape)
            slc[rho > inner_radius] = 0
    
        mask = np.minimum(mask, slc)

    return mask


def rectangle(shape, width, height, shift=(0,0), antialias=True):
    """Draw a rectangle

    Parameters
    ----------
    shape : array_like
        Size of output in pixels (nrows, ncols)
    width : float
        Width of rectangle in pixels
    height : float
        Height of rectangle in pixels
    shift : tuple of floats, optional
        How far to shift center in (rows, cols). Default is (0, 0).
    antialias : bool, optional
        If True (default), the shape edges are antialiased.
        
    Returns
    -------
    ndarray
    
    """
    shape = np.broadcast_to(shape, (2,))
    rr, cc = lentil.helper.mesh(shape, shift)
    rect = np.ones(shape)

    width_clip = np.clip(0.5 + (width/2) - np.abs(cc), 0, 1)
    height_clip = np.clip(0.5 + (height/2) - np.abs(rr), 0, 1)

    rect = np.minimum(np.minimum(rect, width_clip), height_clip)

    if not antialias:
        rect[rect > 0] = 1

    return rect
