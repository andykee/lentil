import numpy as np

import lentil
import lentil._mesh


def circle(shape, radius, shift=(0, 0), antialias=True, indexing='xy'):
    """Draw a circle

    Parameters
    ----------
    shape : array_like
        Size of output in pixels (nrows, ncols)
    radius : float
        Radius of circle in pixels
    shift : (2,) array_like, optional
        How far to shift center of circle. Default is (0, 0).
    antialias : bool, optional
        If True (default), the shape edges are antialiased.
    indexing : {'xy', 'ij'}, optional
        Cartesian ('xy', default) or matrix ('ij') or indexing of shift.

    Returns
    -------
    ndarray

    """
    x1, x2 = lentil._mesh.mesh(shape, shift=shift, indexing=indexing)
    r = np.sqrt(x1**2 + x2**2)
    mask = np.clip(radius + 0.5 - r, 0.0, 1.0)
    if not antialias:
        mask[mask > 0] = 1
    return mask


def hexagon(shape, radius, shift=(0, 0), rotate=True, antialias=True,
            indexing='xy'):
    """Draw a hexagon

    Parameters
    ----------
    shape : array_like
        Size of output in pixels (nrows, ncols)
    radius : float
        Radius of outscribing circle (which also equals the side length) in
        pixels.
    shift : tuple of floats, optional
        How far to shift center of hexagon. Default is (0, 0).
    rotate : bool, optional
        If True (default), the mask is rotated so that flat sides are aligned
        with the X direction instead with the Y direction.
    antialias : bool, optional
        If True (default), the shape edges are antialiased.
    indexing : {'xy', 'ij'}, optional
            Cartesian ('xy', default) or matrix ('ij') or indexing of shift.
        
    Returns
    -------
    ndarray

    """
    shape = np.broadcast_to(shape, (2,))
    r, c = lentil._mesh.mesh(shape, shift=shift, indexing=indexing)
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


def rectangle(shape, width, height, shift=(0,0), angle=0, antialias=True,
              indexing='xy'):
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
        How far to shift center of rectangle. Default is (0, 0).
    angle : float, optional
        Rotation of rectangle in degrees counterclockwise from horizontal.
        Default is 0.
    antialias : bool, optional
        If True (default), the shape edges are antialiased.
    indexing : {'xy', 'ij'}, optional
            Cartesian ('xy', default) or matrix ('ij') or indexing of shift.
        
    Returns
    -------
    ndarray
    
    """
    shape = np.broadcast_to(shape, (2,))
    rr, cc = lentil._mesh.mesh(shape, shift=shift, angle=angle,
                               indexing=indexing)
    rect = np.ones(shape)

    width_clip = np.clip(0.5 + (width/2) - np.abs(cc), 0, 1)
    height_clip = np.clip(0.5 + (height/2) - np.abs(rr), 0, 1)

    rect = np.minimum(np.minimum(rect, width_clip), height_clip)

    if not antialias:
        rect[rect > 0] = 1

    return rect


def spider(shape, width, angle=0, shift=(0,0), antialias=True,
           indexing='xy'):
    """Draw a spider

    Parameters
    ----------
    shape : array_like
        Size of output in pixels (nrows, ncols)
    width : float
        Width of rectangle in pixels
    angle : float, optional
        Rotation of spider in degrees counterclockwise from horizontal.
        Default is 0.
    shift : tuple of floats, optional
        How far to shift center of spider. Default is (0, 0).
    antialias : bool, optional
        If True (default), the spider edges are antialiased.
    indexing : {'xy', 'ij'}, optional
            Cartesian ('xy', default) or matrix ('ij') or indexing of shift.
    
    Returns
    -------
    ndarray

    """
    shape = np.broadcast_to(shape, (2,))
    len = 2*np.max(shape)
    shift_dist = len / 2
    shift_row = -shift_dist * np.sin(np.deg2rad(angle))
    shift_col = shift_dist * np.cos(np.deg2rad(angle))
    if indexing == 'xy':
        shift = [-shift[1], shift[0]]
    shift = (shift[0] + shift_row, shift[1] + shift_col)
    return 1 - lentil.rectangle(shape, width=len, height=width, shift=shift,
                                angle=angle, antialias=antialias, 
                                indexing='ij')
