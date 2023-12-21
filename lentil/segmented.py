import collections

import numpy as np

import lentil


def hex_segments(rings, seg_radius, seg_gap, rotate=False, antialias=True,
                 flatten=False, pad=2, drop=(0,)):
    """Draw a segmented aperture made up of hexagonal segments

    Parameters
    ----------
    rings : int
        Number of segment rings. Must be >= 1.
    seg_radius : float
        Segment outscribing radius in pixels
    seg_gap : float
        Spacing between adjacent segments in pixels
    rotate : bool, optional
        Rotate segments so that flat sides are aligned with the Y direction 
        instead of the default orientation which is aligned with the X 
        direction.
    antialias : bool, optional
        If True (default), the segment edges are antialiased.
    flatten : bool, optional
        If True, the individual segment masks are flattened into a single
        global mask. Default is False.
    pad : int, optional
        Number of additional pixels to zero-pad output by. Default is 2.
    drop : list_like, optional
        Segments to exclude from output. 0 represents the central segment. 
    
    Returns
    -------
    ndarray

    Examples
    --------
    .. plot::
        :include-source:
        :context: reset
        :scale: 50

        >>> import matplotlib.pyplot as plt
        >>> import lentil
        >>> mask = lentil.hex_segments(rings=2, seg_radius=64, seg_gap=2)
        >>> fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(5, 2))
        >>> ax[0].imshow(np.sum(mask, axis=0))
        >>> ax[0].set_title('Full aperture')
        >>> ax[1].imshow(mask[0])
        >>> ax[1].set_title('Segment 1')

    References
    ----------
    [1] Hexagonal Grids (https://www.redblobgames.com/grids/hexagons/)

    """
    inner_radius = seg_radius * np.sqrt(3)/2
    size = np.ceil((rings * 2 + 1) * inner_radius * 2 + (rings * 2) * seg_gap + pad * 2).astype(int)
    shape = np.broadcast_to(size, (2,))

    mask = []

    if 0 not in drop:
        mask.append(lentil.hexagon(shape, seg_radius, shift=(0, 0), antialias=antialias, rotate=rotate))

    seg = 1
    for ring in range(1, rings+1):
        for h in hex_ring(ring):
            r, c = hex_to_rc(h, seg_radius + seg_gap/2, rotate)
            if seg not in drop:
                mask.append(lentil.hexagon(shape, seg_radius, shift=(r, c), antialias=antialias, rotate=rotate))
            seg += 1
    
    mask = np.asarray(mask)
    if flatten:
        mask = np.sum(mask, axis=0)
    
    return mask


# the functions below for working with hex grids was copied almost
# exactly from http://www.redblobgames.com/grids/hexagons/

_Hex = collections.namedtuple("Hex", ["q", "r", "s"])
def Hex(q, r, s):
    assert not (round(q + r + s) != 0), "q + r + s must be 0"
    return _Hex(q, r, s)


def hex_add(a, b):
    return Hex(a.q + b.q, a.r + b.r, a.s + b.s)


hex_directions = [Hex(1, 0, -1), Hex(1, -1, 0), Hex(0, -1, 1), 
                  Hex(-1, 0, 1), Hex(-1, 1, 0), Hex(0, 1, -1)]


def hex_direction(direction):
    return hex_directions[direction]


def hex_neighbor(hex, direction):
    return hex_add(hex, hex_direction(direction))


def hex_ring(radius):
    results = []
    hex = Hex(-radius, radius, 0)

    for i in range(6):
        for j in range(radius):
            results.append(hex)
            hex = hex_neighbor(hex, i)
    return results


def hex_to_xy(hex, radius, rotate=False):
    if rotate:
        x = radius * (np.sqrt(3) * hex.q + np.sqrt(3)/2 * hex.r)
        y = radius * (3/2 * hex.r)
    else:
        x = radius * (3/2 * hex.q)
        y = radius * (np.sqrt(3)/2 * hex.q + np.sqrt(3) * hex.r)
        
    return x, y


def hex_to_rc(hex, radius, rotate=False):
    x, y = hex_to_xy(hex, radius, rotate)
    return -y, x

