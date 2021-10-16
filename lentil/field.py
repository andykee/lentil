import copy
import numpy as np

class Field:
    """A Field object represents a two-dimensional discretely sampled complex
    field.

    Parameters
    ----------
    data : array_like
        Array containing the sampled field. Note that real-valued data is
        cast to complex.
    pixelscale : float
        Data array spatial sampling
    offset : (2,) array_like of ints
        Shift in (row, column) from (0, 0)

    Attributes
    ----------
    shape : tuple of inte
        Tuple of Field dimensions
    extent : typle of ints
        Tuple of Field extent given as (row min, row max, col min, col max)

    """
    
    __slots__ = ('data', 'pixelscale', 'offset')

    def __init__(self, data, pixelscale, offset=None):
        self.data = np.asarray(data, dtype=complex)
        self.pixelscale = pixelscale
        self.offset = offset if offset is not None else [0,0]

    @property
    def shape(self):
        return self.data.shape

    @property
    def ndim(self):
        return self.data.ndim

    @property
    def size(self):
        return self.data.size

    @property
    def extent(self):
        return extent(self.shape, self.offset)

    def __mul__(self, other):
        return multiply(self, other)

    def overlap(self, other):
        """Returns True if two fields overlap

        Parameters
        ----------
        other : `~lentil.Field`
            Other field to check for overlap.

        Returns
        -------
        overlap : bool
            Returns True if the fields overlap.
        
        """
        return overlap(self, other)


def extent(shape, offset):
    shape = np.asarray(shape)
    offset = np.asarray(offset)

    if np.prod(shape) == 1:
        shape = np.array([1, 1])

    rmin, cmin = (-(shape//2) + offset).astype(int)
    rmax, cmax = ([rmin, cmin] + shape - 1).astype(int)

    return rmin, rmax, cmin, cmax


def overlap(a, b):
    """True if two Fields overlap, otherwise False

    Parameters
    ----------
    a1, a2 : `~lentil.Field`
        Input fields.

    Returns
    -------
    overlap : bool
        Returns True of the fields overlap.

    """
    armin, armax, acmin, acmax = a.extent
    brmin, brmax, bcmin, bcmax = b.extent
    return armin <= brmax and armax >= brmin and acmin <= bcmax and acmax >= bcmin


def boundary(*fields):
    rmin, rmax, cmin, cmax = [], [], [], []

    for field in fields:
        frmin, frmax, fcmin, fcmax = field.extent
        rmin.append(frmin)
        rmax.append(frmax)
        cmin.append(fcmin)
        cmax.append(fcmax)

    return min(rmin), max(rmax), min(cmin), max(cmax)


def merge(a, b, check_overlap=True):
    if a.pixelscale != b.pixelscale:
        raise ValueError(f'pixelscales {a.pixelscale} and '\
                         f'{b.pixelscale} are not equal')
    if check_overlap and not overlap(a, b):
        raise ValueError("can't merge non-overlapping fields")

    out = np.zeros(_merge_shape(a, b), dtype=complex)
    a_slc, b_slc = _merge_slices(a, b)
    out[a_slc] = a.data
    out[b_slc] += b.data
    return Field(data=out, pixelscale=a.pixelscale,
                 offset=_merge_offset(a, b))


def _merge_shape(a, b):
    # shape required to merge a and b
    rmin, rmax, cmin, cmax = boundary(a, b)
    return rmax - rmin + 1, cmax - cmin + 1


def _merge_slices(a, b):
    # slices in the merged array where a and b go
    armin, armax, acmin, acmax = a.extent
    brmin, brmax, bcmin, bcmax = b.extent
    rmin, rmax, cmin, cmax = boundary(a, b)

    arow = slice(armin-rmin, armax-rmin+1)
    acol = slice(acmin-cmin, acmax-cmin+1)
    brow = slice(brmin-rmin, brmax-rmin+1)
    bcol = slice(bcmin-cmin, bcmax-cmin+1)

    return (arow, acol), (brow, bcol)

def _merge_offset(a, b):
    # offset of the resulting merged array
    rmin, rmax, cmin, cmax = boundary(a, b)
    nrow = rmax - rmin + 1
    ncol = cmax - cmin + 1
    return rmin + nrow//2, cmin + ncol//2


def multiply(x1, x2):
    """Multiply fields element-wise.
    
    Parameters
    ----------
    x1, x2 : `~lentil.Field`
        Input fields to be multiplied. If ``x1.shape != x2.shape``, they must
        be broadcastable to a common shape (which becomes the shape of the
        output).

    Returns
    -------
    y : `~lentil.Field`
        The product of `x1` and `x2`, appropriately sized and offset.

    Notes
    -----
    The output Field data and offset are computed according to the following 
    rules:
        * If both operands are scalars, the result is 0 unless the operands 
          share the same offset.
        * A single scalar operand is broadcast to the other operand's shape and
          inherits the other operand's offset. This enables default planar 
          fields (``data = 1+0j``, offset = [0, 0]) to automatically grow in 
          size when multiplied by a "real" field.

    """
    if x1.pixelscale != x2.pixelscale:
        raise ValueError(f'pixelscales {x1.pixelscale} and '\
                         f'{x2.pixelscale} are not equal')

    a_data, b_data = x1.data, x2.data
    a_offset, b_offset = x1.offset, x2.offset
    pixelscale = x1.pixelscale

    if a_data.size == 1 and b_data.size == 1:
        data, offset = _mul_scalar(a_data, a_offset, b_data, b_offset)
    else:
        a_data, a_offset, b_data, b_offset = _broadcast_arrays(a_data, a_offset, 
                                                               b_data, b_offset)
        data, offset = _mul_array(a_data, a_offset, b_data, b_offset)

    return Field(data=data, pixelscale=pixelscale, offset=offset)


def _mul_scalar(a_data, a_offset, b_data, b_offset):
    if np.array_equal(a_offset, b_offset):
        data = a_data * b_data
        offset = a_offset
    else:
        data = 0
        offset = [0, 0]
    return data, offset


def _broadcast_arrays(a_data, a_offset, b_data, b_offset):
    # if one array is a scalar, it is broadcast to a compatible shape
    # with the other array. as a part of this operation, the broadcasted
    # Field inherits the other Field's offset as well
    if a_data.shape != b_data.shape:
        if a_data.size == 1:
            a_data = np.broadcast_to(a_data, b_data.shape)
            a_offset = b_offset
        if b_data.size == 1:
            b_data = np.broadcast_to(b_data, a_data.shape)
            b_offset = a_offset
    return a_data, a_offset, b_data, b_offset


def _mul_array(a_data, a_offset, b_data, b_offset):
    a_slice, b_slice = _mul_slices(a_data.shape, a_offset, 
                                   b_data.shape, b_offset)
    offset = _mul_offset(a_data.shape, a_offset,
                         b_data.shape, b_offset)
    data = a_data[a_slice] * b_data[b_slice]

    # if a_data and b_data are arrays, but don't overlap, [] is returned
    # TODO: consider returning None instead?
    if not np.any(data):
        data = np.array(0, dtype=complex)
        offset = [0, 0]

    return data, offset


def _mul_boundary(a_shape, a_offset, b_shape, b_offset):
    #bounding  array indices to be multiplied
    armin, armax, acmin, acmax = extent(a_shape, a_offset)
    brmin, brmax, bcmin, bcmax = extent(b_shape, b_offset)

    rmin, rmax = max(armin, brmin), min(armax, brmax)
    cmin, cmax = max(acmin, bcmin), min(acmax, bcmax)

    return rmin, rmax, cmin, cmax


def _mul_slices(a_shape, a_offset, b_shape, b_offset):
    armin, armax, acmin, acmax = extent(a_shape, a_offset)
    brmin, brmax, bcmin, bcmax = extent(b_shape, b_offset)
    
    rmin, rmax = max(armin, brmin), min(armax, brmax)
    cmin, cmax = max(acmin, bcmin), min(acmax, bcmax)
    
    arow = slice(rmin-armin, rmax-armin+1)
    acol = slice(cmin-acmin, cmax-acmin+1)
    brow = slice(rmin-brmin, rmax-brmin+1)
    bcol = slice(cmin-bcmin, cmax-bcmin+1)

    return (arow, acol), (brow, bcol)


def _mul_offset(a_shape, a_offset, b_shape, b_offset):
    rmin, rmax, cmin, cmax = _mul_boundary(a_shape, a_offset, b_shape, b_offset)
    nrow = rmax - rmin + 1
    ncol = cmax - cmin + 1
    return rmin + nrow//2, cmin + ncol//2


#def _mul_boundary(a, b):
#    #bounding  array indices to be multiplied
#    armin, armax, acmin, acmax = a.extent
#    brmin, brmax, bcmin, bcmax = b.extent
#
#    rmin, rmax = max(armin, brmin), min(armax, brmax)
#    cmin, cmax = max(acmin, bcmin), min(acmax, bcmax)
#
#    return rmin, rmax, cmin, cmax

#def _mul_slices(a, b):
#    #
#    armin, armax, acmin, acmax = a.extent
#    brmin, brmax, bcmin, bcmax = b.extent
#    rmin, rmax, cmin, cmax = _mul_boundary(a, b)
#
#    arow = slice(rmin-armin, rmax-armin+1)
#    acol = slice(cmin-acmin, cmax-acmin+1)
#    brow = slice(rmin-brmin, rmax-brmin+1)
#    bcol = slice(cmin-bcmin, cmax-bcmin+1)
#
#    return (arow, acol), (brow, bcol)


#def _mul_offset(a, b):
#    rmin, rmax, cmin, cmax = _mul_boundary(a, b)
#    nrow = rmax - rmin + 1
#    ncol = cmax - cmin + 1
#    return rmin + nrow//2, cmin + ncol//2

