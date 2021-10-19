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

    __slots__ = ('data', 'pixelscale', 'offset', 'tilt')

    def __init__(self, data, pixelscale, offset=None, tilt=None):
        self.data = np.asarray(data, dtype=complex)
        self.pixelscale = pixelscale
        self.offset = offset if offset is not None else [0, 0]
        self.tilt = tilt if tilt else []

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

    def shift(self, z, wavelength, pixelscale, oversample, indexing='xy'):
        """Compute Field shift.
        The center may be shifted due to objects that implement the
        :class:`~lentil.Tilt` interface in Field.tilt.
        Parameters
        ----------
        z : float
            Propagation distance
        wavelength : float
            Wavelength
        pixelscale : float
            Output plane spatial sampling
        Returns
        -------
        x, y : float
            (x,y) location of the center of the field
        """
        x, y = 0, 0

        for tilt in self.tilt:
            x, y = tilt.shift(xs=x, ys=y, z=z, wavelength=wavelength)

        if indexing == 'xy':
            return x/pixelscale*oversample, y/pixelscale*oversample
        elif indexing == 'ij':
            return -y/pixelscale*oversample, x/pixelscale*oversample


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


def insert(field, out, intensity=False, indexing='xy'):
    output_shape = np.asarray(out.shape)
    field_shape = np.asarray(field.shape)

    if indexing == 'xy':
        field_offset = np.flip(field.offset)
        # field_offset[0] = -field_offset[0]
    elif indexing == 'ij':
        field_offset = np.asarray(field.offset)
    else:
        raise ValueError(f"Unknown indexing {indexing}. indexing must be 'ij' or 'xy'.")

    # Output coordinates of the upper left corner of the shifted data array
    field_shifted_ul = (output_shape // 2) - (field_shape // 2) + field_offset

    # Field slice indices
    field_top = int(0)
    field_bottom = int(field_shape[0])
    field_left = int(0)
    field_right = int(field_shape[1])

    # Output insertion slice indices
    output_top = int(field_shifted_ul[0])
    output_bottom = int(field_shifted_ul[0] + field_shape[0])
    output_left = int(field_shifted_ul[1])
    output_right = int(field_shifted_ul[1] + field_shape[1])

    # reconcile the field and output insertion indices
    if output_top < 0:
        field_top = -1 * output_top
        output_top = 0

    if output_bottom > output_shape[0]:
        field_bottom -= output_bottom - output_shape[0]
        output_bottom = output_shape[0]

    if output_left < 0:
        field_left = -1 * output_left
        output_left = 0

    if output_right > output_shape[1]:
        field_right -= output_right - output_shape[1]
        output_right = output_shape[1]

    if intensity:
        out[output_top:output_bottom, output_left:output_right] = np.abs(field.data[field_top:field_bottom, field_left:field_right]**2)
    else:
        out[output_top:output_bottom, output_left:output_right] += field.data[field_top:field_bottom, field_left:field_right]
    return out



# TODO: update insert() above with this syntax
def test():
    # field slices
    field_rmin, field_cmin = 0, 0
    field_rmax, field_cmax = field.shape

    # slices defining where field goes in out
    out_rmin = int(out.shape[0]/2 - field.shape[0]/2 + field.offset[0])
    out_rmax = out_rmin + field.shape[0]
    out_cmin = int(out.shape[1]/2 - field.shape[1]/2 + field.offset[1])
    out_cmax = out_cmin + field.shape[1]

    # figure out how much of the field is in out and
    # where to trim if it falls off an edge
    if out_rmin < 0:
        field_rmin = -1 * out_rmin
        out_rmin = 0

    if out_rmax > out.shape[0]:
        field_rmax -= out_rmax - out.shape[0]
        out_rmax = out.shape[0]

    if out_cmin < 0:
        field_cmin = -1 * out_rmin
        out_cmin = 0

    if out_cmax > out.shape[1]:
        field_cmax -= out_cmax - out.shape[1]
        out_cmax = out.shape[1]

    if not np.any(np.array([out_rmax, field_rmax, out_cmax, field_cmax]) < 0):
        if intensity:
            out[out_rmin:out_rmax, out_cmin:out_cmax] = np.abs(field.data[field_rmin:field_rmax, field_cmin:field_cmax])**2
        else:
            out[out_rmin:out_rmax, out_cmin:out_cmax] += field.data[field_rmin:field_rmax, field_cmin:field_cmax]

    return out


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


def _merge_shape(*fields):
    # shape required to merge a and b
    rmin, rmax, cmin, cmax = boundary(*fields)
    return rmax - rmin + 1, cmax - cmin + 1


def _merge_slices(*fields):
    # slices in the merged array where a and b go
    rmin, rmax, cmin, cmax = boundary(*fields)

    out = []
    for field in fields:
        frmin, frmax, fcmin, fcmax = field.extent
        row = slice(frmin-rmin, frmax-rmin+1)
        col = slice(fcmin-cmin, fcmax-cmin+1)
        out.append((row, col))
    return out


def _merge_offset(*fields):
    # offset of the resulting merged array
    rmin, rmax, cmin, cmax = boundary(*fields)
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
    # bounding  array indices to be multiplied
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


class NDField:
    __slots__ = 'fields'

    # NDField looks like a standard field but may actually contain multiple
    # possibly overlapping fiuelds that aren't combined until the data attribute
    # is accessed
    def __init__(self, field):
        self.fields = [field]

    @property
    def data(self):
        # can use _merge shape and merge slices
        out = np.zeros(self.shape, dtype=complex)
        slices = _merge_slices(*self.fields)
        for field, slc in zip(self.fields, slices):
            out[slc] += field.data
        return out

    @property
    def offset(self):
        return _merge_offset(*self.fields)

    @property
    def extent(self):
        return boundary(*self.fields)

    @property
    def shape(self):
        return _merge_shape(*self.fields)

    def append(self, field):
        self.fields.append(field)

    def overlap(self, other):
        return overlap(self, other)
