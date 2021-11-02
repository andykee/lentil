from itertools import combinations
import numpy as np

import lentil.helper


class Field:
    """
    Two-dimensional discretely sampled complex field.

    Parameters
    ----------
    data : array_like
        Array containing the sampled field. Note that real-valued data is
        cast to complex.
    pixelscale : float or None, optional
        Data array spatial sampling. If None (default), the field is assumed to
        be broadcastable to any legal shape without interpolation.
    offset : (2,) array_like of ints or None, optional.
        Shift of the Field in (row, column) from (0, 0). If None (default),
        offset = (0, 0).
    tilt : list, optional
        List of objects which implement a ``shift`` method. This method should
        accept the following parameters:

        ``shift(xs, ys, z, wavelength)``

        and return an updated x and y shift. If None (default), tilt = [].

    """
    __slots__ = ('data', 'offset', 'tilt', '_pixelscale')

    def __init__(self, data, pixelscale=None, offset=None, tilt=None):
        self.data = np.asarray(data, dtype=complex)
        self.pixelscale = pixelscale
        self.offset = offset if offset is not None else [0, 0]
        self.tilt = tilt if tilt else []

    @property
    def pixelscale(self):
        """
        Physical (row, col) sampling of each pixel in the Field.

        If ``None``, the Field does not have a pixelscale and is assumed to be
        broadcastable to any legal shape without interpolation.

        Returns
        -------
        pixelscale : ndarray or None
        """
        return self._pixelscale

    @pixelscale.setter
    def pixelscale(self, value):
        if value is not None:
            self._pixelscale = lentil.helper.sanitize_shape(value)
        else:
            self._pixelscale = None

    @property
    def shape(self):
        """
        Tuple of Field dimensions

        Returns
        -------
        shape : tuple
        """
        return self.data.shape

    @property
    def ndim(self):
        """
        Number of Field dimensions

        Returns
        -------
        ndim : int
        """
        return self.data.ndim

    @property
    def size(self):
        """
        Number of elements in the Field

        Returns
        -------
        size : int
        """
        return self.data.size

    @property
    def extent(self):
        """
        The array indices defining the extent of the shifted Field.

        Returns
        -------
        rmin, rmax, cmin, cmax : int

        Notes
        -----
        To use the values returned by ``extent()`` in a slice, ``rmax`` and ``cmax``
        should be increased by 1.
        """
        return extent(self.shape, self.offset)

    def __mul__(self, other):
        return multiply(self, other)

    def overlap(self, other):
        """
        Returns True if two fields overlap

        Parameters
        ----------
        other : `~lentil.Field`
            Other Field to check for overlap.

        Returns
        -------
        overlap : bool
        """
        return overlap(self, other)

    def shift(self, z, wavelength, pixelscale, oversample, indexing='ij'):
        """
        Compute Field shift due to associated :class:`~lentil.Tilt` objects

        Parameters
        ----------
        z : float
            Propagation distance
        wavelength : float
            Wavelength
        pixelscale : (2,) float
            Output plane spatial sampling
        oversample : int
            Output plane oversampling factor
        indexing : {'ij', 'xy'}, optional
            Cartesian ('xy') or matrix ('ij', default) indexing of output.

        Returns
        -------
        shift : tuple of floats
            Shifted location of the center of the field

        Notes
        -----
        This function is capable of doing the conversion between the (x, y) coordinates
        of a Tilt and the (r, c) coordinates of a Field.

        """
        if indexing not in ('xy', 'ij'):
            raise ValueError("Valid values for `indexing` are 'xy' and 'ij'")

        if pixelscale is None:
            raise ValueError('pixelscale must be defined to compute shift')

        x, y = 0, 0

        for tilt in self.tilt:
            x, y = tilt.shift(xs=x, ys=y, z=z, wavelength=wavelength)

        pixelscale = lentil.helper.sanitize_shape(pixelscale)
        out = x/pixelscale[0] * oversample, y/pixelscale[1] * oversample

        if indexing == 'ij':
            out = -out[1], out[0]

        return out


def extent(shape, offset):
    """
    Compute shifted array extent

    Parameters
    ----------
    shape : array_like
        Array shape
    offset : array_like
        Array offset

    Returns
    -------
    rmin, rmax, cmin, cmax : int
        Indices that define the span of the shifted array.

    Notes
    -----
    To use the values returned by ``extent()`` in a slice, ``rmax`` and ``cmax``
    should be increased by 1.

    See Also
    --------
    lentil.field.boundary : compute the extent around a number of Fields

    """
    shape = np.asarray(shape)
    offset = np.asarray(offset)

    if np.prod(shape) == 1:
        shape = np.array([1, 1])

    rmin, cmin = (-(shape//2) + offset).astype(int)
    rmax, cmax = ([rmin, cmin] + shape - 1).astype(int)

    return rmin, rmax, cmin, cmax


def overlap(a, b):
    """
    True if two Fields overlap, otherwise False

    Parameters
    ----------
    a, b : `~lentil.Field`
        Input fields.

    Returns
    -------
    overlap : bool
    """
    armin, armax, acmin, acmax = a.extent
    brmin, brmax, bcmin, bcmax = b.extent
    return armin <= brmax and armax >= brmin and acmin <= bcmax and acmax >= bcmin


def boundary(*fields):
    """
    Compute the bounding extents around a number of Field objects.

    Parameters
    ----------
    fields : :class:`~lentil.field.Field` objects
        Fields to compute boundary around

    Returns
    -------
    rmin, rmax, cmin, cmax : int
        Indices that define the extent of the shifted arrays.

    Notes
    -----
    To use the values returned by ``boundary()`` in a slice, ``rmax`` and ``cmax``
    should be increased by 1.

    See Also
    --------
    lentil.field.extent : Compute the extent of a Field

    """
    rmin, rmax, cmin, cmax = [], [], [], []

    for field in fields:
        frmin, frmax, fcmin, fcmax = field.extent
        rmin.append(frmin)
        rmax.append(frmax)
        cmin.append(fcmin)
        cmax.append(fcmax)

    return min(rmin), max(rmax), min(cmin), max(cmax)


def insert(field, out, intensity=False, weight=1, indexing='ij'):
    """

    Parameters
    ----------
    field : :class:`~lentil.field.Field`
        Field to insert into ``out``
    out : ndarray
        Array to insert ``field`` into
    intensity : bool, optional
        If True, compute intensity of ``field`` before placing it in ``out``. Default is
        False.
    weight : float, optional
        Weight to apply to ``field`` before it is inserted into ``out``
    indexing : 'ij', 'xy', optional
        Cartesian ('xy') or matrix ('ij', default) indexing of output

    Returns
    -------
    out : ndarray

    """
    if indexing not in ('xy', 'ij'):
        raise ValueError("Valid values for `indexing` are 'xy' and 'ij'")

    if field.shape == out.shape and np.array_equal(field.offset, [0, 0]):
        field_slice = Ellipsis
        output_slice = Ellipsis
    else:

        output_shape = np.asarray(out.shape)
        field_shape = np.asarray(field.shape)
        field_offset = np.asarray(field.offset)

        if indexing == 'xy':
            field_offset = -field_offset[1], field_offset[0]

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

        output_slice = slice(output_top, output_bottom), slice(output_left, output_right)
        field_slice = slice(field_top, field_bottom), slice(field_left, field_right)

    if intensity:
        out[output_slice] += (np.abs(field.data[field_slice]**2) * weight)
    else:
        out[output_slice] += (field.data[field_slice] * weight)
    return out

# # TODO: update insert() above with this syntax
# def test():
#     # field slices
#     field_rmin, field_cmin = 0, 0
#     field_rmax, field_cmax = field.shape
#
#     # slices defining where field goes in out
#     out_rmin = int(out.shape[0]/2 - field.shape[0]/2 + field.offset[0])
#     out_rmax = out_rmin + field.shape[0]
#     out_cmin = int(out.shape[1]/2 - field.shape[1]/2 + field.offset[1])
#     out_cmax = out_cmin + field.shape[1]
#
#     # figure out how much of the field is in out and
#     # where to trim if it falls off an edge
#     if out_rmin < 0:
#         field_rmin = -1 * out_rmin
#         out_rmin = 0
#
#     if out_rmax > out.shape[0]:
#         field_rmax -= out_rmax - out.shape[0]
#         out_rmax = out.shape[0]
#
#     if out_cmin < 0:
#         field_cmin = -1 * out_rmin
#         out_cmin = 0
#
#     if out_cmax > out.shape[1]:
#         field_cmax -= out_cmax - out.shape[1]
#         out_cmax = out.shape[1]
#
#     if not np.any(np.array([out_rmax, field_rmax, out_cmax, field_cmax]) < 0):
#         if intensity:
#             out[out_rmin:out_rmax, out_cmin:out_cmax] = np.abs(field.data[field_rmin:field_rmax, field_cmin:field_cmax])**2
#         else:
#             out[out_rmin:out_rmax, out_cmin:out_cmax] += field.data[field_rmin:field_rmax, field_cmin:field_cmax]
#
#     return out


def merge(a, b, check_overlap=True):
    """
    Merge two Fields

    Parameters
    ----------
    a, b : :class:`~lentil.field.Field` objects
        Fields to merge
    check_overlap : bool, optional
        If True (default), a ``ValueError`` is raised if ``a`` and ``b`` do not overlap

    Returns
    -------
    out : :class:`~lentil.field.Field`
        New ``Field`` that represents the union of ``a`` and ``b``
    """
    if not np.all(a.pixelscale == b.pixelscale):
        raise ValueError(f'pixelscales {a.pixelscale} and '
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
    # shape required to merge fields
    rmin, rmax, cmin, cmax = boundary(*fields)
    if not np.any([rmin, rmax, cmin, cmax]):  # evaluates to True if all zeros
        return ()
    else:
        return rmax - rmin + 1, cmax - cmin + 1


def _merge_slices(*fields):
    # slices in the merged array where fields go
    rmin, rmax, cmin, cmax = boundary(*fields)
    out = []
    if not np.any([rmin, rmax, cmin, cmax]):  # evaluates to True if all zeros
        out.append(Ellipsis)
    else:
        for field in fields:
            frmin, frmax, fcmin, fcmax = field.extent
            row = slice(frmin-rmin, frmax-rmin+1)
            col = slice(fcmin-cmin, fcmax-cmin+1)
            out.append((row, col))
    return out


def _merge_offset(*fields):
    # common offset of the resulting merged array
    rmin, rmax, cmin, cmax = boundary(*fields)
    nrow = rmax - rmin + 1
    ncol = cmax - cmin + 1
    return rmin + nrow//2, cmin + ncol//2


def multiply(x1, x2):
    """
    Multiply fields element-wise.

    Parameters
    ----------
    x1, x2 : :class:`~lentil.field.Field`
        Fields to be multiplied. If ``x1.shape != x2.shape``, they must
        be broadcastable to a common shape (which becomes the shape of the
        output).

    Returns
    -------
    y : :class:`~lentil.field.Field`
        The product of ``x1`` and ``x2``, appropriately sized and offset.

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
    a_data, b_data = x1.data, x2.data
    a_offset, b_offset = x1.offset, x2.offset
    pixelscale = multiply_pixelscale(x1.pixelscale, x2.pixelscale)
    tilt = x1.tilt + x2.tilt

    if a_data.size == 1 and b_data.size == 1:
        data, offset = _mul_scalar(a_data, a_offset, b_data, b_offset)
    else:
        a_data, a_offset, b_data, b_offset = _broadcast_arrays(a_data, a_offset,
                                                               b_data, b_offset)
        data, offset = _mul_array(a_data, a_offset, b_data, b_offset)

    return Field(data=data, pixelscale=pixelscale, offset=offset, tilt=tilt)


def multiply_pixelscale(a_pixelscale, b_pixelscale):
    # pixelscale reduction for multiplication
    if np.all(a_pixelscale == b_pixelscale):
        out = a_pixelscale
    else:
        if a_pixelscale is None:
            out = b_pixelscale
        elif b_pixelscale is None:
            out = a_pixelscale
        else:
            if not np.all(a_pixelscale == b_pixelscale):
                raise ValueError(f"can't multiply with inconsistent pixelscales: {a_pixelscale} != {b_pixelscale}")
            out = a_pixelscale

    return out


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
    # possibly overlapping fields that aren't combined until the data attribute
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


def reduce(*fields):
    """
    Reduce a number of disjoint Fields into a potentially smaller set where overlapping
    Fields are merged.

    Parameters
    ----------
    fields :class:`~lentil.field.Field` or :class:`~lentil.field.NDField` objects
        Fields to reduce

    Returns
    -------
    fields : list
        List of :class:`~lentil.field.NDField` objects
    """
    fields = [NDField(field) if not isinstance(field, NDField) else field for field in fields]
    for m, n in combinations(range(len(fields)), 2):
        if overlap(fields[m], fields[n]):
            fields[m].append(fields[n])
            fields.pop(n)
            return reduce(*fields)

    return fields
