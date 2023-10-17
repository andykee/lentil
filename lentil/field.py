import sys
from itertools import combinations
import numpy as np

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

    Attributes
    ----------
    extent : tuple
        Array indices defining the extent of the offset Field.
    """
    __slots__ = ('data', 'offset', 'tilt', 'pixelscale', 'extent')

    def __init__(self, data, pixelscale=None, offset=None, tilt=None):
        self.data = np.asarray(data, dtype=complex)
        self.pixelscale = pixelscale
        self.offset = offset if offset is not None else [0, 0]
        self.tilt = tilt if tilt else []
        self.extent = extent(self.shape, self.offset)

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
    def size(self):
        """
        Number of elements in the Field

        Returns
        -------
        size : int
        """
        return self.data.size

    def __mul__(self, other):
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
            * A single scalar operand is broadcast to the larger operand's shape and
            inherits the larger operand's offset. This enables default planar
            fields (``data = 1+0j``, offset = [0, 0]) to automatically grow in
            size when multiplied by a "real" field.
        """
        tilt = self.tilt + other.tilt

        if self.size == 1 and other.size == 1:
            data, offset = self._mul_scalar(other)
        else:
            # Note that _mul_array is optimized to also handle scalar * array
            data, offset = self._mul_array(other)

        return Field(data=data, offset=offset, tilt=tilt)
    
    def _mul_scalar(self, other):
        """
        Multiply two Fields when both are scalars
        """
        if np.array_equal(self.offset, other.offset):
            data = self.data * other.data
            offset = self.offset
        else:
            data = 0
            offset = [0, 0]
        return data, offset
    
    def _mul_array(self, other):
        """
        Multiply two fields when at least one is an array
        """
        self_data, self_offset, other_data, other_offset = _mul_broadcast(
            self.data, self.offset, other.data, other.offset
        )
        self_extent = extent(self_data.shape, self_offset)
        other_extent = extent(other_data.shape, other_offset)

        self_slice, other_slice = _mul_slices(self_extent, other_extent)
        offset = _mul_offset(self_extent, other_extent)
        data = self_data[self_slice] * other_data[other_slice]

        if data.size == 0:
            data = np.array(0, dtype=complex)
            offset = [0, 0]
        return data, offset

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

        pixelscale = np.broadcast_to(pixelscale, (2,))
        out = x/pixelscale[0] * oversample, y/pixelscale[1] * oversample

        if indexing == 'ij':
            out = -out[1], out[0]

        return out

def boundary(fields):
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
    rmin, rmax, cmin, cmax = sys.maxsize, 0, sys.maxsize, 0

    for field in fields:
        frmin, frmax, fcmin, fcmax = field.extent
        rmin = frmin if frmin < rmin else rmin
        rmax = frmax if frmax > rmax else rmax
        cmin = fcmin if fcmin < cmin else cmin
        cmax = fcmax if fcmax > cmax else cmax

    return rmin, rmax, cmin, cmax


def extent(shape, offset):
    """
    Compute the extent of a shifted array.

    Note: To use the values returned by ``extent()`` in a slice, 
    ``rmax`` and ``cmax`` should be increased by 1.
    """
    if len(shape) < 2:
        shape = (1, 1)

    rmin = int(-(shape[0]//2) + offset[0])
    cmin = int(-(shape[1]//2) + offset[1])
    rmax = int(rmin + shape[0] - 1)
    cmax = int(cmin + shape[1] - 1)

    return rmin, rmax, cmin, cmax


def insert(field, out, intensity=False, weight=1):
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
    
    Returns
    -------
    out : ndarray

    """
    #if indexing not in ('xy', 'ij'):
    #    raise ValueError("Valid values for `indexing` are 'xy' and 'ij'")

    if field.shape == out.shape and np.array_equal(field.offset, [0, 0]):
        field_slice = Ellipsis
        out_slice = Ellipsis
    else:

        out_shape = np.asarray(out.shape)
        field_shape = np.asarray(field.shape)
        field_offset = np.asarray(field.offset)

        #if indexing == 'xy':
        #    field_offset = -field_offset[1], field_offset[0]

        # Output coordinates of the upper left corner of the shifted data array
        field_shifted_ul = (out_shape // 2) - (field_shape // 2) + field_offset

        # Field slice indices
        field_rmin = int(0)
        field_rmax = int(field_shape[0])
        field_cmin = int(0)
        field_cmax = int(field_shape[1])

        # Output insertion slice indices
        out_rmin = int(field_shifted_ul[0])
        out_rmax = int(field_shifted_ul[0] + field_shape[0])
        out_cmin = int(field_shifted_ul[1])
        out_cmax = int(field_shifted_ul[1] + field_shape[1])

        # reconcile the field and output insertion indices
        if out_rmin < 0:
            field_rmin = -1 * out_rmin
            out_rmin = 0

        if out_rmax > out_shape[0]:
            field_rmax -= out_rmax - out_shape[0]
            out_rmax = out_shape[0]

        if out_cmin < 0:
            field_cmin = -1 * out_cmin
            out_cmin = 0

        if out_cmax > out_shape[1]:
            field_cmax -= out_cmax - out_shape[1]
            out_cmax = out_shape[1]

        out_slice = slice(out_rmin, out_rmax), slice(out_cmin, out_cmax)
        field_slice = slice(field_rmin, field_rmax), slice(field_cmin, field_cmax)

    if intensity:
        out[out_slice] += (np.abs(field.data[field_slice]**2) * weight)
    else:
        out[out_slice] += (field.data[field_slice] * weight)
    return out


def merge(a, b, enforce_overlap=True):
    """
    Merge two Fields

    Parameters
    ----------
    a, b : :class:`~lentil.field.Field` objects
        Fields to merge
    enforce_overlap : bool, optional
        If True (default), a ``ValueError`` is raised if ``a`` and ``b`` do not overlap

    Returns
    -------
    out : :class:`~lentil.field.Field`
        New ``Field`` that represents the union of ``a`` and ``b``
    """
    if enforce_overlap and not overlap((a, b)):
        raise ValueError("Can't merge non-overlapping fields")

    return _merge((a, b))
# TODO: throw an error if either field has tilts


def _merge(fields):
    """
    Merge fields into a new Field regardless of whether the supplied fields 
    overlap in space
    """
    if not np.all([f.pixelscale == fields[0].pixelscale for f in fields]):
        raise ValueError("Can't merge: pixelscales must be equal")

    out = np.zeros(_merge_shape(fields), dtype=complex)
    slices = _merge_slices(fields)
    for field, slc in zip(fields, slices):
        out[slc] += field.data

    return Field(data=out,
                 pixelscale=fields[0].pixelscale,
                 offset=_merge_offset(fields))


def _merge_shape(fields):
    """
    Return the shape required to hold merged fields
    """
    rmin, rmax, cmin, cmax = boundary(fields)
    # faster than np.any([rmin, rmax, cmin, cmax])
    if rmin == 0 and rmax == 0 and cmin == 0 and cmax == 0:
        return ()
    else:
        return rmax - rmin + 1, cmax - cmin + 1


def _merge_slices(fields):
    """
    Return the slices in a merged array for field insertion
    """
    rmin, rmax, cmin, cmax = boundary(fields)
    out = []
    # faster than np.any([rmin, rmax, cmin, cmax])
    if rmin == 0 and rmax == 0 and cmin == 0 and cmax == 0:
        out.append(Ellipsis)
    else:
        for field in fields:
            frmin, frmax, fcmin, fcmax = field.extent
            row = slice(frmin-rmin, frmax-rmin+1)
            col = slice(fcmin-cmin, fcmax-cmin+1)
            out.append((row, col))
    return out


def _merge_offset(fields):
    """
    Return the new offset of a merged field
    """
    rmin, rmax, cmin, cmax = boundary(fields)
    nrow = rmax - rmin + 1
    ncol = cmax - cmin + 1
    return rmin + nrow//2, cmin + ncol//2


def overlap(fields):
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
    #return _overlap(a.extent, b.extent)

    if len(fields) == 2:
        return _overlap(fields[0].extent, fields[1].extent)
    else:
        fields = _reduce(fields)
        if len(fields) > 1:
            return False
        else:
            return True

    
def reduce(fields):
    """
    Reduce a number of Fields into a potentially smaller disjoint set where 
    overlapping Fields are merged.

    Parameters
    ----------
    fields :class:`~lentil.field.Field` objects
        Fields to reduce

    Returns
    -------
    fields : list
        List of :class:`~lentil.field.Field` objects
    
    """
    fields = _reduce(fields)

    out = []

    for f in fields:
        if len(f['field']) > 1:
            out.append(_merge(f['field']))
        else:
            out.append(f['field'][0])

    return out


def _reduce(fields):
    """
    Function to identify overlapping fields and reorganize them into
    a dict for easier merging.
    """
    fields = [{'field': [f], 'extent': f.extent} for f in fields]
    return _disjoint(fields)


def _disjoint(fields):
    """
    Return fields as a disjoint set.
    """
    for m, n in combinations(range(len(fields)), 2):
        if _overlap(fields[m]['extent'], fields[n]['extent']):
            fields[m]['field'].extend(fields[n]['field'])
            fields[m]['extent'] = boundary(fields[m]['field'])
            fields.pop(n)
            return _disjoint(fields)
    return fields


def _overlap(a_extent, b_extent):
    """
    Return True if two extents overlap, otherwise False
    """
    armin, armax, acmin, acmax = a_extent
    brmin, brmax, bcmin, bcmax = b_extent
    return armin <= brmax and armax >= brmin and acmin <= bcmax and acmax >= bcmin


def _mul_broadcast(a_data, a_offset, b_data, b_offset):
    """
    Broadcast for multiplication. 

    If one array is a scalar, it is broadcast to a compatible shape with the
    other array. As a part of this operation, the broadcasted field inherits
    the larger array's offset
    """
    # if one array is a scalar, it is broadcast to a compatible shape
    # with the other array. as a part of this operation, the broadcasted
    # Field inherits the other Field's offset as well
    #a_data, a_offset = a.data, a.offset
    #b_data, b_offset = b.data, b.offset
    if a_data.shape != b_data.shape:
        if a_data.size == 1:
            a_data = np.broadcast_to(a_data, b_data.shape)
            a_offset = b_offset
        if b_data.size == 1:
            b_data = np.broadcast_to(b_data, a_data.shape)
            b_offset = a_offset
    return a_data, a_offset, b_data, b_offset


def _mul_boundary(a_extent, b_extent):
    # bounding  array indices to be multiplied
    armin, armax, acmin, acmax = a_extent
    brmin, brmax, bcmin, bcmax = b_extent

    rmin, rmax = max(armin, brmin), min(armax, brmax)
    cmin, cmax = max(acmin, bcmin), min(acmax, bcmax)

    return rmin, rmax, cmin, cmax


def _mul_slices(a_extent, b_extent):
    rmin, rmax, cmin, cmax = _mul_boundary(a_extent, b_extent)

    armin, armax, acmin, acmax = a_extent
    brmin, brmax, bcmin, bcmax = b_extent

    arow = slice(rmin-armin, rmax-armin+1)
    acol = slice(cmin-acmin, cmax-acmin+1)
    brow = slice(rmin-brmin, rmax-brmin+1)
    bcol = slice(cmin-bcmin, cmax-bcmin+1)

    return (arow, acol), (brow, bcol)


def _mul_offset(a_extent, b_extent):
    rmin, rmax, cmin, cmax = _mul_boundary(a_extent, b_extent)
    nrow = rmax - rmin + 1
    ncol = cmax - cmin + 1
    return rmin + nrow//2, cmin + ncol//2
