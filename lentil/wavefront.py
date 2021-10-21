import copy
from itertools import combinations

import numpy as np

import lentil.field
from lentil.field import Field, NDField
import lentil.fourier


class Wavefront:
    """A class representing a monochromatic wavefront. :class:`Wavefront` is
    used internally by Lentil to perform diffraction propagation calculations.

    Parameters
    ----------
    wavelength : float
        Wavelength in meters
    pixelscale : float, optional
        Wavefront array spatial sampling in meters/pixel
    data : list_like, optional
        Wavefront data. Default is [1+0j] (a plane wave).

    Attributes
    ----------
    focal_length : float or np.inf
        Wavefront focal length. A plane wave (default) has an infinite focal
        length (``np.inf``).

    tilt : list
        List of objects which implement a ``shift`` method. This method should
        accept the following parameters:

        ``shift(xs, ys, z, wavelength)``

        and return an updated x and y shift.

    """

    __slots__ = ('wavelength', 'pixelscale', 'focal_length',
                 'data', 'shape')

    def __init__(self, wavelength, pixelscale=None, data=None):

        self.wavelength = wavelength
        self.pixelscale = pixelscale

        if data is None:
            self.data = [Field(data=1, pixelscale=pixelscale, offset=[0, 0], tilt=[])]
        else:
            self.data = [*data]

        # Wavefront focal length (which is infinity for a plane wave)
        self.focal_length = np.inf

    def __mul__(self, plane):
        return plane.multiply(self)

    def __rmul__(self, other):
        return self.__mul__(other)

    @property
    def depth(self):
        """Number of Fields in :attr:`data`"""
        return len(self.data)

    @property
    def field(self):
        out = np.zeros(self.shape, dtype=complex)
        for field in self.data:
            out = lentil.field.insert(field, out)
        return out

    @property
    def intensity(self):
        out = np.zeros(self.shape, dtype=float)
        for field in lentil.field.reduce(*self.data):
            out = lentil.field.insert(field, out, intensity=True)
        return out
        #return np.abs(self.field**2)

    def copy(self):
        return copy.deepcopy(self)

    def propagate_image(self, pixelscale, npix, npix_prop=None, oversample=2):
        out = Wavefront(wavelength=self.wavelength, data=[])

        npix = _sanitize_shape(npix)
        npix_prop = _sanitize_shape(npix_prop, default=npix)

        out.shape = npix * oversample
        prop_shape = npix_prop * oversample

        for field in self.data:
            # compute the field shift from any embedded tilts. note the return value
            # is specified in terms of (r, c)
            shift = field.shift(z=self.focal_length, wavelength=self.wavelength,
                                pixelscale=pixelscale, oversample=oversample,
                                indexing='ij')

            fix_shift = np.fix(shift)
            dft_shift = shift - fix_shift

            if _overlap(prop_shape, fix_shift, out.shape):
                alpha = _dft_alpha(dx=self.pixelscale, du=pixelscale,
                                   wave=self.wavelength, z=self.focal_length,
                                   oversample=oversample)
                data = lentil.fourier.dft2(f=field.data, alpha=alpha,
                                           npix=prop_shape, shift=dft_shift,
                                           offset=field.offset, unitary=True)
                out.data.append(Field(data=data, pixelscale=pixelscale/oversample,
                                      offset=fix_shift))
        return out


def _overlap(field_shape, field_shift, output_shape):
    # Return True if there's any overlap between a shifted field and the
    # output shape
    output_shape = np.asarray(output_shape)
    field_shape = np.asarray(field_shape)
    field_shift = np.asarray(field_shift)

    # Output coordinates of the upper left corner of the shifted data array
    field_shifted_ul = (output_shape / 2) - (field_shape / 2) + field_shift

    if field_shifted_ul[0] > output_shape[0]:
        return False
    if field_shifted_ul[0] + field_shape[0] < 0:
        return False
    if field_shifted_ul[1] > output_shape[1]:
        return False
    if field_shifted_ul[1] + field_shape[1] < 0:
        return False
    return True


def _sanitize_shape(shape, default=()):
    if shape is None:
        shape = default
    shape = np.asarray(shape)
    if shape.shape == ():
        shape = np.append(shape, shape)
    return shape


def _dft_alpha(dx, du, wave, z, oversample):
    return (dx*du)/(wave*z*oversample)


class imtile:
    def __init__(self, data, global_slice, data_slice):
        self._data = [data]
        self._data_slice = [data_slice]
        self._global_slice = [global_slice]
        self._slice = global_slice

    def join(self, imtile):
        self._data = [*self._data, *imtile._data]
        self._data_slice = [*self._data_slice, *imtile._data_slice]
        self._global_slice = [*self._global_slice, *imtile._global_slice]

        # update the bounding slice
        self._slice = bounding_slice(self._slice, imtile.slice)  # the full output tile slice

    @property
    def slice(self):
        return self._slice

    @property
    def _local_slice(self):
        return [local_slice(self._slice, slc) for slc in self._global_slice]

    @property
    def data(self):
        # this creates the output array and places everything

        data = np.zeros((self.slice[0].stop - self.slice[0].start,
                         self.slice[1].stop - self.slice[1].start),
                        dtype=complex)

        # compute local slices based on slice and global slice

        _local_slice = self._local_slice

        for n in range(len(self._data)):
            data[_local_slice[n]] += self._data[n][self._data_slice[n]]

        return data


def overlap(a, b):
    return a[0].start <= b[0].stop and a[0].stop >= b[0].start and a[1].start <= b[1].stop and a[1].stop >= b[1].start


def bounding_slice(a, b):
    rmin = min(a[0].start, b[0].start)
    rmax = max(a[0].stop, b[0].stop)
    cmin = min(a[1].start, b[1].start)
    cmax = max(a[1].stop, b[1].stop)
    return slice(rmin, rmax), slice(cmin, cmax)


def local_slice(output_slice, chip_slice):
    row_slc = chip_slice[0].start - output_slice[0].start, chip_slice[0].stop - output_slice[0].start
    col_slc = chip_slice[1].start - output_slice[1].start, chip_slice[1].stop - output_slice[1].start
    return slice(*row_slc), slice(*col_slc)


def consolidate(tiles):
    for m, n in combinations(range(len(tiles)), 2):
        if overlap(tiles[m].slice, tiles[n].slice):
            tiles[m].join(tiles[n])
            tiles.pop(n)
            return consolidate(tiles)

    return tiles
