import copy
from itertools import combinations
from operator import itemgetter

import numpy as np

from lentil import util
from lentil import fourier
from lentil.plane import Plane, Pupil, Image, Tilt
from lentil.wavefront import Wavefront

__all__ = ['propagate_image']


def propagate_image(wavefront, pixelscale, npix, npix_prop=None, 
                    oversample=2, place_tiles=True):
    """Propagate a :class:`~lentil.Wavefront` from a :class:`~lentil.Pupil` 
    plane to an image plane using the Fraunhoffer diffraction approximation.

    Parameters
    ----------
    wavefront : :class:`~lentil.Wavefront`
        The wavefront to propagate
    pixelscale : float or (2,) array_like
        Physical sampling of output plane in meters. If `pixelscale` is a
        scalar, the sampling is assumed to be uniform in x and y.
    npix : int or (2,) array_like
        Output plane shape. If `npix` is a scalar, the output plane is assumed
        to be square with shape (npix,npix).
    npix_prop : int, (2,) array_like, or None, optional
        Propagation plane shape. If None (default), `npix_prop` = `npix`. If
        `npix_prop` is a scalar, the propagation plane is assumed to be square
        with shape (npix_prop,npix_prop).
    oversample : int, optional
        Number of times to oversample the propagation. Default is 2.
    place_tiles : bool, optional
        ?
    
    Returns
    -------
    wavefront : :class:`~lentil.Wavefront`
        A new wavefront propagated to the specified image plane

    """
    
    # TODO: this should return a NEW wavefront

    npix = _sanitize_shape(npix)
    npix_prop = _sanitize_shape(npix_prop, default=npix)

    out_shape = npix * oversample
    prop_shape = npix_prop * oversample

    center = wavefront.center(pixelscale, oversample)
    data = wavefront.data
    wavefront.data = [np.zeros(out_shape, dtype=complex)]

    # integer portion of the shift that will be accounted for later
    fix_center = np.fix(center)
    wavefront.shift = list(fix_center)

    # subpixel portion of the shift passed to the DFT
    subpixel_center = center - fix_center

    alpha = _dft_alpha(dx=wavefront.pixelscale, du=pixelscale,
                       wave=wavefront.wavelength, z=wavefront.focal_length,
                       oversample=oversample)

    if center.shape[0] != len(data):
        raise ValueError('dimension mismatch between center and wavefront depth')

    for d in range(len(data)):
        data[d] = fourier.dft2(data[d], alpha, prop_shape, subpixel_center[d],
                               offset=wavefront.offset[d], unitary=True)

    ### place tiles
    # TODO: collapse this into the above for loop
    tiles = []
    for d in range(len(data)):
        # The array center is given in terms of (x,y) but we place the chip in
        # terms of (r,c)
        center = np.flip(wavefront.shift[d], axis=0)

        # Compute the chip location
        data_slice, chip_slice = _chip_insertion_slices(out_shape, 
                                                        (data[d].shape[0], 
                                                         data[d].shape[1],),
                                                        center)

        if data_slice:
            tiles.append(imtile(data[d], data_slice, chip_slice))
        
        wavefront.tiles = consolidate(tiles)

        # TODO: this is where we would skip a step if place_tiles=False
        # Note also that we should move away from storing tiles in wf.tiles and 
        # store them in wf.data?
        for tile in wavefront.tiles:
            wavefront.data[0][tile.slice] = tile.data

    return wavefront


def _sanitize_shape(shape, default=()):
    if shape is None:
        shape = default
    shape = np.asarray(shape)
    if shape.shape == ():
        shape = np.append(shape, shape)
    return shape


def _dft_alpha(dx, du, wave, z, oversample):
    return (dx*du)/(wave*z*oversample)
    

def _standardize_shape(shape, default=()):
    if shape is None:
        shape = default
    shape = np.asarray(shape)
    if shape.shape == ():
        shape = np.append(shape, shape)
    return shape


def _standardize_bandpass(vec, default=()):
    if vec is None:
        vec = default
    vec = np.asarray(vec)
    if vec.shape == ():
        vec = vec[np.newaxis, ...]
    return vec


class _iterate_planes:
    def __init__(self, planes):
        self.planes = planes
        self.length = len(planes)
        self.n = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.n < self.length-1:
            plane = self.planes[self.n]
            next_plane = self.planes[self.n+1]
            self.n += 1
            return plane, next_plane
        else:
            raise StopIteration()


class _iterate_planes_reverse:
    def __init__(self, planes):
        self.planes = planes
        self.length = len(planes)
        self.n = 1

    def __iter__(self):
        return self

    def __next__(self):
        if self.n < self.length:
            plane = self.planes[-(self.n+1)]
            next_plane = self.planes[-self.n]
            self.n += 1
            return plane, next_plane
        else:
            raise StopIteration()


def _chip_insertion_slices(npix_canvas, npix_chip, shift):
    npix_canvas = np.asarray(npix_canvas)
    npix_chip = np.asarray(npix_chip)

    # Canvas coordinates of the upper left corner of the shifted chip
    chip_shifted_ul = (npix_canvas / 2) - (npix_chip / 2) + shift

    # Chip slice indices
    chip_top = int(0)
    chip_bottom = int(npix_chip[0])
    chip_left = int(0)
    chip_right = int(npix_chip[1])

    # Canvas insertion slice indices
    canvas_top = int(chip_shifted_ul[0])
    canvas_bottom = int(chip_shifted_ul[0] + npix_chip[0])
    canvas_left = int(chip_shifted_ul[1])
    canvas_right = int(chip_shifted_ul[1] + npix_chip[1])

    # reconcile the chip and canvas insertion indices
    if canvas_top < 0:
        chip_top = -1 * canvas_top
        canvas_top = 0

    if canvas_bottom > npix_canvas[0]:
        chip_bottom -= canvas_bottom - npix_canvas[0]
        canvas_bottom = npix_canvas[0]

    if canvas_left < 0:
        chip_left = -1 * canvas_left
        canvas_left = 0

    if canvas_right > npix_canvas[1]:
        chip_right -= canvas_right - npix_canvas[1]
        canvas_right = npix_canvas[1]

    if np.any(np.array([canvas_bottom, chip_bottom, canvas_right, chip_right]) < 0):
        return None, None
    else:
        return (slice(canvas_top, canvas_bottom), slice(canvas_left, canvas_right)), \
               (slice(chip_top, chip_bottom), slice(chip_left, chip_right))


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
                        dtype=np.complex128)

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
