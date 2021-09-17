import copy
from itertools import combinations

import numpy as np

import lentil


class Wavefront:
    """A class representing a monochromatic wavefront. :class:`Wavefront` is
    used internally by Lentil to perform diffraction propagation calculations.

    Parameters
    ----------
    wavelength : float
        Wavelength in meters

    shape : array_like
        Wavefront shape

    pixelscale : float, optional
        Wavefront array spatial sampling in meters/pixel

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

    __slots__ = ('wavelength', 'pixelscale', 'focal_length', 'offset',
                 'data', 'shift', 'tiles')

    def __init__(self, wavelength, shape=(), pixelscale=None):

        self.wavelength = wavelength
        self.pixelscale = pixelscale

        shape = np.asarray(shape)
        if shape.size < 3:
            self.data = [np.ones(shape, dtype=complex)]
        else:
            self.data = [np.ones((shape[1], shape[2]), dtype=complex) for d in shape[0]]

        # Wavefront focal length (which is infinity for a plane wave)
        self.focal_length = np.inf

        self.shift = []  # List of pre-propagation shifts
        self.offset = []  # List of (r,c) offsets from master array center for cropped DFTs

    def __mul__(self, plane):
        return plane.multiply(self)

    def __rmul__(self, other):
        return self.__mul__(other)

    @property
    def shape(self):
        """Size of data array"""
        return self.data[0].shape

    @property
    def depth(self):
        """Number of individual Wavefront arrays in :attr:`data`"""
        return len(self.data)

    @property
    def field(self):
        pass

    @property
    def intensity(self):
        # TODO: return np.abs(self.field**2)
        if self.tiles:
            out = np.zeros(self.shape, dtype=float)
            for tile in self.tiles:
                out[tile.slice] = np.abs(tile.data**2)
        else:
            out = np.abs(self.data**2)
        return out

    def copy(self):
        return copy.deepcopy(self)

    def center(self, pixelscale, oversample):
        """Compute the Wavefront center which may be shifted due to WavefrontShift
        objects in Wavefront.shift.

        This is a somewhat tricky method. Fundamentally it iterates over the
        :attr:`~lentil.Wavefront.tilt` list and computes the resulting shift in
        terms of number of pixels in oversampled space. This calculation is
        complicated by the fact that in some cases, an element in
        :attr:`~lentil.Wavefront.tilt` will itself be a list. In this case, the
        shift should be tracked individually for each entry in the list. All
        ensuing calculations should be done in parallel (i.e. the
        multi-dimensional shift array should not be recollapsed. This behavior
        allows SegmentedPupil to handle segment tilts individually.

        Parameters
        ----------
        pixelscale : float
            Image plane spatial sampling in meters/pixel

        oversample : int
            Oversampling factor

        Returns
        -------
        shift : (depth, 2) ndarray
            Image plane shift in number of (possibly oversampled) pixels

        """

        # Example:
        # tilt = [Shift(10,10), [Shift(100,100), Shift(200,200)], Shift(50,50)]
        # Beginning shift = [0,0]
        # After first tilt:
        #   shift = [10,10]
        # After second tilt, shift is duplicated before each shift is applied:
        #   shift = [[110,110], [210,210]]
        # All successive tilts are now applied in parallel:
        #   shift = [[160,160], [260,260]]

        center = np.zeros((1, 2))

        for shift in self.shift:
            if isinstance(shift, list):

                # Reshape the center array. It should have shape (len(tilt), 2).
                # If center.shape is (1,2), we'll duplicate center along the 0
                # axis so that it has shape (len(shift),2). If center.shape is
                # anything else, we can assume that the above duplication has
                # already happened so we'll just verify that the sizes have
                # remained consistent.
                if center.shape[0] == 1:
                    center = np.repeat(center, len(shift), axis=0)
                else:
                    assert center.shape[0] == len(shift)

                # Now we can iterate over the shifts
                for dim, shft in enumerate(shift):
                    center[dim, 0], center[dim, 1] = shft.shift(xs=center[dim, 0],
                                                                ys=center[dim, 1],
                                                                z=self.focal_length,
                                                                wavelength=self.wavelength)
            else:
                for dim in np.arange(center.shape[0]):
                    center[dim, 0], center[dim, 1] = shift.shift(xs=center[dim, 0],
                                                                 ys=center[dim, 1],
                                                                 z=self.focal_length,
                                                                 wavelength=self.wavelength)

        center = (center/pixelscale) * oversample

        if center.shape[0] != self.depth:
            center = np.repeat(center, self.depth, axis=0)

        return center

    def propagate_image(self, pixelscale, npix, npix_prop=None, oversample=2):
        """Propagate a :class:`~lentil.Wavefront` from a :class:`~lentil.Pupil`
        plane to an image plane using the Fraunhoffer diffraction approximation.

        Parameters
        ----------
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

        center = self.center(pixelscale, oversample)
        data = self.data
        self.data = [np.zeros(out_shape, dtype=complex)]

        # integer portion of the shift that will be accounted for later
        fix_center = np.fix(center)
        self.shift = list(fix_center)

        # subpixel portion of the shift passed to the DFT
        subpixel_center = center - fix_center

        alpha = _dft_alpha(dx=self.pixelscale, du=pixelscale,
                           wave=self.wavelength, z=self.focal_length,
                           oversample=oversample)

        if center.shape[0] != len(data):
            raise ValueError('dimension mismatch between center and wavefront depth')

        for d in range(len(data)):
            data[d] = lentil.fourier.dft2(data[d], alpha, prop_shape, subpixel_center[d],
                                          offset=self.offset[d], unitary=True)

        ### place tiles
        # TODO: collapse this into the above for loop
        tiles = []
        for d in range(len(data)):
            # The array center is given in terms of (x,y) but we place the chip in
            # terms of (r,c)
            center = np.flip(self.shift[d], axis=0)

            # Compute the chip location
            data_slice, chip_slice = _chip_insertion_slices(out_shape,
                                                            (data[d].shape[0],
                                                             data[d].shape[1],),
                                                            center)

            if data_slice:
                tiles.append(imtile(data[d], data_slice, chip_slice))

            self.tiles = consolidate(tiles)

            # TODO: this is where we would skip a step if place_tiles=False
            # Note also that we should move away from storing tiles in wf.tiles and
            # store them in wf.data?
            for tile in self.tiles:
                self.data[0][tile.slice] = tile.data

        return self


def _sanitize_shape(shape, default=()):
    if shape is None:
        shape = default
    shape = np.asarray(shape)
    if shape.shape == ():
        shape = np.append(shape, shape)
    return shape


def _dft_alpha(dx, du, wave, z, oversample):
    return (dx*du)/(wave*z*oversample)


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
