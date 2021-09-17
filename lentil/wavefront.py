import copy
from itertools import combinations

import numpy as np



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
