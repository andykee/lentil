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

    __slots__ = ('wavelength', 'pixelscale', 'focal_length', 'tilt', 'offset',
                 'data', 'center', '_planetype')

    def __init__(self, wavelength, shape=None, pixelscale=None, planetype=None):

        self.wavelength = wavelength
        self.pixelscale = pixelscale
        self.planetype = planetype

        # All new Wavefront objects represent a perfect plane wave
        # TODO: default shape=1, make sure we create separate list entries of shape
        # has ndim = 3
        if shape:
            self.data = [np.ones((shape[0], shape[1]), dtype=np.complex)]
        else:
            self.data = [np.array([1], dtype=np.complex)]

        # Wavefront focal length (which is infinity for a plane wave)
        self.focal_length = np.inf

        self.tilt = []  # List of pre-propagation tilts
        self.offset = []  # List of (r,c) offsets from master array center for cropped DFTs
        self.center = []  # List of (r,c) center shifts of each item in data

    @property
    def depth(self):
        """Number of individual Wavefront arrays in :attr:`data`"""
        return len(self.data)

    @property
    def planetype(self):
        return self._planetype

    @planetype.setter
    def planetype(self, value):
        assert value in {'pupil', 'image', None}
        self._planetype = value

    def shift(self, pixelscale, oversample):
        """Compute image plane shift due to wavefront tilt.

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

        shift = np.zeros((1, 2))

        for tilt in self.tilt:
            if isinstance(tilt, list):

                # Reshape the shift array. It should have shape (len(tilt), 2).
                # If shift.shape is (1,2), we'll duplicate shift along the 0
                # axis so that it has shape (len(tilt),2). If shift.shape is
                # anything else, we can assume that the above duplication has
                # already happened so we'll just verify that the sizes have
                # remained consistent.
                if shift.shape[0] == 1:
                    shift = np.repeat(shift, len(tilt), axis=0)
                else:
                    assert shift.shape[0] == len(tilt)

                # Now we can iterate over the tilts
                for d, t in enumerate(tilt):
                    shift[d, 0], shift[d, 1] = t.shift(xs=shift[d, 0],
                                                       ys=shift[d, 1],
                                                       z=self.focal_length,
                                                       wavelength=self.wavelength)
            else:
                for d in np.arange(shift.shape[0]):
                    shift[d, 0], shift[d, 1] = tilt.shift(xs=shift[d, 0],
                                                          ys=shift[d, 1],
                                                          z=self.focal_length,
                                                          wavelength=self.wavelength)

        shift = (shift/pixelscale) * oversample

        if shift.shape[0] != self.depth:
            shift = np.repeat(shift, self.depth, axis=0)

        return shift
