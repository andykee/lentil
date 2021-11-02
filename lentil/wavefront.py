import copy

import numpy as np

import lentil.field
from lentil.field import Field
import lentil.fourier
import lentil.helper


class Wavefront:
    """A class representing a monochromatic wavefront.

    Parameters
    ----------
    wavelength : float
        Wavelength in meters
    pixelscale : float, optional
        Wavefront array spatial sampling in meters/pixel
    shape : (2,) array_like, optional
        Wavefront shape. If ``shape`` is None (default), the wavefront is
        assumed to be infinite (broadcastable to any shape).
    data : list_like, optional
        Wavefront data. Default is [1+0j] (a plane wave).
    focal_length : float or np.inf
        Wavefront focal length. A plane wave (default) has an infinite focal
        length (``np.inf``).

    """
    __slots__ = ('wavelength', '_pixelscale', 'focal_length',
                 'data', 'shape', 'planetype')

    def __init__(self, wavelength, pixelscale=None, shape=None, planetype=None,
                 data=None, focal_length=None):

        self.wavelength = wavelength
        self.pixelscale = pixelscale
        self.shape = shape

        if data is None:
            self.data = [Field(data=1, pixelscale=pixelscale, offset=[0, 0], tilt=[])]
        else:
            self.data = [*data]

        self.focal_length = focal_length if focal_length else np.inf
        self.planetype = planetype

    def __mul__(self, plane):
        return plane.multiply(self)

    def __rmul__(self, other):
        return self.__mul__(other)

    @property
    def pixelscale(self):
        return self._pixelscale

    @pixelscale.setter
    def pixelscale(self, value):
        if value is not None:
            self._pixelscale = lentil.helper.sanitize_shape(value)
        else:
            self._pixelscale = None

    @property
    def field(self):
        shape = self.shape if self.shape is not None else ()
        out = np.zeros(shape, dtype=complex)
        for field in self.data:
            out = lentil.field.insert(field, out)
        return out

    @property
    def intensity(self):
        shape = self.shape if self.shape is not None else ()
        out = np.zeros(shape, dtype=float)
        for field in lentil.field.reduce(*self.data):
            out = lentil.field.insert(field, out, intensity=True)
        return out

    def copy(self):
        return copy.deepcopy(self)

    def propagate_image(self, pixelscale, npix, npix_prop=None, oversample=2):
        """Propagate the Wavefront from a Pupil to an Image plane using
        Fraunhofer diffraction.

        Parameters
        ----------
        pixelscale : float or (2,) float
            Physical sampling of output (image) plane. If a single value is supplied,
            the output is assumed to be uniformly sampled in both x and y.
        npix : int or (2,) tuple of ints
            Shape of output plane.
        npix_prop : int or (2,) tuple of ints, optional
            Shape of propagation output plane. If None (default),
            ``npix_prop = npix``. If ``npix_prop != npix``, the propagation
            result is placed in the appropriate location in the output plane.
            npix_prop cannot be larger than npix.
        oversample : int, optional
            Number of times to oversample the output plane. Default is 2.

        Returns
        -------
        wavefront : :class:`~lentil.Wavefront`
            A new Wavefront propagated to the specified image plane

        """
        if self.planetype != 'pupil':
            raise ValueError("Wavefront must have planetype 'pupil'")

        npix = lentil.helper.sanitize_shape(npix)
        npix_prop = lentil.helper.sanitize_shape(npix_prop, default=npix)
        prop_shape = npix_prop * oversample

        out = Wavefront(wavelength=self.wavelength, data=[],
                        pixelscale=pixelscale*oversample, shape=npix*oversample,
                        planetype='image')

        for field in self.data:
            # compute the field shift from any embedded tilts. note the return value
            # is specified in terms of (r, c)
            shift = field.shift(z=self.focal_length, wavelength=self.wavelength,
                                pixelscale=pixelscale, oversample=oversample,
                                indexing='ij')

            fix_shift = np.fix(shift)
            dft_shift = shift - fix_shift

            if _overlap(prop_shape, fix_shift, out.shape):
                alpha = lentil.helper.dft_alpha(dx=self.pixelscale, du=pixelscale,
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
