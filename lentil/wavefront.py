import copy

import numpy as np

import lentil
from lentil import Tilt
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
        Physical sampling of wavefront
    diameter: float, optional
        Wavefront diameter. Default is None
    focal_length : float or np.inf, optional
        Wavefront focal length. A plane wave (default) has an infinite focal
        length (``np.inf``).
    tilt: list_like, optional

    shape : (2,) array_like
        Wavefront shape. If ``shape`` is None (default), the wavefront is
        assumed to be infinite (broadcastable to any shape).
    ptype : lentil.ptype
        Plane type.

    Attributes
    ----------
    data : list_like
        Wavefront data. Default is [1+0j] (a plane wave).
    
    """
    __slots__ = ('wavelength', 'pixelscale', 'focal_length', 'diameter',
                 'focal_length', 'ptype', 'shape', 'data')

    def __init__(self, wavelength, pixelscale=None, diameter=None, focal_length=None,
                 tilt=None, ptype=None):

        self.wavelength = wavelength
        self.pixelscale = None if pixelscale is None else np.broadcast_to(pixelscale, (2,))
        self.focal_length = focal_length if focal_length else np.inf
        self.diameter = diameter
        self.ptype = lentil.ptype(ptype)
        self.shape = ()

        if tilt is not None:
            if len(tilt) != 2:
                raise ValueError('tilt must be specified as [rx, ry]')
            tilt = [Tilt(x=tilt[0], y=tilt[1])]

        self.data = [Field(data=np.array(1, dtype=complex),
                           offset=None,
                           tilt=tilt)]

    def __mul__(self, plane):
        return plane.multiply(self, inplace=False)

    def __imul__(self, plane):
        return plane.multiply(self, inplace=True)

    def __rmul__(self, other):
        return self.__mul__(other)

    @property
    def field(self):
        """Wavefront complex field"""
        out = np.zeros(self.shape, dtype=complex)
        for field in self.data:
            out = lentil.field.insert(field, out)
        return out

    @property
    def intensity(self):
        """Wavefront intensity"""
        out = np.zeros(self.shape, dtype=float)
        for field in lentil.field.reduce(self.data):
            out = lentil.field.insert(field, out, intensity=True)
        return out

    @classmethod
    def empty(cls, wavelength, pixelscale=None, diameter=None, focal_length=None,
              tilt=None, shape=None, ptype=None):
        w = cls(wavelength=wavelength, pixelscale=pixelscale, diameter=diameter,
                focal_length=focal_length, tilt=tilt, ptype=ptype)
        w.data = []
        w.shape = () if shape is None else shape
        return w
        
    def copy(self):
        return copy.deepcopy(self)

    def insert(self, out, weight=1):
        """Directly insert wavefront intensity data into an output array.

        This method can avoid repeatedly allocating large arrays of zeros
        when accumulating :attr:`intensity`.

        Parameters
        ----------
        out : ndarray
            Array to insert wavefront data into
        weight : float
            Scale factor applied to wavefront data

        Returns
        -------
        out : ndarray
            Array with wavefront data inserted into it at the appropriate location
        """
        for field in lentil.field.reduce(*self.data):
            out = lentil.field.insert(field, out, intensity=True, weight=weight)
        return out

    def propagate_image(self, pixelscale, npix, npix_prop=None, oversample=2,
                        inplace=True):
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
        inplace : bool, optional
            If True (default) the wavefront is propagated in-place, otherwise
            a copy is created and propagated.

        Returns
        -------
        wavefront : :class:`~lentil.Wavefront`
            A Wavefront propagated to the specified image plane

        """
        if self.ptype != lentil.pupil:
            raise ValueError("Wavefront must have planetype 'pupil'")

        npix = np.asarray(lentil.sanitize_shape(npix))
        npix_prop = npix if npix_prop is None else np.asarray(lentil.sanitize_shape(npix_prop))
        prop_shape = npix_prop * oversample

        dx = self.pixelscale
        du = np.asarray(lentil.sanitize_shape(pixelscale))
        z = self.focal_length
        data = self.data

        if inplace:
            out = self
            out.data = []
            out.pixelscale = du / oversample
            out.shape = npix * oversample
            out.focal_length = np.inf
            out.ptype = lentil.image
        else:
            out = Wavefront.empty(wavelength=self.wavelength,
                                  pixelscale=du/oversample,
                                  shape=npix*oversample,
                                  ptype=lentil.image)

        for field in data:
            # compute the field shift from any embedded tilts. note the return value
            # is specified in terms of (r, c)
            shift = field.shift(z=z, wavelength=self.wavelength,
                                pixelscale=du, oversample=oversample,
                                indexing='ij')

            fix_shift = np.fix(shift)
            dft_shift = shift - fix_shift

            if _overlap(prop_shape, fix_shift, out.shape):
                alpha = lentil.helper.dft_alpha(dx=dx, du=du,
                                                wave=self.wavelength, z=z,
                                                oversample=oversample)
                data = lentil.fourier.dft2(f=field.data, alpha=alpha,
                                           npix=prop_shape, shift=dft_shift,
                                           offset=field.offset, unitary=True)
                out.data.append(Field(data=data, pixelscale=du/oversample,
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
