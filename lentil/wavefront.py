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
    tilt: (2,) array_like, optional
        Radians of wavefront tilt about the x and y axes provided as 
        ``[rx, ry]``. Default is ``[0, 0]`` (no tilt).
    ptype : lentil.ptype, optional
        Plane type. Default is ``lentil.none``.
    
    """
    def __init__(self, wavelength, pixelscale=None, diameter=None, focal_length=None,
                 tilt=None, ptype=None):
        
        #: float: Wavefront focal length
        self.focal_length = focal_length if focal_length else np.inf

        #: float: Wavefront diameter
        self.diameter = diameter

        #: tuple of ints: Wavefront shape
        self.shape = ()

        self._wavelength = wavelength
        self._pixelscale = None if pixelscale is None else np.broadcast_to(pixelscale, (2,))
        self.ptype = lentil.ptype(ptype)

        if tilt is not None:
            if len(tilt) != 2:
                raise ValueError('tilt must be specified as [rx, ry]')
            tilt = [Tilt(x=tilt[0], y=tilt[1])]

        self.data = [Field(data=np.array(1, dtype=complex),
                           offset=None,
                           tilt=tilt)]

    def __mul__(self, plane):
        return plane.multiply(self)

    def __rmul__(self, other):
        return self.__mul__(other)

    @property
    def wavelength(self):
        """Wavefront wavelength
        
        Returns
        -------
        float
        """
        return self._wavelength

    @property
    def pixelscale(self):
        """Physical sampling of wavefront
        
        Returns
        -------
        tuple of floats
        """
        return self._pixelscale

    @property
    def ptype(self):
        """Wavefront plane type
        
        Returns
        -------
        ptype object
        """
        return self._ptype
    
    @ptype.setter
    def ptype(self, value):
        if lentil.ptype(value) not in (lentil.none, lentil.pupil, lentil.image):
            raise TypeError(f"invalid ptype '{value}' for Wavefront")
        else:
            self._ptype = lentil.ptype(value)

    @property
    def field(self):
        """Wavefront complex field
        
        Returns
        -------
        ndarray
        """
        out = np.zeros(self.shape, dtype=complex)
        for field in self.data:
            out = lentil.field.insert(field, out)
        return out

    @property
    def intensity(self):
        """Wavefront intensity
        
        Returns
        -------
        ndarray
        """
        out = np.zeros(self.shape, dtype=float)
        for field in lentil.field.reduce(self.data):
            out = lentil.field.insert(field, out, intensity=True)
        return out

    @classmethod
    def empty(cls, wavelength, pixelscale=None, diameter=None, focal_length=None,
              tilt=None, shape=None, ptype=None):
        """Create an empty Wavefront

        The resulting wavefront will have an empty :attr:`data` attribute.

        Parameters
        ----------
        
        """
        w = cls(wavelength=wavelength, pixelscale=pixelscale, diameter=diameter,
                focal_length=focal_length, tilt=tilt, ptype=ptype)
        w.data = []
        w.shape = () if shape is None else shape
        return w
    

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
        for field in lentil.field.reduce(self.data):
            out = lentil.field.insert(field, out, intensity=True, weight=weight)
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
