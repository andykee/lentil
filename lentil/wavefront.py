import copy

import numpy as np

import lentil
from lentil import Tilt
import lentil.field
from lentil.field import Field
import lentil.fourier
from lentil.fresnel import GaussianBeam
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
    focal_length : float or None, optional
        Wavefront focal length. A plane wave (default) has an infinite focal
        length (``None``).
    tilt: (2,) array_like, optional
        Radians of wavefront tilt about the x and y axes provided as
        ``[rx, ry]``. Default is ``[0, 0]`` (no tilt).
    ptype : lentil.ptype, optional
        Plane type. Default is ``lentil.none``.
    z : float, optional
        Axial position of the wavefront. Default is 0.
    pilot : :class:`~lentil.fresnel.GaussianBeam`, optional
        Pilot beam used to steer near-field propagation decisions. If None
        (default) and ``diameter`` is provided, a pilot beam is
        automatically created with ``waist_radius = diameter/2`` and its
        waist located at ``z``.

    """
    def __init__(self, wavelength, pixelscale=None, diameter=None, focal_length=None,
                 tilt=None, ptype=None, z=0, pilot=None):

        #: float: Axial position of the wavefront
        self.z = z

        #: str: Reference surface the field data is defined against
        #: ('planar' or 'spherical')
        self.reference = 'planar'

        #: float: Accumulated optical path in meters (path ledger)
        self.path = 0.0

        # Curvature state is stored as the axial location of the analytic
        # focus. focal_length is derived from it and z.
        self._z_focus = z + focal_length if focal_length else None

        #: float: Wavefront diameter
        self.diameter = diameter

        if pilot is not None:
            #: GaussianBeam or None: Pilot beam
            self.pilot = pilot
        elif diameter is not None:
            self.pilot = GaussianBeam(diameter/2, wavelength, waist_position=z)
        else:
            self.pilot = None

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
        return plane.__mul__(self)

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
    def z_focus(self):
        """Axial location of the analytic focus

        A plane wave has ``z_focus = None``. Fixed between lenses; only
        multiplication by a curvature-bearing plane changes it.

        Returns
        -------
        float or None
        """
        return self._z_focus

    @property
    def focal_length(self):
        """Distance from the wavefront's current position to the analytic
        focus

        Derived from :attr:`z_focus` and :attr:`z`. A plane wave (default)
        has an infinite focal length (``None``).

        Returns
        -------
        float or None
        """
        if self._z_focus is None:
            return None
        return self._z_focus - self.z

    @focal_length.setter
    def focal_length(self, value):
        if value is None or np.isinf(value):
            self._z_focus = None
        else:
            self._z_focus = self.z + value

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
              tilt=None, shape=None, ptype=None, z=0, pilot=None,
              reference='planar', path=0.0):
        """Create an empty Wavefront

        The resulting wavefront will have an empty :attr:`data` attribute.

        Parameters
        ----------

        """
        w = cls(wavelength=wavelength, pixelscale=pixelscale, diameter=diameter,
                focal_length=focal_length, tilt=tilt, ptype=ptype, z=z,
                pilot=pilot)
        w.reference = reference
        w.path = path
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
