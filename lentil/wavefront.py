import copy

import numpy as np

import lentil
from lentil import Tilt
import lentil.field
from lentil.field import Field
import lentil.fourier
from lentil.fresnel import GaussianBeam, _chirp
import lentil.helper

# sentinel for Wavefront.derive() deltas. None can't serve: it is a
# meaningful value for several kinds of Wavefront state (z_focus=None is
# a plane wave, pilot=None is no pilot, ...)
_UNSET = object()

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
        # focus. focal_length is derived from it and z. Note focal_length=0
        # is meaningful (the wavefront is at its focus).
        if focal_length is None or np.isinf(focal_length):
            self._z_focus = None
        else:
            self._z_focus = z + focal_length

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

    def phasor(self):
        """Materialize the wavefront's total complex field as a sampled
        array.

        The returned array is the stored field data with all analytically
        bookkept phase multiplied back in:

        * per-field tilt (:class:`~lentil.Tilt` objects)
        * the reference-surface quadratic phase (for a
          spherical-reference wavefront, the sphere centered on
          :attr:`z_focus`)
        * the accumulated path ledger piston ``exp(j*k*path)``

        .. warning::

            Lentil bookkeeps these terms analytically precisely because
            they may not be representable at the array's sampling. The
            materialized phasor can alias badly — a fast beam's
            reference sphere or a large tilt can wrap many times per
            pixel. Use for diagnostics, external handoff, and
            small-phase regimes; never as an input to Lentil's own
            propagators.

        Returns
        -------
        ndarray

        See Also
        --------
        Wavefront.field : the raw (reference-relative) field data
        """
        k = 2*np.pi/self.wavelength

        out = np.zeros(self.shape, dtype=complex)
        for field in self.data:
            if field.tilt:
                data = field.data * self._tilt_phasor(field)
                field = Field(data=data, pixelscale=field.pixelscale,
                              offset=field.offset)
            out = lentil.field.insert(field, out)

        if self.reference == 'spherical':
            if self.pixelscale is None:
                raise ValueError("can't materialize spherical reference "
                                 "phase with pixelscale = None")
            # the array is defined relative to a sphere centered on
            # z_focus; the absorbed phase is q(dz_ref) with
            # dz_ref = z - z_focus = -focal_length
            out = out * _chirp(out.shape, (0, 0), self.pixelscale,
                               self.wavelength, self.z - self._z_focus)

        return out * np.exp(1j*k*self.path)

    def _tilt_phasor(self, field):
        # reconstruct the sampled phasor of a field's analytic Tilt
        # objects, evaluated at the field's own (offset) coordinates.
        # This is the exact inverse of Plane.fit_tilt's OPD removal:
        # opd_tilt = t1*(r*dr) + t2*(-c*dc) with Tilt(x=t1, y=t2)
        # stored as T.x = t2, T.y = t1
        if self.pixelscale is None:
            raise ValueError("can't materialize tilt with pixelscale = None")
        for tilt in field.tilt:
            if type(tilt) is not Tilt:
                raise NotImplementedError(
                    f"can't materialize phase for tilt object of type "
                    f"'{type(tilt).__name__}'")

        n0, n1 = field.shape
        r = (np.arange(n0) - n0//2 + field.offset[0])*self.pixelscale[0]
        c = (np.arange(n1) - n1//2 + field.offset[1])*self.pixelscale[1]
        opd = np.zeros((n0, n1))
        for tilt in field.tilt:
            opd = opd + tilt.y*r[:, np.newaxis] - tilt.x*c[np.newaxis, :]
        return np.exp(2j*np.pi*opd/self.wavelength)

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
    

    def derive(self, *, dz=0, dpath=None, wavelength=_UNSET, pixelscale=_UNSET,
               shape=_UNSET, ptype=_UNSET, z_focus=_UNSET, pilot=_UNSET,
               reference=_UNSET):
        """Create a new Wavefront derived from this one.

        A derived wavefront inherits **all** Wavefront-level bookkeeping
        state (wavelength, sampling, diameter, position, analytic
        curvature, pilot beam, reference surface, path ledger, ptype)
        unless a delta explicitly states otherwise, and always has an
        empty :attr:`data` attribute — Field-level state (data, offset,
        tilt) is the caller's responsibility. This is the required way
        for propagators and planes to construct their output wavefronts:
        state that is not explicitly changed cannot be accidentally
        dropped.

        The inheritance is structural (the entire instance state is
        copied before deltas are applied), so attributes added to
        Wavefront in the future are inherited automatically.

        Parameters
        ----------
        dz : float, optional
            Signed axial distance to advance the wavefront. Advances
            both :attr:`z` and the :attr:`path` ledger. Default is 0.
        dpath : float, optional
            Amount added to the path ledger, overriding the default
            advance of ``dz``. Used for ledger corrections beyond the
            geometric distance (e.g. the dropped Fresnel ``1/j`` piston
            of a waist transform). There is deliberately no absolute
            path setter.
        wavelength, pixelscale, shape, ptype, z_focus, pilot, reference :
            Optional replacement values. Note ``z_focus`` (the analytic
            curvature state) transfers directly; there is deliberately
            no ``focal_length`` delta — the ``z + focal_length``
            round-trip is error-prone and ``focal_length`` remains a
            derived property.

        Returns
        -------
        :class:`~lentil.Wavefront`

        Notes
        -----
        The copy is shallow: in particular, :attr:`pilot` is shared with
        the source wavefront. This is safe because ``GaussianBeam``
        objects are treated as immutable throughout Lentil (lenses
        replace the pilot rather than mutating it).
        """
        out = self.__class__.__new__(self.__class__)
        out.__dict__.update(self.__dict__)
        out.data = []

        out.z = self.z + dz
        out.path = self.path + (dz if dpath is None else dpath)

        if wavelength is not _UNSET:
            out._wavelength = wavelength
        if pixelscale is not _UNSET:
            out._pixelscale = None if pixelscale is None else np.broadcast_to(pixelscale, (2,))
        if shape is not _UNSET:
            out.shape = () if shape is None else shape
        if ptype is not _UNSET:
            out.ptype = ptype
        if z_focus is not _UNSET:
            out._z_focus = z_focus
        if pilot is not _UNSET:
            out.pilot = pilot
        if reference is not _UNSET:
            out.reference = reference

        # state invariant: a spherical reference surface is a sphere
        # centered on the analytic focus, which must therefore exist
        if out.reference == 'spherical' and out._z_focus is None:
            raise ValueError("spherical-reference wavefront requires a "
                             "defined z_focus")

        return out

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
