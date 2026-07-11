"""Near-field (Fresnel) propagation support.

This module implements the near-field propagation machinery described in
[1], including the pilot (surrogate) Gaussian beam used to select reference
surfaces and propagator composition.

Sign conventions
----------------
Lentil adopts the following conventions for beams traveling in the +z
direction:

* The phase radius of curvature ``R(z) = (z - z_w0) + z_R**2/(z - z_w0)``
  is positive for a diverging beam (waist behind, at smaller z), negative
  for a converging beam (waist ahead, at larger z), and infinite at the
  waist.
* A thin lens with focal length ``f`` (positive = converging) transforms
  the phase radius according to ``1/R2 = 1/R1 - 1/f``.
* A wavefront's ``focal_length`` (distance to the analytic focus) relates
  to the local phase radius by ``focal_length = -R``.

References
----------
[1] Lawrence, G. N. Optical Modeling, in Applied Optics and Optical
    Engineering v11, ch. 3 (1992)

"""

import numpy as np


class GaussianBeam:
    """Gaussian beam parameters.

    Primarily used as the *pilot beam* carried by a
    :class:`~lentil.Wavefront` (as ``wavefront.pilot``) to steer near-field
    propagation decisions: reference surface selection and propagator
    composition. In that role the pilot is a bookkeeping surrogate for the
    actual beam, not physics — it only needs to be representative enough to
    pick sampling regimes sensibly [1].

    Parameters
    ----------
    waist_radius : float
        Beam waist radius (1/e field amplitude radius) in meters
    wavelength : float
        Wavelength in meters
    waist_position : float, optional
        Axial location of the beam waist in meters. Default is 0.

    References
    ----------
    [1] Lawrence, G. N. Optical Modeling, in Applied Optics and Optical
        Engineering v11, ch. 3 (1992)

    """
    __slots__ = ('waist_radius', 'wavelength', 'waist_position')

    def __init__(self, waist_radius, wavelength, waist_position=0):
        if waist_radius <= 0:
            raise ValueError('waist_radius must be positive')
        if wavelength <= 0:
            raise ValueError('wavelength must be positive')

        #: float: Beam waist radius
        self.waist_radius = waist_radius

        #: float: Wavelength
        self.wavelength = wavelength

        #: float: Axial location of the beam waist
        self.waist_position = waist_position

    def __repr__(self):
        return (f'{self.__class__.__name__}(waist_radius='
                f'{self.waist_radius:.4g}, wavelength='
                f'{self.wavelength:.4g}, waist_position='
                f'{self.waist_position:.4g})')

    @property
    def rayleigh_distance(self):
        """Rayleigh distance :math:`z_R = \\pi w_0^2/\\lambda`

        Returns
        -------
        float
        """
        return np.pi * self.waist_radius**2 / self.wavelength

    @property
    def divergence(self):
        """Far-field divergence half-angle :math:`\\theta = w_0/z_R`
        in radians

        Returns
        -------
        float
        """
        return self.waist_radius / self.rayleigh_distance

    def radius(self, z):
        """Beam radius :math:`w(z) = w_0\\sqrt{1 + ((z-z_{w0})/z_R)^2}`

        Parameters
        ----------
        z : float
            Axial position

        Returns
        -------
        float
        """
        return self.waist_radius * np.sqrt(1 + ((z - self.waist_position)/self.rayleigh_distance)**2)

    def phase_radius(self, z):
        """Phase radius of curvature
        :math:`R(z) = (z-z_{w0}) + z_R^2/(z-z_{w0})`

        ``R > 0`` for a diverging beam (waist behind), ``R < 0`` for a
        converging beam (waist ahead). Returns ``inf`` at the waist.

        Parameters
        ----------
        z : float
            Axial position

        Returns
        -------
        float
        """
        dz = z - self.waist_position
        if dz == 0:
            return np.inf
        return dz + self.rayleigh_distance**2/dz

    def gouy(self, z):
        """Gouy phase shift :math:`\\theta(z) = \\arctan((z-z_{w0})/z_R)`
        in radians

        Parameters
        ----------
        z : float
            Axial position

        Returns
        -------
        float
        """
        return np.arctan((z - self.waist_position)/self.rayleigh_distance)

    def inside(self, z):
        """True if ``z`` is inside the Rayleigh distance of the waist

        A position inside the Rayleigh distance uses a planar reference
        surface; a position outside uses a spherical reference surface.

        Parameters
        ----------
        z : float
            Axial position

        Returns
        -------
        bool
        """
        return np.abs(z - self.waist_position) <= self.rayleigh_distance

    @classmethod
    def from_radius(cls, radius, phase_radius, wavelength, z=0):
        """Create a GaussianBeam from beam radius and phase radius known
        at some axial position (waist finding).

        Parameters
        ----------
        radius : float
            Beam radius :math:`w` at ``z``
        phase_radius : float
            Phase radius of curvature :math:`R` at ``z``. May be ``inf``
            (the beam has its waist at ``z``).
        wavelength : float
            Wavelength
        z : float, optional
            Axial position where ``radius`` and ``phase_radius`` are
            defined. Default is 0.

        Returns
        -------
        :class:`GaussianBeam`
        """
        if np.isinf(phase_radius):
            return cls(radius, wavelength, waist_position=z)

        # Lawrence Eq. 56, rewritten in Lentil's sign convention:
        #   z - z_w0 = R/(1 + (lambda*R/(pi*w**2))**2)
        #   w_0 = w/sqrt(1 + (pi*w**2/(lambda*R))**2)
        a = wavelength*phase_radius/(np.pi*radius**2)
        dz = phase_radius/(1 + a**2)
        waist_radius = radius/np.sqrt(1 + 1/a**2)
        return cls(waist_radius, wavelength, waist_position=z-dz)

    def lens(self, focal_length, z):
        """Return a new GaussianBeam transformed by a thin lens at ``z``.

        The beam radius is unchanged by the lens; the phase radius is
        transformed according to ``1/R2 = 1/R1 - 1/f`` (Lawrence Eq. 57 in
        Lentil's sign convention) and the new waist is found from the
        transformed parameters.

        Parameters
        ----------
        focal_length : float
            Lens focal length. Positive is converging. May be ``inf``
            (no optical power).
        z : float
            Axial position of the lens

        Returns
        -------
        :class:`GaussianBeam`
        """
        if np.isinf(focal_length):
            return GaussianBeam(self.waist_radius, self.wavelength,
                                waist_position=self.waist_position)

        w = self.radius(z)
        r1 = self.phase_radius(z)

        if np.isinf(r1):
            r2 = -focal_length
        elif r1 == focal_length:
            # collimating lens - output waist is at the lens
            return GaussianBeam(w, self.wavelength, waist_position=z)
        else:
            r2 = r1*focal_length/(focal_length - r1)

        return GaussianBeam.from_radius(w, r2, self.wavelength, z=z)
