from copy import deepcopy

import numpy as np
from scipy import ndimage
import scipy.integrate
import scipy.optimize

from lentil import util

__all__ = ['Plane', 'Pupil', 'Image', 'DispersiveShift', 'Grism',
           'LensletArray', 'Tilt', 'Rotate', 'Flip']


class Plane:
    """Base class for representing a finite geometric plane.

    Parameters
    ----------
    pixelscale : float or (2,) array_like
        Physical sampling (in meters) of each pixel in the plane. If :attr:`pixelscale`
        is a scalar, uniform sampling in x and y is assumed.

    amplitude : array_like, optional
        Electric field amplitude transmission. Amplitude should be normalized
        with :func:`~lentil.util.normalize_power` if conservation of power
        through a diffraction propagation is required. If not specified, a
        default amplitude is created which has no effect on wavefront
        propagation.

    phase : array_like, optional
        Phase change caused by plane. If not specified, a default phase is
        created which has no effect on wavefront propagation.

    mask : array_like, optional
        Binary mask. If not specified, a mask is created from the amplitude.

    Attributes
    ----------
    tilt : list

    """

    def __init__(self, pixelscale=None, amplitude=1, phase=0, mask=None):
        # We directly set the local attributes here in case a subclass has redefined
        # the property (which could cause an weird behavior and will throw an
        # AttributeError if the subclass hasn't defined an accompanying getter
        self.pixelscale = pixelscale
        self._amplitude = np.asarray(amplitude)
        self._phase = np.asarray(phase)
        self._mask = np.asarray(mask) if mask is not None else None

        self.tilt = []

    def __init_subclass__(cls):
        # we have to define default values to avoid AttributeErrors in case
        # subclasses don't call super().__init__ and don't explicitly define
        # an attribute
        cls._pixelscale = None
        cls._amplitude = np.array(1)
        cls._phase = np.array(0)
        cls._mask = None

        cls.tilt = []

    def __repr__(self):
        return f'{self.__class__.__name__}()'

    @property
    def pixelscale(self):
        """Physical sampling (in meters) of each pixel in the plane.

        Returns
        -------
        pixelscale : ndarray
            (dx,dy) sampling of the Plane
        """
        return self._pixelscale

    @pixelscale.setter
    def pixelscale(self, value):
        if value is not None:
            value = np.asarray(value)
            if value.shape == ():
                value = np.append(value, value)
        self._pixelscale = value

    @property
    def amplitude(self):
        """Electric field amplitude transmission.

        Returns
        -------
        amplitude : ndarray

        """
        return self._amplitude

    @amplitude.setter
    def amplitude(self, value):
        self._amplitude = np.asarray(value)

    @property
    def phase(self):
        """Electric field phase shift.

        Returns
        -------
        phase : ndarray

        """
        return self._phase

    @phase.setter
    def phase(self, value):
        self._phase = np.asarray(value)

    @property
    def mask(self):
        """Binary mask

        Returns
        -------
        mask : ndarray

        """
        if self._mask is not None:
            return self._mask
        else:
            mask = np.copy(self.amplitude)
            mask[mask != 0] = 1
            return mask

    @mask.setter
    def mask(self, value):
        if value is not None:
            self._mask = np.asarray(value)
        else:
            self._mask = None

    @property
    def global_mask(self):
        if self.depth < 2:
            return self.mask
        else:
            return np.sum(self.mask, axis=0)

    @property
    def shape(self):
        """Plane dimensions computed from :attr:`mask` Returns None if :attr:`mask` 
        is None.
        """
        if self.depth == 1:
            return self.mask.shape
        else:
            return self.mask.shape[1], self.mask.shape[2]

    @property
    def depth(self):
        """Number of independent masks (segments) in self.mask"""
        if self.mask.ndim in (0, 1, 2):
            return 1
        else:
            return self.mask.shape[0]

    @staticmethod
    def rescale_amplitude(amplitude, scale):
        """Rescale plane amplitude

        This method uses 3rd order interpolation and preserves the power of the amplitude
        array. If custom interpolation is required, this method should be redefined.

        Parameters
        ----------
        amplitude : array_like

        scale : float
            Scaling factor. Scale factors less than 1 will shrink the amplitude. Scale
            factors greater than 1 will grow the amplitude.

        Returns
        -------
        amplitude : ndarray
            Rescaled amplitude

        """
        if scale == 1:
            return amplitude
        else:
            return util.rescale(amplitude, scale=scale, shape=None, mask=None,
                                order=3, mode='nearest', unitary=False)/scale**2

    @staticmethod
    def rescale_phase(phase, scale):
        """Rescale plane phase

        This method uses 3rd order interpolation. If custom interpolation is required,
        this method should be redefined.

        Parameters
        ----------
        phase : ndarray

        scale : float
            Scaling factor. Scale factors less than 1 will shrink the amplitude. Scale
            factors greater than 1 will grow the amplitude.

        Returns
        -------
        phase : ndarray
            Rescaled phase

        """
        if scale == 1:
            return phase
        else:
            return util.rescale(phase, scale=scale, shape=None, mask=None,
                                order=3, mode='nearest', unitary=False)

    @staticmethod
    def rescale_mask(mask, scale):
        """Rescale plane mask

        This method uses linear interpolation. If custom interpolation is required,
        this method should be redefined.

        Parameters
        ----------
        mask : ndarray

        scale : float
            Scaling factor. Scale factors less than 1 will shrink the amplitude. Scale
            factors greater than 1 will grow the amplitude.

        Returns
        -------
        mask : ndarray
            Rescaled mask

        """
        mask = util.rescale(mask, scale=scale, shape=None, mask=None, order=0,
                            mode='constant', unitary=False)
        # mask[mask < np.finfo(mask.dtype).eps] = 0
        mask[np.nonzero(mask)] = 1
        return mask.astype(np.int)

    @property
    def ptt_vector(self):
        """2D vector representing piston and tilt in x and y. Planes with no mask
        have ptt_vector = None.

        Returns
        -------
        ptt_vector : ndarray or None

        """

        # if there's no mask, we just set ptt_vector to None and move on
        if self.shape == () or self.shape is None:
            ptt_vector = None
        else:
            # compute unmasked piston, tip, tilt vector
            x, y = util.mesh(self.shape)
            unmasked_ptt_vector = np.einsum('ij,i->ij', [np.ones(x.size), x.ravel(), y.ravel()],
                                            [1, self.pixelscale[0], self.pixelscale[1]])

            if self.depth == 1:
                ptt_vector = np.einsum('ij,j->ij', unmasked_ptt_vector, self.mask.ravel())
            else:
                # prepare empty ptt_vector
                ptt_vector = np.empty((self.depth * 3, np.prod(self.shape)))

                # loop over the masks and fill in the masked ptt_vectors
                for mask in np.arange(self.depth):
                    ptt_vector[3*mask:3*mask+3] = unmasked_ptt_vector * self.mask[mask].ravel()

        return ptt_vector

    def fit_tilt(self, phase=None, ptt_vector=None):
        """Fit and remove tilt from a phase via least squares.

        Parameters
        ----------
        phase : array_like or None
            Phase to fit. If None (default), the Plane's :attr:`~lentil.Plane.phase`
            attribute it used.

        ptt_vector : array_like or None
            Piston, tip, tilt basis used for fitting. If None (default), the Plane's
            :attr:`~lentil.Plane.ptt_vector` attribute it used.

        Returns
        -------
        phase_no_tilt : ndarray
            ``phase`` with tilt fitted and removed. If ``len(ptt_vector)//3 > 1``, the
            returned phase array will contain one entry for each triad in ``ptt_vector``.

        tilt : list
            List of :class:`~lentil.Tilt` objects representing the fit and removed tilt.

        Note
        ----
        This method is called during :func:`~lentil.propagate` when ``tilt='angle'``. If
        custom tilt-removal behavior is required, this method should be modified in a
        subclass.

        If tilt-removal should `always` be skipped, this method should be modified in a
        subclass to: ``return self.phase, []``

        """
        if phase is None:
            phase = self.phase

        if ptt_vector is None:
            ptt_vector = self.ptt_vector

        # There are a couple of cases where we don't have enough information to remove the
        # tilt, so we just return the Plane as-is
        if ptt_vector is None or phase.size == 1:
            return

        nseg = len(ptt_vector)//3

        if nseg == 1:
            t = np.linalg.lstsq(ptt_vector.T, phase.ravel(), rcond=None)[0]
            phase_tilt = np.einsum('ij,i->j', ptt_vector[1:3], t[1:3])
            phase_no_tilt = phase - phase_tilt.reshape(phase.shape)

            # 01da13b transitioned from specifying tilt in terms of image plane
            # coordinates to about the Tilt plane axes. We can transform from
            # notional image plane to Tilt plane with x_tilt = -y_img, y_tilt = x_img
            tilt = [Tilt(x=-t[2], y=t[1])]

        else:
            t = np.empty((nseg, 3))
            phase_no_tilt = np.empty((nseg, phase.shape[0], phase.shape[1]))

            # iterate over the segments and compute the tilt term
            for seg in np.arange(nseg):
                t[seg] = np.linalg.lstsq(ptt_vector[3 * seg:3 * seg + 3].T, phase.ravel(),
                                         rcond=None)[0]
                seg_tilt = np.einsum('ij,i->j', ptt_vector[3 * seg + 1:3 * seg + 3], t[seg, 1:3])
                phase_no_tilt[seg] = (phase - seg_tilt.reshape(phase.shape)) * self.mask[seg]

            # 01da13b transitioned from specifying tilt in terms of image plane
            # coordinates to about the Tilt plane axes. We can transform from
            # notional image plane to Tilt plane with x_tilt = -y_img, y_tilt = x_img
            tilt = [Tilt(x=-t[seg, 2], y=t[seg, 1]) for seg in range(nseg)]

            phase_no_tilt = np.sum(phase_no_tilt, axis=0)

        self.phase = phase_no_tilt
        self.tilt.append(tilt)

    def slice(self, mask=None):
        """Compute slices defined by the data extent in ``mask``.

        Parameters
        ----------
        mask : array_like or None, optional
            Mask to use for identifying data extents. If None, the Plane's
            :attr:`mask` is used. If :attr:`mask` is None, ``Ellipsis`` is 
            returned (slice returns all data).

        Returns
        -------
        slices : list
            List of slices corresponding to the data extent defined by ``mask``.

        See also
        --------
        * :func:`~lentil.util.boundary_slice`
        * :func:`~lentil.Plane.slice_offset`

        """

        if mask is None:
            mask = self.mask

        # self.mask may still return None so we catch that here
        if mask.ndim < 2 or mask is None:
            # np.s_[...] = Ellipsis -> returns the whole array
            s = [np.s_[...]]
        elif mask.ndim == 2:
            s = [util.boundary_slice(mask)]
        elif mask.ndim == 3:
            s = [util.boundary_slice(m) for m in mask]
        return s

    def slice_offset(self, slices, shape=None, indexing='xy'):
        """Compute (r,c) offsets of the centers of each slice in a list of slices.

        Parameters
        ----------
        slices : list
            List of slices

        shape : (2,) array_like or None, optional
            Shape of containing array each slice is taken from. If None (default),
            ``shape = Plane.shape``.

        indexing : {'xy', 'ij'}, optional
            Offset ordering. Default is 'xy'.

        Returns
        -------
        offsets : list
           List of center offsets corresponding to the provided slices

        See also
        --------
        * :func:`~lentil.util.slice_offset`
        * :func:`~lentil.Plane.slice`

        """

        if shape is None:
            shape = self.shape

        offsets = []
        for s in slices:
            offset = util.slice_offset(slice=s, shape=shape, indexing=indexing)
            if offset:
                offsets.append(offset)
        if offsets:
            return offsets
        else:
            return [None]

    def multiply(self, wavefront):
        """Multiply with a wavefront

        Parameters
        ----------
        wavefront : :class:`~lentil.wavefront.Wavefront` object
            Wavefront to be multiplied

        Note
        ----
        It is possible to customize the way multiplication is performed by
        creating a subclass and overloading its ``multiply`` method.

        Returns
        -------
        wavefront : :class:`~lentil.wavefront.Wavefront` object
            Updated wavefront

        """

        if wavefront.pixelscale:
            # We may interpolate in the future, but for now we'll enforce equality            
            assert all(wavefront.pixelscale == self.pixelscale) 
        else:
            wavefront.pixelscale = self.pixelscale
            
        # broadcast existing wavefront data to Plane shape and clear the 
        # wavefront.data list
        data = [np.broadcast_to(d, self.shape) for d in wavefront.data]
        wavefront.data = []

        slc = self.slice()
        depth = self.depth

        # TODO: need to develop a strategy for handling pre-existing offsets
        # Still need to figure out how to deal with the fact that we may have an
        # incoming offset that may or may not be replaced by something new.
        #
        # this could be handled by a simple rule:
        #   new offset always takes prescedence over old offset UNLESS new offset is None and old offset is not none


        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # maybe enforce a simple rule
        # wavefront with offsets can only be multiplied with Planes that either introduce
        # NO offset or have the SAME offset

        offset = []

        for d in data:
            for n, slc in enumerate(self.slice()):
                mask = self.mask if depth == 1 else self.mask[n]
                p = _phasor(self.amplitude, self.phase, mask, wavefront.wavelength, slc)
                
                wavefront.data.append(d[slc] * p)

                ofst = util.slice_offset(slc, self.amplitude.shape)

                if ofst:
                    offset.append(ofst)
                else:
                    offset.append(None)

        if not wavefront.offset or not all(wavefront.offset):
            # of wavefront.offset is [] or a list of [None]:
             wavefront.offset = offset

        # Use extend instead of append so that nothing happens if self.shift = []
        wavefront.shift.extend(self.tilt)

        return wavefront


def _phasor(amplitude, phase, mask, wavelength, slc=Ellipsis):
    amp = amplitude if amplitude.size == 1 else amplitude[slc]
    ph = phase if phase.size == 1 else phase[slc]
    msk = mask if mask.size == 1 else mask[slc]
    return (amp * util.expc(2*np.pi*ph/wavelength)) * msk


class Pupil(Plane):
    """Class for representing a pupil plane.

    Parameters
    ----------
    diameter : float
        Diameter in meters

    focal_length : float
        Focal length in meters

    pixelscale : float
        Physical sampling (in meters) of each pixel in the pupil

    amplitude : array_like, optional
        Electric field amplitude transmission. Amplitude should be normalized
        with :func:`~lentil.util.normalize_power` if conservation of power
        through a diffraction propagation is required. If not specified, a
        default amplitude is created which has no effect on wavefront
        propagation.

    phase : array_like, optional
        Phase change caused by plane. If not specified, a default phase is
        created which has no effect on wavefront propagation.

    mask : array_like, optional
        Binary mask. If not specified, a mask is created from the amplitude.

    Note
    ----
    By definition, a pupil is represented by a spherical wavefront. Any
    aberrations in the optical system appear as deviations from this perfect
    sphere. The primary use of :class:`Pupil` is to represent these aberrations

    """
    def __init__(self, diameter=None, focal_length=None, pixelscale=None, amplitude=1,
                 phase=0, mask=None):

        super().__init__(pixelscale=pixelscale, amplitude=amplitude, phase=phase,
                         mask=mask)

        # We directly set the local attributes here in case a subclass has redefined
        # the property (which could cause an weird behavior and will throw an
        # AttributeError if the subclass hasn't defined an accompanying getter
        self._diameter = diameter
        self._focal_length = focal_length

    def __init_subclass__(cls):
        cls._diameter = None
        cls._focal_length = None

    @property
    def diameter(self):
        """Optical system diameter in meters."""
        return self._diameter

    @diameter.setter
    def diameter(self, value):
        self._diameter = value

    @property
    def focal_length(self):
        """Optical system focal length in meters."""
        return self._focal_length

    @focal_length.setter
    def focal_length(self, value):
        self._focal_length = value

    @property
    def f_number(self):
        """F-number."""
        return self.focal_length/self.diameter

    def multiply(self, wavefront):

        wavefront = super().multiply(wavefront)

        # we inherit the plane's focal length as the wavefront's focal length
        wavefront.focal_length = self.focal_length

        return wavefront


class Image(Plane):
    """Class for representing an image plane.

    Parameters
    ----------
    pixelscale : float, optional
        Pixel size in meters. Pixels are assumed to be square. Default is None.

    shape : {int, (2,) array_like}, optional
        Number of pixels as (rows, cols). If a single value is provided,
        :class:`Image` is assumed to be square with nrows = ncols = shape.
        Default is None.

    """

    def __init__(self, pixelscale=None, shape=None):

        super().__init__(pixelscale=pixelscale)

        # Ensure shape is a (2,) tuple or None
        if shape:
            shape = np.asarray(shape)
            if shape.size == 1:
                shape = np.append(shape, shape)
            elif shape.size > 2:
                raise IndexError('too many values')
            shape = tuple(shape)

        self._shape = shape

    @property
    def shape(self):
        """Number of pixels as (rows, cols)."""
        return self._shape

    @shape.setter
    def shape(self, value):
        self._shape = value

    def fit_tilt(self, *args, **kwargs):
        return self

    def slice(self, *args, **kwargs):
        # np.s_[...] = Ellipsis -> returns the whole array
        return [np.s_[...]]

class DispersivePhase(Plane):

    def multiply(self, wavefront):
        # NOTE: we can handle wavelength-dependent phase terms here (e.g. chromatic
        # aberrations). Since the phase will vary by wavelength, we can't fit out the
        # tilt pre-propagation and apply the same tilt for each wavelength like we can
        # with run of the mill tilt
        raise NotImplementedError


class DispersiveShift(Plane):

    def shift(self, wavelength, x0, y0, **kwargs):
        raise NotImplementedError

    def multiply(self, wavefront):
        wavefront = super().multiply(wavefront)
        wavefront.shift.extend([self])
        return wavefront


class Grism(DispersiveShift):
    r"""Class for representing a grism.

    A grism is an optical element that can be inserted into a collimated beam
    to disperse incoming light according to its wavelength.

    Light is dispersed along a line called the the spectral trace. The position
    along the trace is determined by the dispersion function of the grism. The
    local origin of the spectral trace is anchored relative to the undispersed
    position of the source. This basic geometry is illustrated in the figure
    below:

    .. image:: /_static/img/grism_geometry.png
        :align: center
        :width: 400px

    The spectral trace is parameterized by a polynomial of the form

    .. math::

        y = a_n x^n + \cdots + a_2 x^2 + a_1 x + a_0

    and should return units of meters on the focal plane provided an input
    in meters on the focal plane.

    Similarly, the wavelength along the trace is parameterized by a
    polynomial of the form

    .. math::

        \lambda = a_n d^n + \cdots + a_2 d^2 + a_1 d + a_0

    and should return units of meters of wavelength provided an input distance
    along the spectral trace.

    Note
    ----
    Lentil supports trace and dispersion functions with any arbitrary polynomial
    order. While a simple analytic solution exists for modeling first-order trace
    and/or dispersion, there is no general solution for higher order functions.

    As a result, trace and/or dispersion polynomials with order > 1 are evaluated
    numerically. Although the effects are small, this approach impacts both the
    speed and precision of modeling grisms with higher order trace and/or
    dispersion functions. In cases where speed or accuracy are extremely important,
    a custom solution may be required.

    Parameters
    ----------
    trace : array_like
        Polynomial coefficients describing the spectral trace produced by the
        grism in decreasing powers (i.e. trace[0] represents the highest order
        coefficient and trace[-1] represents the lowest).

    dispersion : array_like
        Polynomial coefficients describing the dispersion produced by the grism
        in decreasing powers (i.e. dispersion[0] represents the highest order
        coefficient and dispersion[-1] represents the lowest.)

    See Also
    --------

    """
    def __init__(self, trace, dispersion, pixelscale=None, amplitude=1,
                 phase=0, mask=None):
        super().__init__(pixelscale=pixelscale, amplitude=amplitude, phase=phase,
                         mask=mask)

        self.trace = np.asarray(trace)
        self._trace_order = self.trace.size - 1
        assert self._trace_order >= 1

        self.dispersion = np.asarray(dispersion)
        self._dispersion_order = self.dispersion.size - 1
        assert self._dispersion_order >= 1

    def shift(self, wavelength, xs=0., ys=0., **kwargs):

        dist = self._dispersion(wavelength)
        x, y = self._trace(dist)

        # Include any incoming shift
        x += xs
        y += ys

        return x, y

    def _dispersion(self, wavelength):
        """Compute distance along spectral trace given wavelength.

        This method automatically selects the optimal computation stragegy
        depending on the dispersion model order.

        Parameters
        ----------
        wavelength: float
            Wavelength to compute dispersed distance in meters.

        Returns
        -------
        distance: float
            Distance along spectral trace relative to self.dispersion[-1]

        """

        if self._dispersion_order == 1:
            # For a first order polynomial we have lambda = dispersion[0] * dist + dispersion[1]
            # Solving for distance gives dist = (lambda - dispersion[1])/dispersion[0]
            return (wavelength - self.dispersion[1]) / self.dispersion[0]
        else:
            # Compute the arc length by numerically integrating from lambda_ref (dispersion[-1])
            # to wavelength
            #return self._arc_len(self._dispersion_dist_func, self.dispersion[-1], wavelength)
            return scipy.optimize.leastsq(self._dist_cost_func, x0=0, args=(wavelength,))[0]

    def _trace(self, dist):

        if self._trace_order == 1:
            # https://www.physicsforums.com/threads/formula-for-finding-point-on-a-line-given-distance-along-a-line.419561/
            x = dist/np.sqrt(1+self.trace[0]**2)
        else:
            # Find x by matching the observed distance along the trace computed by
            # self._arc_len(self._trace_dist_func(x), 0, x) with the known distance
            # along the trace dist as provided from the wavelength (via self._dispersion)
            x = scipy.optimize.leastsq(self._trace_cost_func, x0=0, args=(dist,))[0]

        y = np.polyval(self.trace, x)

        return x, y

    def _dispersion_wavelength_func(self, dist):
        # Integrand for computing wavelength (via numerical integration of arc length formula)
        # along the polynomial defined by coefficients in self.dispersion
        #return np.sqrt(1 + np.polyval(np.polyder(self.dispersion), dist)**2)
        return np.polyval(self.dispersion, dist)

    def _dist_cost_func(self, x, wavelength):
        return wavelength - self._dispersion_wavelength_func(x)

    def _trace_dist_func(self, x):
        # Integrand for computing distance along the trace (via numerical integration  of arc
        # length formula) along the polynomial defined by coefficients in self.trace
        return np.sqrt(1 + np.polyval(np.polyder(self.trace), x)**2)

    def _trace_cost_func(self, x, dist):
        # Compute difference between dist and distance computed along trace given x
        return dist - self._arc_len(self._trace_dist_func, 0, x)

    @staticmethod
    def _arc_len(dist_func, a, b):
        """Compute arc length by numerically evaluating the following expression for arc length:

        s = integrate sqrt(1+f'(x)^2) dx = integrate dist_func dx

        where dist_func is either _dispersion_dist_func or _trace_dist_func and f(x) is the
        polynomial defining either the trace or dispersion.

        Parameters
        ----------
        dist_func: func
            Function defining the integrand to be numerically integrated
        a, b: float
            Lower and upper integration bounds

        Reference
        ---------
        https://en.wikipedia.org/wiki/Arc_length#Finding_arc_lengths_by_integrating

        """
        return scipy.integrate.quad(dist_func, a, b)[0]


class LensletArray(Plane):
    pass


class Tilt(Plane):
    """Object for representing tilt in terms of an angle

    Parameters
    ----------
    x : float
        Radians of tilt about the x-axis

    y : float
        Radians of tilt about the y-axis

    """
    def __init__(self, x, y):
        super().__init__()

        self.x = y  # y tilt is about the x-axis.
        self.y = -x  # x tilt is about the y axis. There's also a sign flip to get the direction right.

    def multiply(self, wavefront):
        wavefront = super().multiply(wavefront)
        wavefront.shift.extend([self])
        return wavefront

    def shift(self, xs=0, ys=0, z=0, **kwargs):
        """Compute image plane shift due to this angular tilt

        Parameters
        ----------
        xs : float
            Incoming x shift in meters. Default is 0.

        ys : float
            Incoming y shift in meters. Default is 0.

        z : float
            Propagation distance

        Returns
        -------
        shift : tuple
            Updated x and y shift terms

        """
        x = xs + (z * self.x)
        y = ys + (z * self.y)
        return x, y


class Rotate(Plane):
    """Rotate a Wavefront by a specified angle

    Parameters
    ----------
    angle : float
        Rotation angle counterclockwise from the horizontal.
    unit : {'degrees', 'radians'}, optional
        Units of angle. Default is 'degrees'.
    order : int
        The order of the spline interpolation (if needed), default is 3. The
        order has to be in the range 0-5.

    Note
    ----
    If the angle is an even multiple of 90 degrees, ``numpy.rot90`` is used to
    perform the rotation rather than ``scipy.ndimage.rotate``. In this case,
    the order parameter is irrelevant because no interpolation occurs.

    """
    def __init__(self, angle=0, unit='degrees', order=3):
        super().__init__()

        if unit == 'radians':
            angle *= 180/np.pi
        self.angle = -angle
        self.order = order

    def multiply(self, wavefront):
        """Multiply with a wavefront

        Parameters
        ----------
        wavefront : :class:`~lentil.wavefront.Wavefront` object
            Wavefront to be multiplied

        Returns
        -------
        wavefront : :class:`~lentil.wavefront.Wavefront` object
            Updated wavefront

        """
        if self.angle % 90 == 0:
            wavefront.data = np.rot90(wavefront.data, k=int(self.angle/90), axes=(1, 2))
        else:
            real = ndimage.rotate(wavefront.data.real, angle=self.angle, reshape=False, order=self.order, axes=(1, 2))
            imag = ndimage.rotate(wavefront.data.imag, angle=self.angle, reshape=False, order=self.order, axes=(1, 2))
            wavefront.data = real + 1j*imag
        return wavefront


class Flip(Plane):
    """Flip a wavefront along specified axis

    Parameters
    ----------
    axis : {int, tuple of ints, None}, optional
        Axis or axes along which to flip over. The default, axis=None, will
        flip over all of the axes of the input array. If axis is negative it
        counts from the last to the first axis.

    """

    def __init__(self, axis=None):
        super().__init__()

        if axis is None:
            # The first dimension of wavefront.data is depth so we actually want
            # to flip of the next two axes
            self.axis = (1, 2)
        else:
            # We convert axis to an array in case the user provides a tuple.
            # Again, because the first dimension of wavefront.data is depth we
            # add one to whatever axes are provided to make the flip operation
            # behave as expected.
            self.axis = np.asarray(axis) + 1

    def multiply(self, wavefront):
        """Multiply with a wavefront

        Parameters
        ----------
        wavefront : :class:`~lentil.wavefront.Wavefront` object
            Wavefront to be multiplied

        Returns
        -------
        wavefront : :class:`~lentil.wavefront.Wavefront` object
            Updated wavefront

        """
        wavefront.data = np.flip(wavefront.data, axis=self.axis)
        return wavefront


class Quadratic(Plane):
    """Base class for representing an optical plane with a quadratic phase
    term.

    """
    pass


class Conic(Plane):
    """Base class for representing an optical plane with a conic phase term.

    """
    pass


#    def Q(self, wave, pixelscale, oversample=1):
#        return (self.f_number*wave*oversample)/pixelscale
#
#    def alpha(self, wave, pixelscale, oversample=1):
#        return (self.pixelscale * pixelscale)/(self.focal_length * wave * oversample)
#
#    def q(self, wave, pixelscale, npix, oversample=1):
#        alpha = self.alpha(wave, pixelscale, oversample)
#        return 1/(alpha * npix)
