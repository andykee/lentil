import copy

import numpy as np
from scipy import ndimage
import scipy.integrate
import scipy.optimize

import lentil
from lentil.field import Field
import lentil.helper


class Plane:
    """
    Base class for representing a finite geometric plane.

    Parameters
    ----------
    amplitude : array_like, optional
        Electric field amplitude transmission. Amplitude should be normalized
        with :func:`~lentil.normalize_power` if conservation of power
        through a diffraction propagation is required. If not specified, a
        default amplitude is created which has no effect on wavefront
        propagation.
    phase : array_like, optional
        Phase change caused by plane. If not specified, a default phase is
        created which has no effect on wavefront propagation.
    mask : array_like, optional
        Binary mask. If not specified, a mask is created from the amplitude.
        If ``mask`` has 2 dimensions, the plane is assumed to be monolithic. If
        ``mask`` has 3 dimensions, the plane is assumed to be segmented with the
        individual segment masks inserted along the first dimension.

        .. image:: /_static/img/segmask.png
            :width: 550px
            :align: center

    pixelscale : float or (2,) array_like
        Physical sampling (in meters) of each pixel in the plane. If :attr:`pixelscale`
        is a scalar, uniform sampling in x and y is assumed.
    """
    def __init__(self, amplitude=1, phase=0, mask=None, pixelscale=None):
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
        cls._tilt = []

    def __repr__(self):
        return f'{self.__class__.__name__}()'

    @property
    def amplitude(self):
        """
        Electric field amplitude transmission

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
        """
        Electric field phase shift

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

        If ``mask`` has 2 dimensions, the plane is assumed to be monolithic. If
        ``mask`` has 3 dimensions, the plane is assumed to be segmented with the
        individual segment masks inserted along the first dimension.

        .. image:: /_static/img/segmask.png
            :width: 550px
            :align: center

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
        """
        Flattened view of :attr:`mask`

        Returns
        -------
        mask : ndarray
        """
        if self.depth < 2:
            return self.mask
        else:
            return np.sum(self.mask, axis=0)

    @property
    def pixelscale(self):
        """
        Physical (row, col) sampling of each pixel in the Plane.

        Returns
        -------
        pixelscale : ndarray
        """
        return self._pixelscale

    @pixelscale.setter
    def pixelscale(self, value):
        if value is not None:
            self._pixelscale = lentil.helper.sanitize_shape(value)
        else:
            self._pixelscale = None

    @property
    def tilt(self):
        """
        List of additional :class:`~lentil.Tilt` terms associated wirth this Plane

        Returns
        -------
        tilt : list

        See Also
        --------
        Tilt : Tilt object
        Plane.fit_tilt : Method for fitting and removing tilt
        """
        return self._tilt

    @tilt.setter
    def tilt(self, value):
        self._tilt = value

    @property
    def shape(self):
        """
        Plane dimensions computed from :attr:`mask`.

        Returns (mask.shape[1], mask.shape[2]) if mask has ``ndim == 3``. Returns
        None if :attr:`mask` is None.
        """
        if self.depth == 1:
            return self.mask.shape
        else:
            return self.mask.shape[1], self.mask.shape[2]

    @property
    def depth(self):
        """
        Number of independent masks (segments) in :attr:`mask`

        Returns
        -------
        depth : int
        """
        if self.mask.ndim in (0, 1, 2):
            return 1
        else:
            return self.mask.shape[0]

    @property
    def ptt_vector(self):
        """
        2D vector representing piston and tilt in x and y.

        Planes with no mask have ``ptt_vector = None``.

        Returns
        -------
        ptt_vector : ndarray or None
        """

        # if there's no mask, we just set ptt_vector to None and move on
        if self.shape == () or self.shape is None:
            ptt_vector = None
        else:
            # compute unmasked piston, x-tilt, y-tilt vector
            r, c = lentil.helper.mesh(self.shape)
            # note c is negative here. this is done to match Lentil's
            # definition of what +x tilt and +y tilt look like. The columns of
            # unmasked_ptt_vector should match the shapes returned by lentil.zernike for
            # modes 1:3
            unmasked_ptt_vector = np.einsum('ij,i->ij', [np.ones(r.size), r.ravel(), -c.ravel()],
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

    def copy(self):
        """
        Make a copy of this object.

        Returns
        -------
        copy : :class:`~lentil.Plane`
        """
        return copy.deepcopy(self)

    def fit_tilt(self, inplace=True):
        """
        Fit and remove tilt from Plane :attr:`phase` via least squares. The
        equivalent angular tilt is bookkept in Plane :attr:`tilt`.

        Parameters
        ----------
        inplace : bool, optional
            Modify the original object in place (True) or create a copy (False,
            default).

        Returns
        -------
        plane : :class:`~lentil.Plane`
        """
        if inplace:
            plane = self
        else:
            plane = self.copy()

        ptt_vector = plane.ptt_vector

        # There are a couple of cases where we don't have enough information to remove the
        # tilt, so we just return the Plane as-is
        if ptt_vector is None or plane.phase.size == 1:
            return plane

        if self.depth == 1:
            t = np.linalg.lstsq(ptt_vector.T, plane.phase.ravel(), rcond=None)[0]
            phase_tilt = np.einsum('ij,i->j', ptt_vector[1:3], t[1:3])
            plane.phase -= phase_tilt.reshape(plane.phase.shape)
            plane.tilt.append(Tilt(x=t[1], y=t[2]))

        else:
            t = np.empty((self.depth, 3))
            phase_no_tilt = np.empty((self.depth, plane.phase.shape[0], plane.phase.shape[1]))

            # iterate over the segments and compute the tilt term
            for seg in np.arange(self.depth):
                t[seg] = np.linalg.lstsq(ptt_vector[3 * seg:3 * seg + 3].T, plane.phase.ravel(),
                                         rcond=None)[0]
                seg_tilt = np.einsum('ij,i->j', ptt_vector[3 * seg + 1:3 * seg + 3], t[seg, 1:3])
                phase_no_tilt[seg] = (plane.phase - seg_tilt.reshape(plane.phase.shape)) * self.mask[seg]

            plane.phase = np.sum(phase_no_tilt, axis=0)
            plane.tilt.extend([Tilt(x=t[seg, 1], y=t[seg, 2]) for seg in range(self.depth)])

        return plane

    def rescale(self, scale, inplace=True):
        """
        Rescale a plane via interpolation.

        The following Plane attributes are resampled:

        * `Plane.amplitude` is rescaled via 3rd order spline interpolation. The
            result is scaled to preserve total power.
        * `Plane.phase` is rescaled via 3rd order spline interpolation.
        * `Plane.mask` is rescaled via 0-order nearest neighbor interpolation.
        * `Plane.pixelscale` is adjusted appropriately.

        Parameters
        ----------
        scale : float
            Scale factor for interpolation. Scale factors less than 1 shrink the
            Plane while scale factors greater than 1 grow it.
        inplace : bool, optional
            Modify the original object in place (True, default) or create a
            copy (False).

        Returns
        -------
        plane : :class:`Plane`

        Note
        ----
        All interpolation is performed via `scipy.ndimage.map_coordinates`

        See Also
        --------
        Plane.resample

        """
        if inplace:
            plane = self
        else:
            plane = self.copy()

        if plane.amplitude.ndim > 1:
            plane.amplitude = lentil.rescale(plane.amplitude, scale=scale, shape=None,
                                                mask=None, order=3, mode='nearest',
                                                unitary=False)/scale

        if plane.phase.ndim > 1:
            plane.phase = lentil.rescale(plane.phase, scale=scale, shape=None, mask=None,
                                         order=3, mode='nearest', unitary=False)

        # because plane.mask is automatically computed from amplitude if it is not
        # provided, we really only want to reinterpolate if a mask was provided (stored
        # in plane._mask) vs computed on the fly (by the plane.mask property)
        if plane._mask is not None:
            if plane.mask.ndim > 1:
                plane.mask = lentil.rescale(plane.mask, scale=scale, shape=None, mask=None, order=0,
                                            mode='constant', unitary=False)
            # mask[mask < np.finfo(mask.dtype).eps] = 0
            plane.mask[np.nonzero(plane.mask)] = 1
            plane.mask = plane.mask.astype(int)

        plane.pixelscale = plane.pixelscale/scale

        return plane

    def resample(self, pixelscale, inplace=True):
        """Resample a plane via interpolation.

        The following Plane attributes are resampled:

        * `Plane.amplitude` is resampled via 3rd order spline interpolation. The
          result is scaled to preserve total power.
        * `Plane.phase` is resampled via 3rd order spline interpolation.
        * `Plane.mask` is resampled via 0-order nearest neighbor interpolation.
        * `Plane.pixelscale` is adjusted appropriately.

        Parameters
        ----------
        pixelscale : float
            Desired Plane pixelscale.
        inplace : bool, optional
            Modify the original object in place (True, default) or create a
            copy (False).

        Returns
        -------
        plane : :class:`Plane`
            Resampled Plane.

        Note
        ----
        All interpolation is performed via `scipy.ndimage.map_coordinates`

        See Also
        --------
        * :func:`Plane.rescale`

        """
        if self.pixelscale[0] != self.pixelscale[1]:
            raise NotImplementedError("Can't resample non-uniformly sampled Plane")

        return self.rescale(scale=self.pixelscale[0]/pixelscale, inplace=inplace)

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
        * :func:`~lentil.helper.boundary_slice`
        * :func:`~lentil.Plane.slice_offset`

        """
        if mask is None:
            mask = self.mask

        # self.mask may still return None so we catch that here
        if mask.ndim < 2 or mask is None:
            # np.s_[...] = Ellipsis -> returns the whole array
            s = [np.s_[...]]
        elif mask.ndim == 2:
            s = [lentil.helper.boundary_slice(mask)]
        elif mask.ndim == 3:
            s = [lentil.helper.boundary_slice(m) for m in mask]
        else:
            raise ValueError('mask has invalid dimensions')
        return s

    # TODO: add out=None to actually handle in-place multiplication
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
        pixelscale = lentil.field.multiply_pixelscale(self.pixelscale, wavefront.pixelscale)
        out = lentil.Wavefront(wavelength=wavefront.wavelength,
                               pixelscale=pixelscale,
                               planetype=wavefront.planetype,
                               focal_length=wavefront.focal_length,
                               data=[])

        slc = self.slice()
        for field in wavefront.data:
            for n, s in enumerate(slc):
                amp = self.amplitude if self.amplitude.size == 1 else self.amplitude[s]
                phase = self.phase if self.phase.size == 1 else self.phase[s]

                # construct complex phasor
                phasor = Field(data=amp*np.exp(-2*np.pi*1j*phase/wavefront.wavelength),
                               pixelscale=self.pixelscale,
                               offset=lentil.helper.slice_offset(s, self.shape),
                               tilt=[self.tilt[n]] if self.tilt else [])

                out.data.append(field * phasor)

        return out


class Pupil(Plane):
    """Class for representing a pupil plane.

    Parameters
    ----------
    focal_length : float
        Focal length
    pixelscale : float
        Physical sampling of each pixel in the pupil
    amplitude : array_like, optional
        Electric field amplitude transmission. Amplitude should be normalized
        with :func:`~lentil.normalize_power` if conservation of power
        through a diffraction propagation is required. If not specified, a
        default amplitude is created which has no effect on wavefront
        propagation.
    phase : array_like, optional
        Phase change caused by plane. If not specified, a default phase is
        created which has no effect on wavefront propagation.
    mask : array_like, optional
        Binary mask. If not specified, a mask is created from the amplitude.
        If ``mask`` has 2 dimensions, the plane is assumed to be monolithic. If
        ``mask`` has 3 dimensions, the plane is assumed to be segmented with the
        individual segment masks inserted along the first dimension.

        .. image:: /_static/img/segmask.png
            :width: 550px
            :align: center

    Note
    ----
    By definition, a pupil is represented by a spherical wavefront. Any
    aberrations in the optical system appear as deviations from this perfect
    sphere. The primary use of :class:`Pupil` is to represent these aberrations

    """
    def __init__(self, focal_length=None, pixelscale=None, amplitude=1,
                 phase=0, mask=None):

        super().__init__(pixelscale=pixelscale, amplitude=amplitude, phase=phase,
                         mask=mask)

        # We directly set the local attributes here in case a subclass has redefined
        # the property (which could cause an weird behavior and will throw an
        # AttributeError if the subclass hasn't defined an accompanying getter
        self._focal_length = focal_length

    def __init_subclass__(cls):
        cls._focal_length = None

    @property
    def focal_length(self):
        """Pupil focal length"""
        return self._focal_length

    @focal_length.setter
    def focal_length(self, value):
        self._focal_length = value

    def multiply(self, wavefront):

        wavefront = super().multiply(wavefront)

        # we inherit the plane's focal length as the wavefront's focal length
        wavefront.focal_length = self.focal_length
        wavefront.planetype = 'pupil'

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

    Note
    ----
    If image plane intensity is desired, significant performance improvements
    can be realized by using a :class:`Detector` plane instead.

    See Also
    --------
    Detector

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

    def multiply(self, wavefront):
        wavefront = super().multiply(wavefront)
        wavefront.planetype = 'image'
        return wavefront


class Detector(Image):
    """Class for representing an image plane that returns intensity.

    The Detector should only be used as the last plane in a propagation. If
    individual wavelength results or access to the complex field is required, an
    :class:`Image` plane should be used instead.

    Parameters
    ----------
    pixelscale : float, optional
        Pixel size in meters. Pixels are assumed to be square. Default is None.
    shape : {int, (2,) array_like}, optional
        Number of pixels as (rows, cols). If a single value is provided,
        :class:`Image` is assumed to be square with nrows = ncols = shape.
        Default is None.

    See Also
    --------
    Image

    """
    def __init__(self, pixelscale=None, shape=None):
        super().__init__(pixelscale=pixelscale, shape=shape)


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

        References
        ----------
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
        self.y = x  # x tilt is about the y axis.

    def multiply(self, wavefront):
        wavefront = super().multiply(wavefront)
        for field in wavefront.data:
            field.tilt.append(self)
        return wavefront

    def shift(self, xs=0, ys=0, z=0, **kwargs):
        """Compute image plane shift due to this angular tilt

        Parameters
        ----------
        xs : float
            Incoming x shift distance. Default is 0.
        ys : float
            Incoming y shift distance. Default is 0.
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
