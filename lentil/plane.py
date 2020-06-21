import numpy as np
from scipy import ndimage

from lentil.cache import Cache
from lentil.modeltools import (cached_property, iterable_amplitude,
                                iterable_mask, iterable_phase, iterable_segmask)
from lentil import util

__all__ = ['Plane', 'Pupil', 'Image', 'Detector', 'DispersiveShift', 'Grism',
           'LensletArray', 'Tilt', 'Rotate', 'Flip']


class Plane:
    """Base class for representing a finite geometric plane.

    Parameters
    ----------
    pixelscale : float
        Physical sampling (in meters) of each pixel in the plane

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

    segmask : (nseg, ...) {array_like, None} optional
        Binary segment mask. Default is None.

    """

    def __init__(self, pixelscale=None, amplitude=1, phase=0, mask=None, segmask=None):

        # We directly set the local attributes here in case a subclass has redefined
        # the property (which could cause an weird behavior and will throw an
        # AttributeError if the subclass hasn't defined an accompanying getter
        self._pixelscale = pixelscale
        self._amplitude = np.asarray(amplitude) if amplitude is not None else None
        self._phase = np.asarray(phase) if phase is not None else None
        self._mask = np.asarray(mask) if mask is not None else None
        self._segmask = np.asarray(segmask) if segmask is not None else None

        self._cache = Cache()

    def __init_subclass__(cls):
        # we have to define default values to avoid AttributeErrors in case
        # subclasses don't call super().__init__ and don't explicitly define
        # an attribute
        cls._pixelscale = None
        cls._amplitude = np.array(1)
        cls._phase = np.array(0)
        cls._mask = None
        cls._segmask = None

        cls._cache = Cache()

    def __repr__(self):
        return self.name + '()'

    @property
    def pixelscale(self):
        """Physical sampling (in meters) of each pixel in the plane"""
        return self._pixelscale

    @pixelscale.setter
    def pixelscale(self, value):
        self._pixelscale = value

    @property
    def cache(self):
        """Cache

        Returns
        -------
        cache : :class:`~lentil.cache.Cache`

        """
        return self._cache

    @cached_property
    def amplitude(self):
        """Electric field amplitude transmission.

        Returns
        -------
        amplitude : ndarray

        Note
        ----
        ``amplitude`` is cacheable. If a cached value exists, it is returned.

        """
        return self._amplitude

    @amplitude.setter
    def amplitude(self, value):
        self.cache.delete('amplitude')
        if value is not None:
            self._amplitude = np.asarray(value)
        else:
            self._amplitude = None

    @cached_property
    def phase(self):
        """Electric field phase shift.

        Returns
        -------
        phase : ndarray

        Note
        ----
        ``phase`` is cacheable. If a cached value exists, it is returned.

        """
        return self._phase

    @phase.setter
    def phase(self, value):
        self.cache.delete('phase')
        if value is not None:
            self._phase = np.asarray(value)
        else:
            self._phase = None

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
            mask = self.amplitude
            mask[mask != 0] = 1
            return mask

    @mask.setter
    def mask(self, value):
        if value is not None:
            self._mask = np.asarray(value)
        else:
            self._mask = None

    @property
    def segmask(self):
        """Binary segment mask

        Returns
        -------
        segmask : (nseg,...) ndarray or None

        """
        return self._segmask

    @segmask.setter
    def segmask(self, value):
        # Note this operation implicitly converts a list or tuple of arrays
        # to an appropriately stacked ndarray
        if value is None:
            segmask = None
        else:
            # Note this operation implicitly converts a list or tuple of arrays
            # to an appropriately stacked ndarray
            segmask = np.asarray(value)
            if segmask.ndim == 2:
                segmask = segmask[np.newaxis, ...]
            elif segmask.ndim > 3:
                raise ValueError(f"don't know what to do with {segmask.ndim}th dimension")
        self._segmask = segmask

    @property
    def nseg(self):
        """Number of segments represented in the plane."""
        if self.segmask is None:
            return 0
        else:
            return self.segmask.shape[0]

    @property
    def shape(self):
        """Plane dimensions computed from :attr:`mask`. Returns None if :attr:`mask`
        is None.
        """
        if self.mask is None:
            return None
        else:
            return self.mask.shape

    @cached_property
    def ptt_vector(self):
        """2D vector representing piston and tilt in x and y. Planes with no mask or
        segmask have ptt_vector = None

        Returns
        -------
        ptt_vector : (3,...) ndarray or None

        """

        # if there's no mask, we just set ptt_vector to None and move on
        if self.shape == () and self.segmask is None:
            ptt_vector = None
        else:
            # compute unmasked piston, tip, tilt vector
            x, y = util.mesh(self.shape)
            unmasked_ptt_vector = np.array([np.ones(x.size), x.ravel(), y.ravel()]) * self.pixelscale

            if self.segmask is None:
                ptt_vector = np.einsum('ij,j->ij', unmasked_ptt_vector, self.mask.ravel())
            else:
                # prepare empty ptt_vector
                ptt_vector = np.empty((self.nseg*3, self.mask.size))

                # loop over the segments and fill in the masked ptt_vectors
                for seg in np.arange(self.nseg):
                    ptt_vector[3*seg:3*seg+3] = unmasked_ptt_vector * self.segmask[seg].ravel()

        return ptt_vector

    def cache_propagate(self):
        """Cache expensive to compute attributes for propagation."""
        self.cache.add('amplitude', self.amplitude)
        self.cache.add('phase', self.phase)
        self.cache.add('ptt_vector', self.ptt_vector)

    def clear_cache_propagate(self):
        """Clear propagation cache values."""
        self.cache.delete('amplitude')
        self.cache.delete('phase')
        self.cache.delete('ptt_vector')

    def multiply(self, wavefront, tilt='phase'):
        """Multiply with a wavefront

        Parameters
        ----------
        wavefront : :class:`~lentil.wavefront.Wavefront` object
            Wavefront to be multiplied

        tilt : {'phase', 'angle'}, optional
            * 'phase' - any tilt present in the Element phase contribution is
              included in the Element*Wavefront product. (Default)
            * 'angle' - any tilt present in the Element phase is removed before
              computing the Element*Wavefront product. The equivalent angular
              tilt is included in the Wavefront's
              :attr:`~lentil.Wavefront.tilt` attribute.

        Note
        ----
        It is possible to customize the way multiplication is performed in this
        class and any subclasses by reimplementing any of the following methods:

        * ``_multiply_phase`` is called when :attr:`tilt` is 'phase'
        * ``_multiply_angle`` is called when :attr:`tilt` is 'angle'

        In each case, the interface is defined as follows:

        ``wavefront = _multiply_method(amplitude, phase, wavefront)``

        Returns
        -------
        wavefront : :class:`~lentil.wavefront.Wavefront` object
            Updated wavefront

        """

        if wavefront.pixelscale is None:
            wavefront.pixelscale = self.pixelscale
        else:
            if self.pixelscale is not None:
                assert wavefront.pixelscale == self.pixelscale

        if tilt == 'phase':
            wavefront = self._multiply_phase(self.amplitude, self.phase, wavefront)
        elif tilt == 'angle':
            wavefront = self._multiply_angle(self.amplitude, self.phase, wavefront)
        else:
            raise AttributeError('Unknown tilt parameter', tilt)

        return wavefront

    @staticmethod
    def _multiply_phase(amplitude, phase, wavefront):
        # Compute and apply the phasor
        phasor = amplitude * np.exp(1j * phase * 2 * np.pi / wavefront.wavelength)
        wavefront.data *= phasor

        return wavefront

    def _multiply_angle(self, amplitude, phase, wavefront):

        # There are a couple of different paths here

        # If ptt_vector is None or the phase has no shape, there's nothing extra to do
        # here so we will just defer to _multiply_phase:
        if self.ptt_vector is None or phase.size == 1:
            wavefront = self._multiply_phase(amplitude, phase, wavefront)

        # If there's no segment mask, we can just do a single step projection,
        # subtraction, and multiplication:
        elif self.segmask is None:
            # TODO: move this into _multiple_angle_one (or something like that)
            ptt_vector = self.ptt_vector
            phase_vector = phase.ravel()
            t = np.linalg.lstsq(ptt_vector.T, phase_vector, rcond=None)[0]
            phase_tilt = \
                np.einsum('ij,i->j', ptt_vector[1:3], t[1:3]).reshape(self.shape)

            # Note that it is tempting to do an in-place subtraction here. Don't do
            # it! In the event that self.opd is returning a static ndarray, phase
            # will actually be a reference to this underlying array in memory and an
            # in-place operation will operate on that array. The in-place
            # operation's effects will persist beyond the scope of this function,
            # causing very strange behavior and lots of head scratching trying to
            # figure out what happened.
            phase = phase - phase_tilt

            # Set the wavefront angle
            wavefront.tilt.extend([Tilt(x=t[1], y=t[2])])

            # Compute and apply the phasor
            wavefront = self._multiply_phase(amplitude, phase, wavefront)

        # Otherwise we have to go segment by segment:
        else:
            # TODO: move this to _multiply_angle_multi
            ptt_vector = self.ptt_vector
            phase_vector = phase.ravel()
            t = np.zeros((self.nseg, 3))
            data = np.copy(wavefront.data)
            wavefront.data = np.zeros((self.nseg, data.shape[1], data.shape[2]),
                                      dtype=np.complex128)

            # iterate over the segments and compute the tilt term
            for seg in np.arange(self.nseg):
                t[seg] = np.linalg.lstsq(ptt_vector[3 * seg:3 * seg + 3].T, phase_vector,
                                         rcond=None)[0]
                seg_tilt = np.einsum('ij,i->j', ptt_vector[3 * seg + 1:3 * seg + 3],
                                     t[seg, 1:3]).reshape(self.shape)

                phasor = amplitude * np.exp(1j * (phase - seg_tilt) * 2 * np.pi / wavefront.wavelength)

                wavefront.data[seg] = data * phasor * self.segmask[seg]

            # Set the tilt term
            # Create a wavefront.Angle object for each segment and put them all in a
            # list
            wavefront.tilt.append([Tilt(x=t[seg, 1], y=t[seg, 2]) for seg in range(self.nseg)])

        return wavefront


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

    segmask : (nseg, ...) {array_like, None} optional
        Binary segment mask. Default is None.

    Note
    ----
    By definition, a pupil is represented by a spherical wavefront. Any
    aberrations in the optical system appear as deviations from this perfect
    sphere. The primary use of :class:`Pupil` is to represent these aberrations

    """
    def __init__(self, diameter=None, focal_length=None, pixelscale=None, amplitude=1,
                 phase=0, mask=None, segmask=None):

        super().__init__(pixelscale=pixelscale, amplitude=amplitude, phase=phase,
                         mask=mask, segmask=segmask)

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

    def multiply(self, wavefront, tilt='phase'):

        wavefront = super().multiply(wavefront, tilt)

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

    def multiply(self, wavefront, *args, **kwargs):
        """Multiply with a :class:`~lentil.wavefront.Wavefront`."""

        # TODO: make sure this is correct
        wavefront.focal_length = np.inf

        wavefront.planetype = 'image'

        return wavefront


class Detector(Image):
    """Class for representing a discretely sampled  :class:`~lentil.Image`
    plane.

    Parameters
    ----------
    pixelscale : float
        Pixel size in meters. Pixels are assumed to be square.

    shape : {int, (2,) array_like}
        Number of pixels as (rows, cols). If a single value is provided, the
        :class:`Detector` is assumed to be square with nrows = ncols = shape.

    """
    def __init__(self, pixelscale, shape):
        super().__init__(pixelscale, shape)


class DispersiveShift(Plane):

    def shift(self, wavelength, x0, y0, **kwargs):
        raise NotImplementedError

    def multiply(self, wavefront, tilt='angle'):
        wavefront = super().multiply(wavefront, tilt)

        wavefront.tilt.extend([self])

        return wavefront


class DispersivePhase(Plane):
    pass


class Grism(DispersiveShift):
    """Class for representing a grism."""

    def __init__(self, trace=None, dispersion=None, pixelscale=None, amplitude=1,
                 phase=0, mask=None, segmask=None):
        super().__init__(pixelscale=pixelscale, amplitude=amplitude, phase=phase,
                         mask=mask, segmask=segmask)

        self._trace = trace
        self._dispersion = dispersion

    def __init_subclass__(cls, **kwargs):
        cls._trace = None
        cls._dispersion = None

    @property
    def trace(self):
        return self._trace

    @trace.setter
    def trace(self, value):
        if value is None:
            self._trace = None
        else:
            self._trace = value

    @property
    def dispersion(self):
        return self._dispersion

    @dispersion.setter
    def dispersion(self, value):
        self._dispersion = value

    def shift(self, wavelength, xs=0., ys=0., **kwargs):

        # first we compute the distance d from the reference point to the
        # requested wavelength
        d = (wavelength - self.dispersion[1])/self.dispersion[0]

        # https://www.physicsforums.com/threads/formula-for-finding-point-on-a-line-given-distance-along-a-line.419561/

        # now we can compute the resulting coords due to the dispersion
        x = d/np.sqrt(1+self.trace[0]**2)
        y = np.polyval(self.trace, x)

        # finally, we include the incoming shift
        x += xs
        y += ys

        return x, y


class LensletArray(Plane):
    pass


class Tilt(Plane):
    """Object for representing tilt in terms of an angle

    Parameters
    ----------
    x : float
        x tilt in radians

    y : float
        y tilt in radians

    """
    def __init__(self, x, y):
        super().__init__()
        self.x = x
        self.y = y

    def multiply(self, wavefront, tilt):
        wavefront = super().multiply(wavefront, tilt)
        wavefront.tilt.extend([self])
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
        self.angle = angle
        self.order = order

    def cache_propagate(self):
        pass

    def clear_cache_propagate(self):
        pass

    def multiply(self, wavefront, *args, **kwargs):
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
            wavefront.data = np.rot90(wavefront.data, k=int(self.angle/90), axes=(1,2))
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

    def cache_propagate(self):
        pass

    def clear_cache_propagate(self):
        pass

    def multiply(self, wavefront, *args, **kwargs):
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
