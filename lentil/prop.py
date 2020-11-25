from operator import itemgetter
import copy

import numpy as np

from lentil import util
from lentil import fourier
from lentil.plane import Plane, Pupil, Image, Tilt
from lentil.wavefront import Wavefront

__all__ = ['propagate']


def propagate(planes, wave, weight=None, npix=None, npix_chip=None, oversample=2,
              rebin=True, tilt='phase', interp_phasor=True, flatten=True,
              use_multiprocessing=False):
    """Compute a polychromatic point spread function using Fraunhofer
    diffraction.

    Parameters
    ----------
    planes : list_like
        List of :class:`~lentil.Plane` objects

    wave : array_like
        Array of wavelengths (in meters)

    weight : array_like, optional
        Weight multiple applied to each wavelength slice in :attr:`wave`.
        :attr:`weight` can be relative (for example, when considering the
        optical transmission through the optical system) or absolute (for
        example, when performing radiometrically accurate propagations where
        wavelength-dependent flux at the image plane is known). Must have
        the same length as :attr:`wave`. If not specified, ones are used.

    npix : int or (2,) tuple of ints, optional
        Shape of output plane. If not specified,
        ``npix = OpticalSystem.planes[-1].shape``.

    npix_chip : int or (2,) tuple of ints, optional
        Shape of propagation output plane. If None (default),
        ``npix_chip = npix``. If ``npix_chip != npix``, the propagation
        result is placed in the appropriate location in the output plane.
        npix_chip cannot be larger than npix.

    oversample : int, optional
        Number of times to oversample the output plane. Default is 2.

    rebin : bool, optional
        If ``True``, return the output plane in the sampling given by
        ``pixelscale``, binning down the output plane by a factor of
        ``oversample`` as needed. Note that this operation preserves power
        in the output plane. Default is ``True``.

    tilt : {'phase', 'angle'}, optional
        Propagation tilt handling strategy

        * 'phase' - any tilt present in the Element phase contribution is
          included in the Element*Wavefront product. (Default)
        * 'angle' - any tilt present in the Element phase is removed before
          computing the Element*Wavefront product. The equivalent angular
          tilt is included in the Wavefront's
          :attr:`~lentil.Wavefront.tilt` attribute.

    interp_phasor : bool, optional
        If True (default), the phasor components will be automatically
        interpolated to avoid aliasing and wraparound in the detector plane.
        If False, no checking or interpolation is performed.

    flatten : bool, optional
        If ``True``, the cube of wavelength-dependent output planes is
        flattened into a single 2D array before being returned. If
        ``False``, a cube of output planes is returned. Default is ``True``.

    Returns
    -------
    psf : ndarray
        Resulting point spread function.

    """

    npix = _standardize_shape(npix, default=planes[-1].shape)
    npix_chip = _standardize_shape(npix_chip, default=npix)
    wave = _standardize_vec(wave)
    weight = _standardize_vec(weight, default=np.ones_like(wave))

    # What should happen here as a part of cache_propagate:
    # -----------------------------------------------------
    # * if fit_tilt = True, the tilt in each plane should be fit and removed
    #   from the phase and bookkept in the plane.tilt attribute
    #   - note that this can be done on the un-sliced data for ease of interaction
    #     with ptt_vector
    #   - note further that ptt_vector should be reworked as a property with a
    #     setter that is called at object creation with the mask in case we need
    #     to do some reinterpolation down the line
    #
    # *  Once this is done we can remove ptt_vector from the cache that is
    #    created in cache_propagate() since it won't be needed anymore. As
    #    a matter of fact, we can just turn this into a method that takes
    #    a mask or segmask and returns the vector -> _ptt_vector(mask)
    #
    # * Somewhere in here (probably before we fit and remove tilt), we need
    #   to sort out the sampling and make sure each plane is adequately
    #   sampled. Write up sampling rules for this.
    #
    #

    # note the actual order of operations is probably:
    # 1. sampling
    # 2. set slice
    # 3. fit_tilt (but this can be done over full arrays)

    _prepare_planes(planes, wave, npix_chip, oversample, tilt, interp_phasor)

    # Create an empty output
    if flatten:
        output = np.zeros((npix[0]*oversample, npix[1]*oversample))
    else:
        output = np.zeros((len(wave), npix[0]*oversample, npix[1]*oversample))

    # We also need a temporary array to place the chips and compute intensity
    _output = np.zeros((npix[0]*oversample, npix[1]*oversample), dtype=np.complex128)

    for n, (wl, wt) in enumerate(zip(wave, weight)):
        if wt > 0:

            _output = _propagate_mono(planes, wl, npix, npix_chip, oversample, tilt, _output)

            # Compute intensity
            if flatten:
                output += np.abs(_output)**2 * wt
            else:
                output[n] = np.abs(_output)**2 * wt

            # Zero out the local output array
            _output[:] = 0

    if rebin:
        output = util.rebin(output, oversample)

    _cleanup_planes(planes)

    return output


def _prepare_planes(planes, wave, npix, oversample, tilt, interp_phasor):
    #shapes = [plane.shape for plane in planes if isinstance(plane, Pupil)]
    #npix_wavefront = (max(shapes, key=itemgetter(0))[0],
    #                  max(shapes, key=itemgetter(1))[1])

    # TODO: we may be able to avoid the deepcopies here on amplitude and
    # phase if they are reshaped

    # Loop over the planes to build caches
    for plane in planes:
        plane.cache['amplitude'] = copy.deepcopy(plane.amplitude)
        plane.cache['tilt'] = copy.deepcopy(plane.tilt)

        if tilt == 'angle':
            phase, fit_tilt = _fit_tilt(plane)
            plane.cache['phase'] = phase
            if fit_tilt is not None:
                plane.cache['tilt'].extend([fit_tilt])
        
        else:
            plane.cache['phase'] = copy.deepcopy(plane.phase)


def _cleanup_planes(planes):
    for plane in planes:
        plane.cache.clear()



def _fit_tilt(plane):

    ptt_vector = plane.ptt_vector
    if ptt_vector is None or plane.phase.size == 1:
        # There's nothing to do so we'll just return the unmodified
        # phase and a no-tilt Tilt object
        return plane.phase, None

    elif plane.segmask is None:
        phase_vector = plane.phase.ravel()

        t = np.linalg.lstsq(ptt_vector.T, phase_vector, rcond=None)[0]
        phase_tilt = np.einsum('ij,i->j', ptt_vector[1:3], t[1:3])

        phase = plane.phase - phase_tilt.reshape(plane.shape)
        # 01da13b transitioned from specifying tilt in terms of image plane
        # coordinates to about the Tilt plane axes. We can transform from
        # notional image plane to Tilt plane with x_tilt = -y_img, y_tilt = x_img
        tilt = [Tilt(x=-t[2], y=t[1])]
        return phase, tilt

    else:
        return _fit_tilt_segmented(plane)


def _fit_tilt_segmented(plane):
    ptt_vector = plane.ptt_vector

    phase_vector = plane.phase.ravel()
    t = np.zeros((plane.nseg, 3))
    phase = np.zeros((plane.nseg, plane.shape[0], plane.shape[1]))

    # iterate over the segments and compute the tilt term
    for seg in np.arange(plane.nseg):
        t[seg] = np.linalg.lstsq(ptt_vector[3 * seg:3 * seg + 3].T, phase_vector,
                                         rcond=None)[0]
        seg_tilt = np.einsum('ij,i->j', ptt_vector[3 * seg + 1:3 * seg + 3], t[seg, 1:3])
        phase[seg] = (plane.phase - seg_tilt.reshape(plane.shape)) * plane.segmask[seg]

    # 01da13b transitioned from specifying tilt in terms of image plane
    # coordinates to about the Tilt plane axes. We can transform from
    # notional image plane to Tilt plane with x_tilt = -y_img, y_tilt = x_img
    #
    # Note the
    tilt = [Tilt(x=-t[seg, 2], y=t[seg, 1]) for seg in range(plane.nseg)]

    return phase, tilt


def _propagate_mono(planes, wavelength, npix, npix_chip, oversample, tilt, out=None):
    """Propagate a monochromatic wavefront from plane to plane through the
    optical system using Fraunhofer diffraction.

    Parameters
    ----------
    w : :class:`~lentil.Wavefront`
        Wavefront object to propagate through the optical system.

    npix : int or (2,) ints
        Shape of output plane.

    oversample : int
       Number of times to oversample the output plane.

    tilt : {'phase', 'angle'}
        * 'phase' - any tilt present in the Element phase contribution is
          included in the Element*Wavefront product.
        * 'angle' - any tilt present in the Element phase is removed
          before computing the Element*Wavefront product. The equivalent
          angular tilt is included in the Wavefront's
          :attr:`~lentil.Wavefront.tilt` attribute.

    Returns
    -------
    field : ndarray
        Resulting complex field propagated though the optical system.

    """
    if out is None:
        out = np.zeros((npix[0]*oversample, npix[1]*oversample), dtype=np.complex128)

    w = Wavefront(wavelength=wavelength, pixelscale=None,
                  shape=planes[0].shape, planetype=None)

    for plane, next_plane in _iterate_planes(planes):

        # Multiply by the current plane
        w = plane.multiply(w)

        # Now, we propagate to the next plane in the optical system

        # TODO: figure out how/when to accumulate tilt

        if (w.planetype == 'pupil') and isinstance(next_plane, Image):
            if next_plane.pixelscale is not None:
                w = _propagate_pti(w, next_plane.pixelscale, npix_chip, oversample)
            else:
                pass
        elif (w.planetype == 'image') and isinstance(next_plane, Pupil):
            pass
        elif (w.planetype == 'pupil') and isinstance(next_plane, Pupil):
            continue
        elif (w.planetype == 'image') and isinstance(next_plane, Image):
            continue
        elif isinstance(next_plane, Plane):
            continue
        else:
            raise TypeError('Unsupported propagation type ', w.planetype, ' to ', next_plane)

    for d in range(w.depth):

        # The shift term is given in terms of (x,y) but we place the chip in
        # terms of (r,c)
        # TODO: I think this will break if the last plane isn't a Detector
        shift = np.flip(w.pixel_shift[d], axis=0)

        # Compute the chip location
        canvas_slice, chip_slice = _chip_insertion_slices(out.shape,
                                                          (w.data.shape[1], w.data.shape[2]),
                                                          shift)

        # Insert the complex result in the output
        if canvas_slice:
            out[canvas_slice] += w.data[d, chip_slice[0], chip_slice[1]]

    return out


def _propagate_pti(w, pixelscale, npix, oversample):

    # TODO: we should only apply the shift if this is the final plane
    shift = w.shift(pixelscale, oversample)
    w.tilt = []

    # Integer portion of the shift that will be accounted for
    # later
    w.pixel_shift = np.fix(shift)

    # Residual subpixel portion of the shift that is passed to
    # the DFT
    res_shift = shift - w.pixel_shift

    npix = npix * oversample

    alpha = (w.pixelscale * pixelscale) / (w.wavelength * w.focal_length * oversample)

    data = np.empty((w.depth, npix[0], npix[1]), dtype=np.complex128)

    assert shift.shape[0] == w.depth, \
        'Dimension mismatch between tilt and wavefront depth'

    for d in range(w.depth):
        # data[d] = fourier.dft2(w.data[d], alpha, npix, res_shift[d])
        data[d] = fourier.dft2(w.data[d], alpha, npix, res_shift[d], unitary= True, out=data[d])

    w.data = data

    return w


def _standardize_shape(shape, default=()):
    if shape is None:
        shape = default
    shape = np.asarray(shape)
    if shape.shape == ():
        shape = np.append(shape, shape)
    return shape

def _standardize_vec(vec, default=[]):
    if vec is None:
        vec = default
    vec = np.asarray(vec)
    if vec.shape == ():
        vec = vec[np.newaxis, ...]
    return vec

class _iterate_planes:
    def __init__(self, planes):
        self.planes = planes
        self.length = len(planes)
        self.n = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.n < self.length-1:
            plane = self.planes[self.n]
            next_plane = self.planes[self.n+1]
            self.n += 1
            return plane, next_plane
        else:
            raise StopIteration()

class _iterate_planes_reverse:
    def __init__(self, planes):
        self.planes = planes
        self.length = len(planes)
        self.n = 1

    def __iter__(self):
        return self

    def __next__(self):
        if self.n < self.length:
            plane = self.planes[-(self.n+1)]
            next_plane = self.planes[-self.n]
            self.n += 1
            return plane, next_plane
        else:
            raise StopIteration()

def _chip_insertion_slices(npix_canvas, npix_chip, shift):
    npix_canvas = np.asarray(npix_canvas)
    npix_chip = np.asarray(npix_chip)

    # Canvas coordinates of the upper left corner of the shifted chip
    chip_shifted_ul = (npix_canvas / 2) - (npix_chip / 2) + shift

    # Chip slice indices
    chip_top = int(0)
    chip_bottom = int(npix_chip[0])
    chip_left = int(0)
    chip_right = int(npix_chip[1])

    # Canvas insertion slice indices
    canvas_top = int(chip_shifted_ul[0])
    canvas_bottom = int(chip_shifted_ul[0] + npix_chip[0])
    canvas_left = int(chip_shifted_ul[1])
    canvas_right = int(chip_shifted_ul[1] + npix_chip[1])

    # reconcile the chip and canvas insertion indices
    if canvas_top < 0:
        chip_top = -1 * canvas_top
        canvas_top = 0

    if canvas_bottom > npix_canvas[0]:
        chip_bottom -= canvas_bottom - npix_canvas[0]
        canvas_bottom = npix_canvas[0]

    if canvas_left < 0:
        chip_left = -1 * canvas_left
        canvas_left = 0

    if canvas_right > npix_canvas[1]:
        chip_right -= canvas_right - npix_canvas[1]
        canvas_right = npix_canvas[1]

    if np.any(np.array([canvas_bottom, chip_bottom, canvas_right, chip_right]) < 0):
        return None, None
    else:
        return (slice(canvas_top, canvas_bottom), slice(canvas_left, canvas_right)), \
               (slice(chip_top, chip_bottom), slice(chip_left, chip_right))
