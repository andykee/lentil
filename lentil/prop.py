from operator import itemgetter
import copy

import numpy as np

from lentil import util
from lentil import fourier
from lentil.plane import Plane, Pupil, Image, Detector, Tilt
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
    wave = _standardize_bandpass(wave)
    weight = _standardize_bandpass(weight, default=np.ones_like(wave))

    # TODO: this is temporarily commented out for v0.5.0 beta testing
    # if (rebin is True) and (not isinstance(planes[-1], Detector)):
    #     raise AttributeError('rebin is only supported when propagating to a Detector plane')

    _prepare_planes(planes, wave, npix_chip, oversample, tilt, interp_phasor)

    # Create an empty output

    # TODO: temporary image = detector workaround for v0.5.0 beta testing
    output_dtype = np.float64
    #output_dtype = np.float64 if isinstance(planes[-1], Detector) else np.complex128
    output_shape = (npix[0]*oversample, npix[1]*oversample)

    if flatten:
        output = np.zeros(output_shape, dtype=output_dtype)
    else:
        output = np.zeros((len(wave), output_shape[0], output_shape[1]), dtype=output_dtype)

    for n, (wl, wt) in enumerate(zip(wave, weight)):
        if wt > 0:

            w = _propagate_mono(planes, wl, npix_chip, oversample)

            if w.planetype == 'image':

                for d in range(w.depth):
                    # The shift term is given in terms of (x,y) but we place the chip in
                    # terms of (r,c)
                    shift = np.flip(w.center[d], axis=0)

                    # Compute the chip location
                    data_slice, chip_slice = _chip_insertion_slices(output_shape,
                                                                    (w.data[d].shape[0], w.data[d].shape[1]),
                                                                    shift)

                    # Insert the complex result in the output
                    if data_slice:
                        if flatten:
                            output[data_slice] += w.data[d][chip_slice] * wt
                        else:
                            output[n][data_slice] = w.data[d][chip_slice] * wt




            # At this point w.data should always have len == 1
            #if flatten:
            #    output += w.data[0] * wt
            #else:
            #    output[n] = w.data[0] * wt

    if rebin:
        output = util.rebin(output, oversample)

    _cleanup_planes(planes)

    return output


def _prepare_planes(planes, wave, npix, oversample, tilt, interp_phasor):

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
    # TODO: rework this
    # Compute the wavefront shape required to support propagation through the
    # provided list of planes

    # TODO: not sure where we should reinterpolate and how it should be done.

    # The current approach is to find the maximum row and column dimensions of
    # only the Pupil elements in self.planes. This whole approach will have to
    # change once we support more complicated multi-plane propagations.
    #
    # I think what we should really do is find the shape of the first Plane
    # of consequence. What we have hwre works for now though
    shapes = [plane.shape for plane in planes if isinstance(plane, Pupil)]
    npix_wavefront = (max(shapes, key=itemgetter(0))[0],
                      max(shapes, key=itemgetter(1))[1])

    # Loop over the planes to build caches
    for plane in planes:
        # handle the tilt first
        plane.cache['tilt'] = copy.deepcopy(plane.tilt)
        if tilt == 'angle':
            phase, fit_tilt = plane.fit_tilt()
            plane.cache['phase'] = phase
            if fit_tilt:
                plane.cache['tilt'].extend([fit_tilt])
        else:
            plane.cache['phase'] = copy.deepcopy(plane.phase)

        plane.cache['amplitude'] = copy.deepcopy(plane.amplitude)
        plane.cache['pixelscale'] = copy.deepcopy(plane.pixelscale)
        plane.cache['npix_wavefront'] = npix_wavefront
        slc = plane.slice()
        plane.cache['slice'] = slc
        plane.cache['offset'] = plane.slice_offset(slices=slc, shape=plane.shape, indexing='xy')

        # reinterpolate as needed
        # NOTE: if we reinterpolate, we'll have to make sure we provide the correct
        # slices and offsets (using the appropriate reinterpolated mask, probably)
        if interp_phasor:
            # We'll have already computed the interpolation scale when
            # figuring out npix_wavefront so we just need to rescale
            # amplitude and phase here
            pass


def _cleanup_planes(planes):
    for plane in planes:
        plane.cache.clear()


def _propagate_mono(planes, wavelength, npix_chip, oversample):
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

    w = Wavefront(wavelength=wavelength, pixelscale=None,
                  shape=planes[0].cache['npix_wavefront'], planetype=None)

    for plane in planes:

        # Propagate to plane
        if w.planetype is None:
            # free space propagation to first plane
            pass
        elif (w.planetype == 'pupil') and isinstance(plane, Image):
            if plane.pixelscale is not None:
                w = _propagate_pti(w, plane.pixelscale, npix_chip, oversample)
            else:
                pass
        elif (w.planetype == 'image') and isinstance(plane, Pupil):
            pass
        elif (w.planetype == 'pupil') and isinstance(plane, Pupil):
            pass
        elif (w.planetype == 'image') and isinstance(plane, Image):
            pass
        elif isinstance(plane, Plane):
            pass
        else:
            raise TypeError('Unsupported propagation type ', w.planetype, ' to ', plane)

        # Multiply by plane
        w = plane.multiply(w)

    # TODO: More advanced use cases may fail this assertion. We should figure out a
    # way to handle them here
    #assert w.depth == 1

    return w


# TODO:
# 1. accept flatten parameter
# 2. accept some sort of chip processing lambda function
# 3. accept data dtype parameter

def _propagate_pti(w, pixelscale, npix, oversample):

    shift = w.shift(pixelscale, oversample)
    w.tilt = []

    # Integer portion of the shift that will be accounted for
    # later
    fix_shift = np.fix(shift)
    w.center = list(fix_shift)

    # Residual subpixel portion of the shift that is passed to
    # the DFT
    res_shift = shift - fix_shift

    #chip = np.empty((npix[0]*oversample, npix[1]*oversample), dtype=np.complex128)

    alpha = (w.pixelscale * pixelscale) / (w.wavelength * w.focal_length * oversample)

    if shift.shape[0] != w.depth:
        raise ValueError('Dimension mismatch between tilt and wavefront depth')

    for d in range(w.depth):
        # data[d] = fourier.dft2(w.data[d], alpha, npix, res_shift[d])
        w.data[d] = fourier.dft2(w.data[d], alpha, (npix[0]*oversample, npix[1]*oversample),
                                 res_shift[d], offset=w.offset[d], unitary=True)



    return w


def _standardize_shape(shape, default=()):
    if shape is None:
        shape = default
    shape = np.asarray(shape)
    if shape.shape == ():
        shape = np.append(shape, shape)
    return shape


def _standardize_bandpass(vec, default=()):
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
