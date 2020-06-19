import numpy as np

from lentil import util
from lentil import fourier
from lentil.plane import Pupil, Image, Detector, Plane
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

    Note
    ----
    This function is a thin wrapper around Monocle's :class:`~lentil.prop.Propagate`
    object. If you need to do anything fancy that isn't provided  by the standard
    :func:`~lentil.propagate` method, consider subclassing or extending
    :class:`~lentil.prop.Propagate`.

    """

    with Propagate(planes) as p:
        psf = p.propagate(wave, weight, npix, npix_chip, oversample, rebin, tilt,
                          flatten)

    return psf


class Propagate:
    """Compute a polychromatic point spread function using Fraunhofer
    diffraction.

    Parameters
    ----------
    planes : list_like
        list of :class:`~lentil.Plane` objects

    Returns
    -------
    :class:`~lentil.prop.Propagate`

    Example
    -------
    ::

        with Propagate(planes) as p:
            psf = p.propagate(wave, weight, npix)

    """
    def __init__(self, planes):
        self.planes = planes

    def __enter__(self):
        # Pre-process and set up any static data we may need during the
        # propagation. Examples of this include setting caches, finding
        # dispersive elements, or really anything else that requires iterating
        # over the list of planes (which would be slow to perform with each
        # monochromatic wavefront propagation.)

        # Loop over the planes to build caches and identify any dispersive
        # elements
        for plane in self.planes:
            plane.cache_propagate()

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Reset all the plane caches to their initial empty states
        for plane in self.planes:
            plane.clear_cache_propagate()

    def propagate(self, wave, weight=None, npix=None, npix_chip=None, oversample=2,
                  rebin=True, tilt='phase', flatten=True):
        """Compute a polychromatic point spread function using Fraunhofer
        diffraction.

        Parameters
        ----------
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

        # TODO: move all this setup into __enter__:
        #   * npix
        #   * npix_chip
        #   * wave
        #   * weight
        if npix is None:
            npix = self.planes[-1].shape
        npix = np.asarray(npix)
        if npix.shape == ():
            npix = np.append(npix, npix)

        if npix_chip is None:
            npix_chip = npix
        npix_chip = np.asarray(npix_chip)
        if npix_chip.shape == ():
            npix_chip = np.append(npix_chip, npix_chip)

        wave = np.asarray(wave)
        if wave.shape == ():
            wave = wave[np.newaxis, ...]

        if weight is None:
            weight = np.ones(wave.shape)
        weight = np.asarray(weight)
        if weight.shape == ():
            weight = weight[np.newaxis, ...]

        # Create an empty output
        oversample_shape = (npix[0]*oversample, npix[1]*oversample)
        if flatten:
            output_shape = oversample_shape
        else:
            output_shape = (len(wave), oversample_shape[0], oversample_shape[1])

        output = np.zeros(output_shape)

        # We also need a temporary array to place the chips and compute intensity
        _output = np.zeros(oversample_shape, dtype=np.complex128)

        for n, (wl, wt) in enumerate(zip(wave, weight)):
            if wt > 0:
                w = _input_wavefront(self.planes, wl)
                w = self._propagate_mono(w, npix_chip, oversample, tilt)

                for d in range(w.depth):

                    # The shift term is given in terms of (x,y) but we place the chip in
                    # terms of (r,c)
                    # TODO: I think this will break if the last plane isn't a Detector
                    shift = np.flip(w.pixel_shift[d], axis=0)

                    # Compute the chip location
                    canvas_slice, chip_slice = _chip_insertion_slices(oversample_shape,
                                                                      (w.data.shape[1], w.data.shape[2]),
                                                                      shift)

                    # Insert the complex result in the output
                    if canvas_slice:
                        _output[canvas_slice] += w.data[d, chip_slice[0], chip_slice[1]]

                # Compute intensity
                if flatten:
                    output += np.abs(_output).real**2 * wt
                else:
                    output[n] = np.abs(_output).real**2 * wt

                # Zero out the local output array
                _output[:] = 0

        if rebin:
            output = util.rebin(output, oversample)

        return output

    def _propagate_mono(self, w, npix, oversample, tilt):
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
        for plane, next_plane in _iterate_planes(self.planes):

            # Multiply by the current plane
            w = plane.multiply(w, tilt)

            # Now, we propagate to the next plane in the optical system

            # TODO: figure out how/when to accumulate tilt

            if (w.planetype == 'pupil') and isinstance(next_plane, Image):
                if isinstance(next_plane, Detector):
                    w = _propagate_pupil_detector(w, next_plane.pixelscale, npix, oversample)
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

        return w


def _propagate_pupil_detector(w, detector_pixelscale, npix, oversample):

    shift = w.shift(detector_pixelscale, oversample)
    w.tilt = []

    # Integer portion of the shift that will be accounted for
    # later
    w.pixel_shift = np.fix(shift)

    # Residual subpixel portion of the shift that is passed to
    # the DFT
    res_shift = shift - w.pixel_shift

    npix = npix * oversample

    alpha = (w.pixelscale * detector_pixelscale) / (w.wavelength * w.focal_length * oversample)

    data = np.zeros((w.depth, npix[0], npix[1]), dtype=np.complex128)

    assert shift.shape[0] == w.depth, \
        'Dimension mismatch between tilt and wavefront depth'

    for d in range(w.depth):
        data[d] = fourier.dft2(w.data[d], alpha, npix, res_shift[d])

    w.data = data

    return w


def _input_wavefront(planes, wave):

    # TODO: need to be robust against the fact that the first plane may be a
    # source or tilt only object and therefor insufficient to fully set up
    # the wavefront object
    #  * let shape be None. If shape is None, Wavefront.data should just be 1

    return Wavefront(wave,
                     pixelscale=planes[0].pixelscale,
                     shape=planes[0].shape,
                     planetype=None)


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
