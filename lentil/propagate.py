import numpy as np

import lentil
import lentil.field
import lentil.helper


def propagate(planes, wave, weight=None, npix=None, npix_prop=None, oversample=2,
              rebin=True, fit_tilt=False, min_q=2):
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
        Shape of output plane.
    npix_prop : int or (2,) tuple of ints, optional
        Shape of propagation output plane. If None (default),
        ``npix_prop = npix``. If ``npix_prop != npix``, the propagation
        result is placed in the appropriate location in the output plane.
        npix_prop cannot be larger than npix.
    oversample : int, optional
        Number of times to oversample the output plane. Default is 2.
    rebin : bool, optional
        If ``True``, return the output plane in the sampling given by
        ``pixelscale``, binning down the output plane by a factor of
        ``oversample`` as needed. Note that this operation preserves power
        in the output plane. Default is ``True``.
    fit_tilt : bool, optional
        If True, tilt is removed from each Plane and bookkept separately. The
        effective tilt is included by directly placing the propagation result
        in the correct location in the output. Setting fit_tilt = True can help
        to eliminate periodic wraparound errors caused by large tilts aliasing
        in the complex phasor. Default is False. See :func:`lentil.Plane.fit_tilt`
        for more details.
    min_q : float, optional

    Returns
    -------
    psf : ndarray
        Resulting point spread function.

    """
    npix = lentil.helper.sanitize_shape(npix)
    npix_prop = npix if npix_prop is None else lentil.helper.sanitize_shape(npix_prop)
    wave = lentil.helper.sanitize_bandpass(wave)
    weight = np.ones_like(wave) if weight is None else lentil.helper.sanitize_bandpass(weight)

    # Create empty output
    out_shape = (npix[0] * oversample, npix[1]*oversample)

    # If the final plane is a Detector, we return intensity, otherwise we return
    # the complex field
    detector = True if isinstance(planes[-1], lentil.Detector) else False

    if detector:
        out = np.zeros(out_shape, dtype=float)
    else:
        out = np.zeros((len(wave), *out_shape), dtype=complex)

    if fit_tilt:
        planes = [plane.fit_tilt(inplace=False) for plane in planes]

    for n, (wl, wt) in enumerate(zip(wave, weight)):
        if wt > 0:
            w = lentil.Wavefront(wl)
            for plane, next_plane in _iterate_planes(planes):
                # propagate through the current plane
                w *= plane
                if w.planetype == 'pupil' and isinstance(next_plane, lentil.Image):
                    w = w.propagate_image(pixelscale=next_plane.pixelscale,
                                          npix=npix, npix_prop=npix_prop,
                                          oversample=oversample)
                else:
                    continue

            # place w
            if detector:
                for field in lentil.field.reduce(*w.data):
                    out = lentil.field.insert(field, out, intensity=True, weight=wt)
            else:
                for field in lentil.field.reduce(*w.data):
                    out[n] = lentil.field.insert(field, out[n], intensity=False,
                                                 weight=wt)
                out = np.squeeze(out)  # squeeze out singleton dimension

    if rebin:
        out = lentil.rebin(out, oversample)

    return out


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


