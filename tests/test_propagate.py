import math

import lentil
import numpy as np
import pytest


def test_amplitude_normalize_power(pupil):
    p = pupil(focal_length=10, diameter=1, shape=(256,256), radius=120, coeffs=None)
    w = lentil.Wavefront(wavelength=650e-9)
    w *= p
    w = lentil.propagate_dft(w, shape=(64,64), pixelscale=5e-6, oversample=1)
    psf = w.intensity
    assert (np.sum(psf) <= 1) and (np.sum(psf) >= 0.95)


def test_amplitude_normalize_power_oversample(pupil):
    p = pupil(focal_length=10, diameter=1, shape=(256,256), radius=120, coeffs=None)
    w = lentil.Wavefront(wavelength=650e-9)
    w *= p
    w = lentil.propagate_dft(w, shape=(64,64), pixelscale=5e-6, oversample=2)
    psf = w.intensity
    assert (np.sum(psf) <= 1) and (np.sum(psf) >= 0.95)


def test_propagate_airy(airy2, pupil):
    # we set oversample = 1 and apply an "oversampling" factor of 25 to pixelscale
    # so that we can more easily enforce an odd npix. having an odd npix enables
    # comparison between numerical propagations via dft2() and the airy2() function
    # since their DC pixels will be coincident.
    p = pupil(focal_length=10, diameter=1, shape=(512,512), radius=250, coeffs=None)

    psf_airy = airy2(diameter=p.diameter, focal_length=p.focal_length, wavelength=650e-9,
                     pixelscale=[5e-6, 5e-6], shape=(511, 511), oversample=1)
    psf_airy = psf_airy/np.max(psf_airy)

    w = lentil.Wavefront(wavelength=650e-9)
    w *= p
    w = lentil.propagate_dft(w, shape=511, pixelscale=5e-6, oversample=1)
    psf = w.intensity
    psf = psf/np.max(psf)
    assert np.all(np.isclose(psf, psf_airy, atol=1e-3))


def test_propagate_fit_tilt(pupil):
    coeffs = 5e-6 * np.random.uniform(low=-0.5, high=0.5, size=3)
    coeffs[0] = 0
    p = pupil(focal_length=10, diameter=1, shape=(256,256), radius=120, coeffs=coeffs)

    w = lentil.Wavefront(650e-9)
    w = w * p
    w = lentil.propagate_dft(w, shape=128, pixelscale=5e-6, oversample=2)
    psf = w.intensity

    p.fit_tilt(inplace=True)
    w_fit_tilt = lentil.Wavefront(650e-9)
    w_fit_tilt = w_fit_tilt * p
    w_fit_tilt = lentil.propagate_dft(w_fit_tilt, shape=128, pixelscale=5e-6, oversample=2)
    psf_fit_tilt = w_fit_tilt.intensity

    # threshold the PSFs so that the centroiding is consistent
    psf[psf < 1e-5] = 0
    psf_fit_tilt[psf_fit_tilt < 1e-5] = 0

    delta = np.abs(np.asarray(lentil.util.centroid(psf)) - np.asarray(lentil.util.centroid(psf_fit_tilt)))
    assert np.all(delta <= 1e-10)


def test_propagate_tilt_analytic(pupil):
    oversample = 10
    pixelscale = 5e-6
    npix = np.array([64, 64])
    
    coeffs = 5e-6 * np.random.uniform(low=-0.5, high=0.5, size=3)
    coeffs[0] = 0
    p = pupil(focal_length=10, diameter=1, shape=(256,256), radius=120, coeffs=coeffs)

    w = lentil.Wavefront(650e-9)
    w = w * p
    w = lentil.propagate_dft(w, shape=npix, pixelscale=pixelscale, oversample=oversample)
    psf = w.intensity

    psf = psf/np.max(psf)
    psf[psf < 0.2] = 0  # threshold for centroiding

    shift = np.asarray(lentil.centroid(psf, mode='center', indexing='xy'))

    pupil_tilt = coeffs[1:3]

    # The factor of 4 gets you from RMS to PV
    analytic_shift = ((pupil_tilt/p.diameter)*p.focal_length/pixelscale) * oversample * 4


    assert np.all((np.abs(shift - analytic_shift)/oversample) < 0.2)


def test_propagate_tilt_fit_tilt_analytic(pupil):
    oversample = 10
    pixelscale = 5e-6
    npix = np.array([64, 64])
    
    coeffs = 5e-6 * np.random.uniform(low=-0.5, high=0.5, size=3)
    coeffs[0] = 0
    p = pupil(focal_length=10, diameter=1, shape=(256,256), radius=120, coeffs=coeffs)

    p.fit_tilt(inplace=True)
    w = lentil.Wavefront(650e-9)
    w = w * p
    w = lentil.propagate_dft(w, shape=npix, pixelscale=pixelscale, oversample=oversample)
    psf = w.intensity

    psf = psf/np.max(psf)
    psf[psf < 0.2] = 0  # threshold for centroiding

    shift = np.asarray(lentil.centroid(psf, mode='center', indexing='xy'))

    pupil_tilt = coeffs[1:3]

    # The factor of 4 gets you from RMS to PV
    analytic_shift = ((pupil_tilt/p.diameter)*p.focal_length/pixelscale) * oversample * 4

    assert np.all((np.abs(shift - analytic_shift)/oversample) < 0.2)


def test_propagate_no_output(pupil):

    p = pupil(focal_length=10, diameter=1, shape=(255,255), radius=120, coeffs=[0, 1e-3])
    p.fit_tilt(inplace=True)
    w = lentil.Wavefront(650e-9)
    w *= p
    w = lentil.propagate_dft(w, shape=(64,64), pixelscale=5e-6, oversample=2)
    assert np.all(w.intensity == 0)


def test_propagate_resample(pupil):

    coeffs = np.random.uniform(low=-200e-9, high=200e-9, size=6)

    p = pupil(focal_length=10, diameter=1, shape=(256,256), radius=120, coeffs=coeffs)
    w = lentil.Wavefront(650e-9)
    w *= p
    wi = lentil.propagate_dft(w, shape=(64,64), pixelscale=5e-6, oversample=10)

    p2 = p.rescale(scale=3)
    w2 = lentil.Wavefront(650e-9)
    w2 *= p2
    w2i = lentil.propagate_dft(w2, shape=(64,64), pixelscale=5e-6, oversample=10)

    # compute cross correlation between wi and w2i
    xc = np.fft.ifftshift(np.conj(np.fft.fft2(wi.intensity)) * np.fft.fft2(w2i.intensity))
    cent = lentil.centroid(np.abs(xc))

    assert np.allclose(cent, [320, 320])
    assert math.isclose(np.sum(wi.intensity), np.sum(w2i.intensity), rel_tol=1e-2)

@pytest.mark.parametrize('shift', [(5, 5), (-3, 7), (0, -4)])
def test_propagate_shift_integer(pupil, shift):
    oversample = 5
    npix = np.array([64, 64])
    p = pupil(focal_length=10, diameter=1, shape=(256, 256), radius=120)

    w = lentil.Wavefront(650e-9) * p
    psf = lentil.propagate_dft(w, shape=npix, pixelscale=5e-6,
                               oversample=oversample).intensity

    w = lentil.Wavefront(650e-9) * p
    psf_shift = lentil.propagate_dft(w, shape=npix, pixelscale=5e-6,
                                     oversample=oversample,
                                     shift=shift, indexing='ij').intensity

    # an integer shift should exactly translate the PSF. compare the region
    # where the shifted and unshifted PSFs overlap
    dr, dc = np.asarray(shift) * oversample
    nr, nc = npix * oversample
    shifted = psf_shift[max(dr, 0):nr+min(dr, 0), max(dc, 0):nc+min(dc, 0)]
    unshifted = psf[max(-dr, 0):nr+min(-dr, 0), max(-dc, 0):nc+min(-dc, 0)]

    assert np.allclose(shifted, unshifted, rtol=0, atol=1e-12*psf.max())


def test_propagate_shift_subpixel(pupil):
    oversample = 10
    npix = np.array([64, 64])
    shift = np.random.uniform(low=-5, high=5, size=2)
    p = pupil(focal_length=10, diameter=1, shape=(256, 256), radius=120)

    w = lentil.Wavefront(650e-9) * p
    w = lentil.propagate_dft(w, shape=npix, pixelscale=5e-6,
                             oversample=oversample, shift=shift)
    psf = w.intensity

    psf = psf/np.max(psf)
    psf[psf < 0.2] = 0  # threshold for centroiding

    center = npix//2 * oversample
    observed_shift = np.asarray(lentil.util.centroid(psf, indexing='xy')) - center

    assert np.all(np.abs(observed_shift - shift*oversample)/oversample < 0.02)


def test_propagate_shift_tilt_additive(pupil):
    oversample = 10
    npix = np.array([64, 64])
    shift = (3, -2)

    coeffs = 5e-6 * np.random.uniform(low=-0.5, high=0.5, size=3)
    coeffs[0] = 0
    p = pupil(focal_length=10, diameter=1, shape=(256, 256), radius=120,
              coeffs=coeffs)
    p.fit_tilt(inplace=True)

    p = pupil(focal_length=10, diameter=1, shape=(256, 256), radius=120,
              coeffs=coeffs)
    p.fit_tilt(inplace=True)

    w = lentil.Wavefront(650e-9) * p
    psf = lentil.propagate_dft(w, shape=npix, pixelscale=5e-6,
                               oversample=oversample).intensity

    w = lentil.Wavefront(650e-9) * p
    psf_shift = lentil.propagate_dft(w, shape=npix, pixelscale=5e-6,
                                     oversample=oversample,
                                     shift=shift).intensity

    psf = psf/np.max(psf)
    psf[psf < 0.2] = 0
    psf_shift = psf_shift/np.max(psf_shift)
    psf_shift[psf_shift < 0.2] = 0

    delta = (np.asarray(lentil.util.centroid(psf_shift, indexing='xy')) -
             np.asarray(lentil.util.centroid(psf, indexing='xy')))

    assert np.all(np.abs(delta - np.asarray(shift)*oversample)/oversample < 0.02)
