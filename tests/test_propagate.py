import math
import numpy as np
import scipy.interpolate
import scipy.special
import lentil


def airy(diameter, focal_length, wavelength, pixelscale, length, oversample=1):
    # https://en.wikipedia.org/wiki/Airy_disk#Mathematical_Formulation
    f_number = focal_length / diameter
    length *= oversample

    c = (length - 1) / 2
    x = np.arange(length, dtype=float)
    x -= c

    q = x * (pixelscale / oversample)
    X = (np.pi * q) / (wavelength * f_number)

    # if length is odd, the center value will be zero which will throw a
    # divide by zero error. To avoid this, we'll set any zeros to machine
    # epsilon (eps)
    X[X == 0] = np.finfo(X.dtype).eps

    return (2 * scipy.special.jn(1, X) / X) ** 2


def airy2(diameter, focal_length, wavelength, pixelscale, shape, oversample=1):
    # https://en.wikipedia.org/wiki/Airy_disk#Mathematical_Formulation

    f_number = focal_length / diameter

    shape = np.asarray(shape)
    shape *= oversample

    c = (shape - 1) / 2
    y, x = np.indices(shape, dtype=float)

    x -= c[1]
    y -= c[0]

    x *= (pixelscale[0] / oversample)
    y *= (pixelscale[1] / oversample)

    q = np.sqrt(x ** 2 + y ** 2)
    X = (np.pi * q) / (wavelength * f_number)

    # if length is odd, the center value will be zero which will throw a
    # divide by zero error. To avoid this, we'll set any zeros to machine
    # epsilon (eps)
    X[X == 0] = np.finfo(X.dtype).eps

    return (2 * scipy.special.jn(1, X) / X) ** 2


def test_airy():
    diameter = 1
    focal_length = 10
    f_number = focal_length / diameter
    wave = 650e-9
    pixelscale = 5e-6 / 1000  # divide by 1000 here to fake oversampling
    length = 2 ** 14  # an appropriately large number
    oversample = 1

    airy1 = airy(diameter, focal_length, wave, pixelscale, length, oversample)

    c = (length - 1) / 2
    x = np.arange(length, dtype=float)
    x -= c
    x *= pixelscale

    airy_interp = scipy.interpolate.interp1d(x, airy1)

    # TODO: move this in to the documentation and just point to that as a reference here
    # derivation of airy FWHM expression
    # ----------------------------------
    # the half-maximum of the central airy disk occurs at x_hwhm ~= 1.61633 [1]
    # since x = (pi*q)/(wave*f_number), q = (x * wave * f_number)/pi where
    # q is the radial distance from the optics axis in the observation
    # (or focal) plane
    #
    # substituting in x_hwhm,
    # q_hwhm = (1.61633/pi) * wave * f_number
    #
    # because we are more interested in FWHM, we multiply by 2
    # q_fwhm = 2*(1.61633/pi) * wave * f_number = 1.028987 * wave * f_number
    #
    # Refs:
    # [1] https://en.wikipedia.org/wiki/Airy_disk#Mathematical_Formulation

    # verify the FWHM occurs at 1.029 * wave * f_number
    j1_hwhm = 1.61633
    fwhm = 2 * (j1_hwhm / np.pi) * wave * f_number  # actual linear distance
    assert (np.abs(airy_interp(fwhm / 2) - 0.5) < 1e-5)

    # verify the first five zeros (as listed on Wikipedia) occur in the
    # right places
    j1_zeros = np.asarray([3.8317, 7.0156, 10.1735, 13.3237, 16.4706])
    for j1_zero in j1_zeros:
        assert (airy_interp(j1_zero / np.pi * wave * f_number) < 1e-7)


def test_airy2():
    # evaluate the 2D Airy code against the 1D Airy code. since we've done
    # some pretty exhaustive testing of the 1D Airy code, this simple
    # comparison should be sufficient

    diameter = 1
    focal_length = 10
    wave = 650e-9
    pixelscale = 5e-6 / 10  # divide by 10 here to fake oversampling
    npix = 511
    oversample = 1
    slice_index = npix // 2

    a1 = airy(diameter, focal_length, wave, pixelscale, npix, oversample)
    a2 = airy2(diameter, focal_length, wave, (pixelscale, pixelscale), (npix, npix), oversample)

    res = a2[:, slice_index] - a1
    assert (np.all(np.abs(res) < 1e-9))


class TiltPupil(lentil.Pupil):
    def __init__(self, npix, diameter=1, coeffs=None):

        amplitude = lentil.normalize_power(lentil.util.circle((npix, npix), npix//2))
        mask = lentil.util.circlemask((npix, npix), npix//2)

        if coeffs is None:
            coeffs = 5e-6 * np.random.uniform(low=-0.5, high=0.5, size=3)
            coeffs[0] = 0

        opd = lentil.zernike_compose(mask=mask, coeffs=coeffs)

        super().__init__(focal_length=10,
                         pixelscale=1/npix,
                         amplitude=amplitude,
                         phase=opd,
                         mask=mask)

        self.diameter = diameter
        self.coeffs = coeffs


def test_amplitude_normalize_power():
    p = TiltPupil(npix=256)
    w = lentil.Wavefront(wavelength=650e-9)
    w *= p
    w = w.propagate_image(pixelscale=5e-6, npix=(64,64), oversample=1)
    psf = w.intensity
    assert (np.sum(psf) <= 1) and (np.sum(psf) >= 0.95)


def test_amplitude_normalize_power_oversample():
    p = TiltPupil(npix=256)
    w = lentil.Wavefront(wavelength=650e-9)
    w *= p
    w = w.propagate_image(pixelscale=5e-6, npix=(64,64), oversample=2)
    psf = w.intensity

    #planes = [TiltPupil(npix=256), BasicDetector()]
    #psf = lentil.propagate(planes, 650e-9, npix=(64, 64), oversample=2, rebin=True)
    assert (np.sum(psf) <= 1) and (np.sum(psf) >= 0.95)


def test_propagate_airy():
    # we set oversample = 1 and apply an "oversampling" factor of 25 to pixelscale
    # so that we can more easily enforce an odd npix. having an odd npix enables
    # comparison between numerical propagations via dft2() and the airy2() function
    # since their DC pixels will be coincident.
    p = TiltPupil(npix=512, diameter=1, coeffs=[0])

    psf_airy = airy2(diameter=p.diameter, focal_length=p.focal_length, wavelength=650e-9,
                     pixelscale=[5e-6, 5e-6], shape=(511, 511), oversample=1)
    psf_airy = psf_airy/np.max(psf_airy)

    w = lentil.Wavefront(wavelength=650e-9)
    w *= p
    w = w.propagate_image(pixelscale=5e-6, npix=511, oversample=1)
    psf = w.intensity
    psf = psf/np.max(psf)
    assert np.all(np.isclose(psf, psf_airy, atol=1e-3))


def test_propagate_tilt_angle():
    #planes = [TiltPupil(npix=256), BasicDetector()]

    p = TiltPupil(npix=256)

    w_phase = lentil.Wavefront(650e-9)
    w_phase = p.multiply(w_phase)
    w_phase = w_phase.propagate_image(pixelscale=5e-6, npix=128, oversample=2)
    psf_phase = w_phase.intensity

    p.fit_tilt()
    w_angle = lentil.Wavefront(650e-9)
    w_angle = w_angle * p
    w_angle = w_angle.propagate_image(pixelscale=5e-6, npix=128, oversample=2)
    psf_angle = w_angle.intensity

    # threshold the PSFs so that the centroiding is consistent
    psf_phase[psf_phase < 1e-5] = 0
    psf_angle[psf_angle < 1e-5] = 0

    delta = np.abs(np.asarray(lentil.util.centroid(psf_phase)) - np.asarray(lentil.util.centroid(psf_angle)))
    assert np.all(delta <= 1e-10)


def test_propagate_tilt_phase_analytic():
    oversample = 10
    pixelscale = 5e-6
    npix = np.array([64, 64])
    pupil = TiltPupil(npix=256, diameter=1)

    w = lentil.Wavefront(650e-9)
    w = pupil.multiply(w)
    w = w.propagate_image(pixelscale=pixelscale, npix=npix, oversample=oversample)
    psf = w.intensity

    psf = psf/np.max(psf)
    psf[psf < 0.2] = 0  # threshold for centroiding

    center = npix//2 * oversample
    centroid = np.asarray(lentil.util.centroid(psf))
    shift = centroid - center

    pupil_tilt = pupil.coeffs[1:3]

    # The factor of 4 gets you from RMS to PV
    analytic_shift = ((pupil_tilt/pupil.diameter)*pupil.focal_length/pixelscale) * oversample * 4

    # Analytic shift is in terms of (x, y) but we really want (r, c) to match the centroid
    analytic_shift = -analytic_shift[1], analytic_shift[0]

    assert np.all((np.abs(shift - analytic_shift)/oversample) < 0.2)


def test_propagate_tilt_angle_analytic():
    oversample = 10
    pixelscale = 5e-6
    npix = np.array([64, 64])
    pupil = TiltPupil(npix=256, diameter=1)

    pupil.fit_tilt()
    w = lentil.Wavefront(650e-9)
    w = w * pupil
    w = w.propagate_image(pixelscale=pixelscale, npix=npix, oversample=oversample)
    psf = w.intensity

    psf = psf/np.max(psf)
    psf[psf < 0.2] = 0  # threshold for centroiding

    center = npix//2 * oversample
    centroid = np.asarray(lentil.util.centroid(psf))
    shift = centroid - center

    pupil_tilt = pupil.coeffs[1:3]

    # The factor of 4 gets you from RMS to PV
    analytic_shift = ((pupil_tilt/pupil.diameter)*pupil.focal_length/pixelscale) * oversample * 4

    # Analytic shift is in terms of (x, y) but we really want (r, c) to match the centroid
    analytic_shift = -analytic_shift[1], analytic_shift[0]

    assert np.all((np.abs(shift - analytic_shift)/oversample) < 0.2)


def test_propagate_resample():

    amp = lentil.circle((256, 256), 120)
    coeffs = np.random.uniform(low=-200e-9, high=200e-9, size=6)
    opd = lentil.zernike_compose(amp, coeffs)

    p = lentil.Pupil(focal_length=10, pixelscale=1 / 240, amplitude=amp, phase=opd)
    w = lentil.Wavefront(650e-9)
    w *= p
    wi = w.propagate_image(pixelscale=5e-6, npix=64, oversample=10)

    p2 = p.rescale(scale=3, inplace=False)
    w2 = lentil.Wavefront(650e-9)
    w2 *= p2
    w2i = w2.propagate_image(pixelscale=5e-6, npix=64, oversample=10)

    # compute cross correlation between wi and w2i
    xc = np.fft.ifftshift(np.conj(np.fft.fft2(wi.intensity)) * np.fft.fft2(w2i.intensity))
    cent = lentil.centroid(np.abs(xc))

    assert np.allclose(cent, [320, 320])
    assert math.isclose(np.sum(wi.intensity), np.sum(w2i.intensity), rel_tol=1e-2)

