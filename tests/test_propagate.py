import numpy as np
import scipy.interpolate
import scipy.special
import lentil


def airy(diameter, focal_length, wavelength, pixelscale, length, oversample=1):
    # https://en.wikipedia.org/wiki/Airy_disk#Mathematical_Formulation
    f_number = focal_length / diameter
    length *= oversample

    c = (length - 1) / 2
    x = np.arange(length, dtype=np.float)
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
    y, x = np.indices(shape, dtype=np.float)

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
    x = np.arange(length, dtype=np.float)
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
    def __init__(self, npix, coeffs=None):

        amplitude = lentil.modeltools.normalize_power(lentil.util.circle((npix, npix), npix//2))
        mask = lentil.util.circlemask((npix, npix), npix//2)

        if coeffs is None:
            coeffs = 5e-6 * np.random.uniform(low=-0.5, high=0.5, size=3)
            coeffs[0] = 0

        opd = lentil.zernike.zernike_compose(mask=mask, coeffs=coeffs)

        super().__init__(diameter=1,
                         focal_length=10,
                         pixelscale=1/npix,
                         amplitude=amplitude,
                         phase=opd,
                         mask=mask)

        self.coeffs = coeffs


class BasicDetector(lentil.Detector):
    def __init__(self, pixelscale=5e-6, shape=(512, 512)):
        super().__init__(pixelscale=pixelscale, shape=shape)


def test_amplitude_normalize_power():
    p = TiltPupil(npix=256)
    w = lentil.Wavefront(wavelength=650e-9)
    w *= p
    w.propagate_image(pixelscale=5e-6, npix=(64,64), oversample=1)
    psf = w.intensity
    assert (np.sum(psf) <= 1) and (np.sum(psf) >= 0.95)


def test_amplitude_normalize_power_oversample():
    p = TiltPupil(npix=256)
    w = lentil.Wavefront(wavelength=650e-9)
    w *= p
    w.propagate_image(pixelscale=5e-6, npix=(64,64), oversample=2)
    psf = w.intensity

    #planes = [TiltPupil(npix=256), BasicDetector()]
    #psf = lentil.propagate(planes, 650e-9, npix=(64, 64), oversample=2, rebin=True)
    assert (np.sum(psf) <= 1) and (np.sum(psf) >= 0.95)


def test_propagate_airy():
    # we set oversample = 1 and apply an "oversampling" factor of 25 to pixelscale
    # so that we can more easily enforce an odd npix. having an odd npix enables
    # comparison between numerical propagations via dft2() and the airy2() function
    # since their DC pixels will be coincident.
    p = TiltPupil(npix=512, coeffs=[0])

    psf_airy = airy2(diameter=p.diameter, focal_length=p.focal_length, wavelength=650e-9,
                     pixelscale=[5e-6, 5e-6], shape=(511, 511), oversample=1)
    psf_airy = psf_airy/np.max(psf_airy)

    w = lentil.Wavefront(wavelength=650e-9)
    w *= p
    w.propagate_image(pixelscale=5e-6, npix=511, oversample=1)
    psf = w.intensity
    psf = psf/np.max(psf)
    assert np.all(np.isclose(psf, psf_airy, atol=1e-3))


def test_propagate_tilt_angle():
    #planes = [TiltPupil(npix=256), BasicDetector()]

    p = TiltPupil(npix=256)
    
    w_phase = lentil.Wavefront(650e-9)
    w_phase *= p
    w_phase.propagate_image(pixelscale=5e-6, npix=128, oversample=2)
    psf_phase = w_phase.intensity

    p.fit_tilt()
    w_angle = lentil.Wavefront(650e-9)
    w_angle *= p
    w_angle.propagate_image(pixelscale=5e-6, npix=128, oversample=2)
    psf_angle = w_angle.intensity

    # threshold the PSFs so that the centroiding is consistent
    psf_phase[psf_phase < 1e-5] = 0
    psf_angle[psf_angle < 1e-5] = 0

    delta = np.abs(np.asarray(lentil.util.centroid(psf_phase)) - np.asarray(lentil.util.centroid(psf_angle)))
    assert np.all(delta <= 1e-10)


# def test_propagate_tilt_angle_poly():
#     planes = [TiltPupil(npix=256), BasicDetector()]
#     wave = np.array([600e-9, 650e-9, 700e-9])
#     psf_phase = lentil.propagate(planes, wave, npix=(128, 128), oversample=2, rebin=False, tilt='phase')
#     psf_angle = lentil.propagate(planes, wave, npix=(128, 128), oversample=2, rebin=False, tilt='angle')

#     # threshold the PSFs so that the centroiding is consistent
#     psf_phase[psf_phase < 1e-5] = 0
#     psf_angle[psf_angle < 1e-5] = 0

#     delta = np.abs(np.asarray(lentil.util.centroid(psf_phase)) - np.asarray(lentil.util.centroid(psf_angle)))
#     assert np.all(delta <= 1e-10)


# def test_shape_invariance_propagate_phase():
#     pupil_sm = TiltPupil(npix=256)
#     pupil_lg = TiltPupil(npix=512, coeffs=pupil_sm.coeffs)

#     planes_sm = [pupil_sm, BasicDetector()]
#     psf_sm = lentil.propagate(planes_sm, 650e-9, npix=(64, 64), oversample=2, rebin=False, tilt='phase')

#     planes_lg = [pupil_lg, BasicDetector()]
#     psf_lg = lentil.propagate(planes_lg, 650e-9, npix=(64, 64), oversample=2, rebin=False, tilt='phase')

#     delta = np.abs(np.asarray(lentil.util.centroid(psf_sm)) - np.asarray(lentil.util.centroid(psf_lg)))
#     assert np.all(delta <= 0.1)  # can't reliably do better than about 1/20th pixel so we''l score against 1/10th


# def test_shape_invariance_propagate_angle():
#     pupil_sm = TiltPupil(npix=256)
#     pupil_lg = TiltPupil(npix=512, coeffs=pupil_sm.coeffs)

#     planes_sm = [pupil_sm, BasicDetector()]
#     psf_sm = lentil.propagate(planes_sm, 650e-9, npix=(512, 512), oversample=2, rebin=True, npix_chip=64, tilt='angle')

#     planes_lg = [pupil_lg, BasicDetector()]
#     psf_lg = lentil.propagate(planes_lg, 650e-9, npix=(512, 512), oversample=2, rebin=True, npix_chip=64, tilt='angle')

#     delta = np.abs(np.asarray(lentil.util.centroid(psf_sm)) - np.asarray(lentil.util.centroid(psf_lg)))
#     assert np.all(delta <= 0.5)  # can't reliably do better than about 1/20th pixel so we''l score against 1/10th


# def test_wavelength_invariance_propagate_phase():
#     planes = [TiltPupil(npix=256), BasicDetector()]
#     psf_650 = lentil.propagate(planes, 650e-9, npix=(128, 128), oversample=2, rebin=False, tilt='phase')
#     psf_900 = lentil.propagate(planes, 900e-9, npix=(128, 128), oversample=2, rebin=False, tilt='angle')

#     psf_650 = psf_650/np.max(psf_650)
#     psf_900 = psf_900/np.max(psf_900)

#     psf_650[psf_650 < 1e-2] = 0
#     psf_900[psf_900 < 1e-2] = 0

#     delta = np.abs(np.asarray(lentil.util.centroid(psf_650)) - np.asarray(lentil.util.centroid(psf_900)))
#     assert np.all(delta <= 5e-2)  # can't reliably do better than about 1/50th pixel so we''l score against 1/20th


#def test_wavelength_invariance_propagate_angle():
#    planes = [TiltPupil(npix=256), BasicDetector()]
#    psf_650 = lentil.propagate(planes, 650e-9, npix=(128, 128), oversample=2, rebin=False, tilt='phase')
#    psf_900 = lentil.propagate(planes, 900e-9, npix=(128, 128), oversample=2, rebin=False, tilt='angle')
#
#    psf_650 = psf_650/np.max(psf_650)
#    psf_900 = psf_900/np.max(psf_900)
#
#    psf_650[psf_650 < 1e-2] = 0
#    psf_900[psf_900 < 1e-2] = 0
#
#    delta = np.abs(np.asarray(lentil.util.centroid(psf_650)) - np.asarray(lentil.util.centroid(psf_900)))
#    assert np.all(delta <= 5e-2)


def test_propagate_tilt_phase_analytic():
    oversample = 10
    pixelscale = 5e-6
    npix = np.array([64, 64])
    pupil = TiltPupil(npix=256)

    w = lentil.Wavefront(650e-9)
    w *= pupil
    w.propagate_image(pixelscale=pixelscale, npix=npix, oversample=oversample)
    psf = w.intensity

    psf = psf/np.max(psf)
    psf[psf < 0.2] = 0  # threshold for centroiding

    center = npix//2 * oversample
    centroid = np.asarray(lentil.util.centroid(psf))
    shift = centroid - center

    pupil_tilt = pupil.coeffs[1:3]

    # The factor of 4 gets you from RMS to PV
    analytic_shift = -((pupil_tilt/pupil.diameter)*pupil.focal_length/pixelscale) * oversample * 4

    assert np.all((np.abs(shift - analytic_shift)/oversample) < 0.2)


def test_propagate_tilt_angle_analytic():
    oversample = 10
    pixelscale = 5e-6
    npix = np.array([64, 64])
    pupil = TiltPupil(npix=256)
    detector = BasicDetector()

    pupil.fit_tilt()
    w = lentil.Wavefront(650e-9)
    w *= pupil
    w.propagate_image(pixelscale=pixelscale, npix=npix, oversample=oversample)
    psf = w.intensity

    #psf = lentil.propagate([pupil, detector], wave=650e-9, npix=npix,
    #                       oversample=oversample, rebin=False, tilt='angle')
    psf = psf/np.max(psf)
    psf[psf < 0.2] = 0  # threshold for centroiding

    center = npix//2 * oversample
    centroid = np.asarray(lentil.util.centroid(psf))
    shift = centroid - center

    pupil_tilt = pupil.coeffs[1:3]

    # The factor of 4 gets you from RMS to PV
    analytic_shift = -((pupil_tilt/pupil.diameter)*pupil.focal_length/pixelscale) * oversample * 4

    assert np.all((np.abs(shift - analytic_shift)/oversample) < 0.2)


#def test_propagate_tilt_angle_flatten():
#    planes = [TiltPupil(npix=256), BasicDetector()]
#    psf = lentil.propagate(planes, [650e-9, 700e-9], npix=(128, 128), rebin=False, tilt='angle', flatten=True)
#    psf_stack = lentil.propagate(planes, [650e-9, 700e-9], npix=(128, 128), rebin=False, tilt='angle', flatten=False)
#
#    assert np.array_equal(psf, np.sum(psf_stack, axis=0))


#def test_propagate_tilt_phase_flatten():
#    planes = [TiltPupil(npix=256), BasicDetector()]
#    psf = lentil.propagate(planes, [650e-9, 700e-9], npix=(128, 128), rebin=False, tilt='phase', flatten=True)
#    psf_stack = lentil.propagate(planes, [650e-9, 700e-9], npix=(128, 128), rebin=False, tilt='phase', flatten=False)
#
#    assert np.array_equal(psf, np.sum(psf_stack, axis=0))
