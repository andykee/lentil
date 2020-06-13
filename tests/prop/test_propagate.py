import numpy as np
import lentil

from tests.airy import airy2


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
    planes = [TiltPupil(npix=256), BasicDetector()]
    psf = lentil.propagate(planes, 650e-9, npix=(64, 64))
    assert (np.sum(psf) <= 1) and (np.sum(psf) >= 0.95)


def test_amplitude_normalize_power_oversample():
    planes = [TiltPupil(npix=256), BasicDetector()]
    psf = lentil.propagate(planes, 650e-9, npix=(64, 64), oversample=2, rebin=False)
    assert (np.sum(psf) <= 1) and (np.sum(psf) >= 0.95)


def test_amplitude_normalize_power_oversample_rebin():
    planes = [TiltPupil(npix=256), BasicDetector()]
    psf = lentil.propagate(planes, 650e-9, npix=(64, 64), oversample=2, rebin=False)
    assert (np.sum(psf) <= 1) and (np.sum(psf) >= 0.95)


def test_propagate_airy():
    # we set oversample = 1 and apply an "oversampling" factor of 25 to pixelscale
    # so that we can more easily enforce an odd npix. having an odd npix enables
    # comparison between numerical propagations via dft2() and the airy2() function
    # since their DC pixels will be coincident.
    p = TiltPupil(npix=512, coeffs=[0])
    d = BasicDetector(pixelscale=5e-6 / 25)

    airy = airy2(diameter=p.diameter, focal_length=p.focal_length, wavelength=650e-9,
                 pixelscale=d.pixelscale, shape=(511, 511), oversample=1)
    airy = airy/np.max(airy)

    planes = [p, d]
    psf = lentil.propagate(planes, wave=650e-9, npix=511, oversample=1)
    psf = psf/np.max(psf)

    assert np.all(np.isclose(psf, airy, atol=1e-3))


def test_propagate_tilt_angle_mono():
    planes = [TiltPupil(npix=256), BasicDetector()]
    psf_phase = lentil.propagate(planes, 650e-9, npix=(128, 128), oversample=2, rebin=False, tilt='phase')
    psf_angle = lentil.propagate(planes, 650e-9, npix=(128, 128), oversample=2, rebin=False, tilt='angle')

    # threshold the PSFs so that the centroiding is consistent
    psf_phase[psf_phase < 1e-5] = 0
    psf_angle[psf_angle < 1e-5] = 0

    delta = np.abs(np.asarray(lentil.util.centroid(psf_phase)) - np.asarray(lentil.util.centroid(psf_angle)))
    assert np.all(delta <= 1e-10)


def test_propagate_tilt_angle_poly():
    planes = [TiltPupil(npix=256), BasicDetector()]
    wave = np.array([600e-9, 650e-9, 700e-9])
    psf_phase = lentil.propagate(planes, wave, npix=(128, 128), oversample=2, rebin=False, tilt='phase')
    psf_angle = lentil.propagate(planes, wave, npix=(128, 128), oversample=2, rebin=False, tilt='angle')

    # threshold the PSFs so that the centroiding is consistent
    psf_phase[psf_phase < 1e-5] = 0
    psf_angle[psf_angle < 1e-5] = 0

    delta = np.abs(np.asarray(lentil.util.centroid(psf_phase)) - np.asarray(lentil.util.centroid(psf_angle)))
    assert np.all(delta <= 1e-10)


def test_shape_invariance_propagate_phase():
    pupil_sm = TiltPupil(npix=256)
    pupil_lg = TiltPupil(npix=512, coeffs=pupil_sm.coeffs)

    planes_sm = [pupil_sm, BasicDetector()]
    psf_sm = lentil.propagate(planes_sm, 650e-9, npix=(64, 64), oversample=2, rebin=False, tilt='phase')

    planes_lg = [pupil_lg, BasicDetector()]
    psf_lg = lentil.propagate(planes_lg, 650e-9, npix=(64, 64), oversample=2, rebin=False, tilt='phase')

    delta = np.abs(np.asarray(lentil.util.centroid(psf_sm)) - np.asarray(lentil.util.centroid(psf_lg)))
    assert np.all(delta <= 0.1)  # can't reliably do better than about 1/20th pixel so we''l score against 1/10th


def test_shape_invariance_propagate_angle():
    pupil_sm = TiltPupil(npix=256)
    pupil_lg = TiltPupil(npix=512, coeffs=pupil_sm.coeffs)

    planes_sm = [pupil_sm, BasicDetector()]
    psf_sm = lentil.propagate(planes_sm, 650e-9, npix=(512, 512), oversample=2, rebin=True, npix_chip=64, tilt='angle')

    planes_lg = [pupil_lg, BasicDetector()]
    psf_lg = lentil.propagate(planes_lg, 650e-9, npix=(512, 512), oversample=2, rebin=True, npix_chip=64, tilt='angle')

    delta = np.abs(np.asarray(lentil.util.centroid(psf_sm)) - np.asarray(lentil.util.centroid(psf_lg)))
    assert np.all(delta <= 0.5)  # can't reliably do better than about 1/20th pixel so we''l score against 1/10th


def test_wavelength_invariance_propagate_phase():
    planes = [TiltPupil(npix=256), BasicDetector()]
    psf_650 = lentil.propagate(planes, 650e-9, npix=(128, 128), oversample=2, rebin=False, tilt='phase')
    psf_900 = lentil.propagate(planes, 900e-9, npix=(128, 128), oversample=2, rebin=False, tilt='angle')

    psf_650 = psf_650/np.max(psf_650)
    psf_900 = psf_900/np.max(psf_900)

    psf_650[psf_650 < 1e-2] = 0
    psf_900[psf_900 < 1e-2] = 0

    delta = np.abs(np.asarray(lentil.util.centroid(psf_650)) - np.asarray(lentil.util.centroid(psf_900)))
    assert np.all(delta <= 5e-2)  # can't reliably do better than about 1/50th pixel so we''l score against 1/20th


def test_wavelength_invariance_propagate_angle():
    planes = [TiltPupil(npix=256), BasicDetector()]
    psf_650 = lentil.propagate(planes, 650e-9, npix=(128, 128), oversample=2, rebin=False, tilt='phase')
    psf_900 = lentil.propagate(planes, 900e-9, npix=(128, 128), oversample=2, rebin=False, tilt='angle')

    psf_650 = psf_650/np.max(psf_650)
    psf_900 = psf_900/np.max(psf_900)

    psf_650[psf_650 < 1e-2] = 0
    psf_900[psf_900 < 1e-2] = 0

    delta = np.abs(np.asarray(lentil.util.centroid(psf_650)) - np.asarray(lentil.util.centroid(psf_900)))
    assert np.all(delta <= 5e-2)


def test_propagate_tilt_phase_analytic():
    oversample = 10
    npix = np.array([64, 64])
    pupil = TiltPupil(npix=256)
    detector = BasicDetector()
    psf = lentil.propagate([pupil, detector], wave=650e-9, npix=npix,
                           oversample=oversample, rebin=False, tilt='phase')
    psf = psf/np.max(psf)
    psf[psf < 0.2] = 0  # threshold for centroiding

    center = npix//2 * oversample
    centroid = np.asarray(lentil.util.centroid(psf))
    shift = centroid - center

    pupil_tilt = pupil.coeffs[1:3]

    # The factor of 4 gets you from RMS to PV
    analytic_shift = -((pupil_tilt/pupil.diameter)*pupil.focal_length/detector.pixelscale) * oversample * 4

    assert np.all((np.abs(shift - analytic_shift)/oversample) < 0.2)


def test_propagate_tilt_angle_analytic():
    oversample = 10
    npix = np.array([64, 64])
    pupil = TiltPupil(npix=256)
    detector = BasicDetector()
    psf = lentil.propagate([pupil, detector], wave=650e-9, npix=npix,
                           oversample=oversample, rebin=False, tilt='angle')
    psf = psf/np.max(psf)
    psf[psf < 0.2] = 0  # threshold for centroiding

    center = npix//2 * oversample
    centroid = np.asarray(lentil.util.centroid(psf))
    shift = centroid - center

    pupil_tilt = pupil.coeffs[1:3]

    # The factor of 4 gets you from RMS to PV
    analytic_shift = -((pupil_tilt/pupil.diameter)*pupil.focal_length/detector.pixelscale) * oversample * 4

    assert np.all((np.abs(shift - analytic_shift)/oversample) < 0.2)


def test_propagate_tilt_angle_flatten():
    planes = [TiltPupil(npix=256), BasicDetector()]
    psf = lentil.propagate(planes, [650e-9, 700e-9], npix=(128, 128), rebin=False, tilt='angle', flatten=True)
    psf_stack = lentil.propagate(planes, [650e-9, 700e-9], npix=(128, 128), rebin=False, tilt='angle', flatten=False)

    assert np.array_equal(psf, np.sum(psf_stack, axis=0))


def test_propagate_tilt_phase_flatten():
    planes = [TiltPupil(npix=256), BasicDetector()]
    psf = lentil.propagate(planes, [650e-9, 700e-9], npix=(128, 128), rebin=False, tilt='phase', flatten=True)
    psf_stack = lentil.propagate(planes, [650e-9, 700e-9], npix=(128, 128), rebin=False, tilt='phase', flatten=False)

    assert np.array_equal(psf, np.sum(psf_stack, axis=0))
