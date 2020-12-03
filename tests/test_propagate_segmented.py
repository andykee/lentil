import numpy as np
import lentil


class SimpleSegmentedPupil(lentil.Pupil):

    def __init__(self, npix=256):
        segmask1 = lentil.util.circlemask((npix, npix), 0.2 * npix, shift=(-0.25 * npix, 0.2 * npix))
        segmask2 = lentil.util.circlemask((npix, npix), 0.2 * npix, shift=(0.25 * npix, 0.2 * npix))
        segmask3 = lentil.util.circlemask((npix, npix), 0.2 * npix, shift=(0, -0.2 * npix))
        segmask = np.stack((segmask1, segmask2, segmask3))
        mask = np.sum(segmask, axis=0)

        zern = np.random.uniform(-0.5, 0.5, size=2) * 2e-6
        seg_zern = np.random.uniform(-0.5, 0.5, size=6) * 4e-6

        opdp = lentil.zernike.zernike_compose(mask, np.concatenate(([0], zern)))
        opd1 = lentil.zernike.zernike_compose(segmask[0], np.concatenate(([0], seg_zern[0:2])))
        opd2 = lentil.zernike.zernike_compose(segmask[1], np.concatenate(([0], seg_zern[2:4])))
        opd3 = lentil.zernike.zernike_compose(segmask[2], np.concatenate(([0], seg_zern[4:6])))
        opd = opdp + opd1 + opd2 + opd3

        super().__init__(diameter=1, focal_length=10, pixelscale=1/npix,
                         amplitude=lentil.modeltools.normalize_power(mask), phase=opd,
                         mask=mask, segmask=segmask)


class SimpleDetector(lentil.Detector):
    def __init__(self, pixelscale=5e-6, shape=(256, 256)):
        super().__init__(pixelscale, shape)


def test_propagate_tilt_angle_mono():
    planes = [SimpleSegmentedPupil(npix=256), SimpleDetector()]
    psf_phase = lentil.propagate(planes, 650e-9, npix=(128, 128), oversample=2, rebin=False, tilt='phase')
    psf_angle = lentil.propagate(planes, 650e-9, npix=(128, 128), oversample=2, rebin=False, tilt='angle')

    # Normalize and threshold the PSFs so that the centroiding is consistent
    psf_phase /= np.max(psf_phase)
    psf_phase[psf_phase < 0.2] = 0

    psf_angle /= np.max(psf_angle)
    psf_angle[psf_angle < 0.2] = 0

    delta = np.abs(np.asarray(lentil.util.centroid(psf_phase)) - np.asarray(lentil.util.centroid(psf_angle)))
    assert np.all(delta <= 0.2)
