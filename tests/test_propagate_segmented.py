import numpy as np
import lentil


class SimpleSegmentedPupil(lentil.Pupil):
    def __init__(self, npix=256):
        segmask1 = lentil.circle((npix, npix), 0.2 * npix, shift=(-0.25 * npix, 0.2 * npix), antialias=False)
        segmask2 = lentil.circle((npix, npix), 0.2 * npix, shift=(0.25 * npix, 0.2 * npix), antialias=False)
        segmask3 = lentil.circle((npix, npix), 0.2 * npix, shift=(0, -0.2 * npix), antialias=False)
        mask = np.stack((segmask1, segmask2, segmask3))
        global_mask = np.sum(mask, axis=0)

        zern = np.random.uniform(-0.5, 0.5, size=2) * 2e-6
        seg_zern = np.random.uniform(-0.5, 0.5, size=6) * 4e-6

        opdp = lentil.zernike_compose(global_mask, np.concatenate(([0], zern)))
        opd1 = lentil.zernike_compose(mask[0], np.concatenate(([0], seg_zern[0:2])))
        opd2 = lentil.zernike_compose(mask[1], np.concatenate(([0], seg_zern[2:4])))
        opd3 = lentil.zernike_compose(mask[2], np.concatenate(([0], seg_zern[4:6])))
        opd = opdp + opd1 + opd2 + opd3

        super().__init__(focal_length=10, pixelscale=1/npix,
                         amplitude=lentil.normalize_power(global_mask), opd=opd,
                         mask=mask)


def test_propagate_tilt_angle_mono():

    p = SimpleSegmentedPupil()

    w1 = lentil.Wavefront(wavelength=650e-9)
    w1 *= p
    w1 = lentil.propagate_dft(w1, shape=(128,128), pixelscale=5e-6)
    psf = w1.intensity

    p.fit_tilt(inplace=True)
    w2 = lentil.Wavefront(wavelength=650e-9)
    w2 *= p
    w2 = lentil.propagate_dft(w2, shape=(128,128), pixelscale=5e-6)
    psf_fit_tilt = w2.intensity

    # Normalize and threshold the PSFs so that the centroiding is consistent
    psf /= np.max(psf)
    psf[psf < 0.2] = 0

    psf_fit_tilt /= np.max(psf_fit_tilt)
    psf_fit_tilt[psf_fit_tilt < 0.2] = 0

    delta = np.abs(np.asarray(lentil.util.centroid(psf)) - np.asarray(lentil.util.centroid(psf_fit_tilt)))
    assert np.all(delta <= 0.2)
