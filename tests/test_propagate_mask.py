import numpy as np
import lentil

def test_propagate_mask(pupil):
    p = pupil(focal_length=10, diameter=1, shape=(256,256), radius=120, coeffs=None)

    mask = lentil.rectangle((256,256), 128,128, shift=(20, -30), antialias=False)

    w = lentil.Wavefront(650e-9)
    w = p.multiply(w)
    w = lentil.propagate_dft(w, shape=128, pixelscale=5e-6, oversample=2)
    psf = w.intensity

    w_mask = lentil.Wavefront(650e-9)
    w_mask = p.multiply(w_mask)
    w_mask = lentil.propagate_dft(w_mask, shape=128, pixelscale=5e-6, oversample=2, mask=mask)
    psf_mask = w_mask.intensity

    assert np.allclose(psf_mask, psf*mask)


def test_propagate_mask_tilt_analytic(pupil):
    p = pupil(focal_length=10, diameter=1, shape=(256,256), radius=120, coeffs=[0,1e-6, 2e-6])

    mask = lentil.rectangle((256,256), 64, 64, shift=(20, -30), antialias=False)

    w = lentil.Wavefront(650e-9)
    w = p.multiply(w)
    w = lentil.propagate_dft(w, shape=128, pixelscale=5e-6, oversample=2)
    psf = w.intensity

    w_mask = lentil.Wavefront(650e-9)
    w_mask = p.multiply(w_mask)
    w_mask = lentil.propagate_dft(w_mask, shape=128, pixelscale=5e-6, oversample=2, mask=mask)
    psf_mask = w_mask.intensity

    assert np.allclose(psf_mask, psf*mask)

