import lentil
import numpy as np


def test_propagate_airy(airy2, pupil):
    diameter = 1
    focal_length = 10
    wavelength = 500e-9
    pixelscale = (5e-6, 5e-6)
    shape = (511, 511)
    oversample = 1

    psf_airy = airy2(diameter, focal_length, wavelength, pixelscale, shape, oversample)
    psf_airy = psf_airy/np.max(psf_airy)

    p = pupil(focal_length=focal_length,
              diameter=diameter,
              shape=(512, 512),
              radius=256,
              coeffs=None)

    w = lentil.Wavefront(wavelength=500e-9)
    w *= p
    w = lentil.propagate_fft(w, shape=(511,511), pixelscale=5e-6, oversample=1)
    psf = w.intensity
    psf = psf/np.max(psf)

    np.all(np.isclose(psf, psf_airy, atol=1e-3))