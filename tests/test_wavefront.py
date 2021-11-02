import pytest
import numpy as np
import lentil


def test_default_wavefront():
    w = lentil.Wavefront(wavelength=500e-9)
    assert np.array_equal(w.field, 1+0j)
    assert np.array_equal(w.intensity, 1)

def test_wavefront_rmul():
    w = lentil.Wavefront(wavelength=500e-9)
    p = lentil.Plane()
    assert p * w

def test_wavefront_propagate_image_non_pupil():
    w = lentil.Wavefront(wavelength=500e-9)
    with pytest.raises(ValueError):
        w.propagate_image(pixelscale=5e-6, npix=64)


@pytest.mark.parametrize('field_shape, field_shift, output_shape', [
    ((5,5), (-25,0), (10,10)),
    ((5,5), (25,0), (10,10)),
    ((5,5), (0,-25), (10,10)),
    ((5,5), (0,25), (10,10))
])
def test_overlap(field_shape, field_shift, output_shape):
    assert lentil.wavefront._overlap(field_shape, field_shift, output_shape) is False
