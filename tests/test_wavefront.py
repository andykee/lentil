import numpy as np
import lentil


def test_default_wavefront():
    w = lentil.Wavefront(wavelength=500e-9)
    assert np.array_equal(w.field, 1+0j)
    assert np.array_equal(w.intensity, 1)
