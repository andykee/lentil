import pytest

import numpy as np
import scipy.special

@pytest.fixture
def airy2():

    def _airy2(diameter, focal_length, wavelength, pixelscale, shape, oversample):
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
    
    return _airy2