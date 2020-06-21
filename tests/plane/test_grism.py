import numpy as np
import lentil


def test_grism_center():
    dispersion = [1, 650e-9]
    trace = [2, 0]
    grism = lentil.Grism(dispersion=dispersion, trace=trace)

    x0 = np.random.uniform(low=-5, high=5)
    y0 = np.random.uniform(low=-5, high=5)
    x, y = grism.shift(wavelength=650e-9, xs=x0, ys=y0)
    assert np.all((x == x0, y == y0))


def test_grism_shift():
    grism = lentil.Grism()
    grism.trace = [1, 1]
    grism.dispersion = [1, 650e-9]
    wave = 900e-9
    x, y = grism.shift(wavelength=wave)

    assert x == (wave - grism.dispersion[1])/np.sqrt(2)
    assert y == 1+x

