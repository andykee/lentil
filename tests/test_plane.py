import math
import numpy as np

import lentil
import lentil.wavefront


size = (10, 10)


class RandomPlane(lentil.Plane):
    def __init__(self):
        super().__init__(pixelscale=1,
                         amplitude=np.random.uniform(size=size),
                         phase=np.random.uniform(size=size))


def test_default_plane():
    # Ensure that a default Plane creates an object that won't have any
    # impact on an optical system (a perfect optic with no phase change and
    # perfect optical and spectral transmission).

    p = lentil.Plane()
    assert p.pixelscale is None
    assert np.all(p.amplitude == 1)
    assert np.all(p.phase == 0)
    assert p.mask == p.amplitude


def test_wavefront_plane_multiply():
    p = RandomPlane()
    w = lentil.wavefront.Wavefront(650e-9)

    w1 = p.multiply(w)

    slc = lentil.util.boundary_slice(p.mask)
    phasor = p.amplitude[slc] * np.exp(-2*np.pi*1j*p.phase[slc]/w.wavelength)

    assert np.array_equal(w1.data[0].data, phasor)


class CircularPupil(lentil.Pupil):
    def __init__(self):
        super().__init__(focal_length=10,
                         pixelscale=2/256,
                         amplitude=lentil.util.circle((256, 256), 128),
                         phase=np.zeros((256, 256)))


def test_wavefront_pupil_multiply():
    p = CircularPupil()
    w = lentil.wavefront.Wavefront(650e-9)
    w = p.multiply(w)
    phasor = p.amplitude * np.exp(1j*p.phase * 2 * np.pi / w.wavelength)

    assert np.array_equal(w.data[0].data, phasor)
    assert w.focal_length == p.focal_length


def test_pupil_rescale_power():
    p = CircularPupil()
    pr = p.rescale(3, inplace=False)

    amp_power = np.sum(np.abs(p.amplitude)**2)
    ampr_power = np.sum(np.abs(pr.amplitude)**2)
    assert math.isclose(amp_power, ampr_power, rel_tol=1e-3)


def test_grism_center():
    dispersion = [1, 650e-9]
    trace = [2, 0]
    grism = lentil.Grism(dispersion=dispersion, trace=trace)

    x0 = np.random.uniform(low=-5, high=5)
    y0 = np.random.uniform(low=-5, high=5)
    x, y = grism.shift(wavelength=650e-9, xs=x0, ys=y0)
    assert np.all((x == x0, y == y0))


def test_grism_shift():
    grism = lentil.Grism(trace=[1,1], dispersion=[1,650e-9])
    wave = 900e-9
    x, y = grism.shift(wavelength=wave)

    assert x == (wave - grism.dispersion[1])/np.sqrt(2)
    assert y == 1+x

