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
    assert p.segmask is None
    assert p.nseg == 0
    assert p.cache_attrs == ['amplitude', 'phase']


def test_set_plane_attrs_none():
    p = lentil.Plane(pixelscale=1, amplitude=1, phase=0, mask=1, segmask=1)
    p.pixelscale = None
    p.amplitude = None
    p.phase = None
    p.mask = None
    p.segmask = None

    assert p.pixelscale is None
    assert p.amplitude is None
    assert p.phase is None
    assert p.mask is None
    assert p.segmask is None
    assert p.shape is None
    assert p.ptt_vector is None


def test_wavefront_plane_multiply():
    p = RandomPlane()
    w = lentil.wavefront.Wavefront(650e-9, p.shape)

    w1 = p.multiply(w)
    phasor = p.amplitude * lentil.mathtools.expc(p.phase * 2 * np.pi / w.wavelength)

    assert np.array_equal(w1.data[0], phasor)


class CircularPupil(lentil.Pupil):
    def __init__(self):
        super().__init__(diameter=2,
                         focal_length=10,
                         pixelscale=2/256,
                         amplitude=lentil.util.circle((256, 256), 128),
                         phase=np.zeros((256, 256)))


def test_wavefront_pupil_multiply():
    p = CircularPupil()
    w = lentil.wavefront.Wavefront(650e-9, shape=p.shape)
    w = p.multiply(w)
    phasor = p.amplitude * np.exp(1j*p.phase * 2 * np.pi / w.wavelength)

    assert np.array_equal(w.data[0], phasor)
    assert w.planetype == 'pupil'
    assert w.focal_length == p.focal_length


def test_f_number():
    p = lentil.Pupil()
    p.diameter = np.random.normal(loc=1)
    p.focal_length = np.random.normal(loc=10)
    assert p.f_number == p.focal_length/p.diameter


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

