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
    phasor = p.amplitude * np.exp(1j*p.phase * 2 * np.pi / w.wavelength)

    assert np.array_equal(w1.data[0], phasor)
