import numpy as np

import lentil
import lentil.wavefront


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
