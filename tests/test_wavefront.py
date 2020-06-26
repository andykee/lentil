from lentil.wavefront import Wavefront


def test_wavefront_none_shape():
    w = Wavefront(650e-9)
    assert w.shape == (1, 1)
