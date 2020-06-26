import numpy as np
import lentil


def test_jitter_scale():
    a = np.zeros((11, 11))
    a[5, 5] = 1

    j = lentil.convolvable.Jitter()
    scale = 2
    b = j(a, scale)

    def fwhm(arr):
        arr = arr / np.max(arr)
        arr[arr < 0.5] = 0
        extents = lentil.util.boundary(arr)
        return np.max([extents[1] - extents[0], extents[3] - extents[2]])

    assert fwhm(b) // 2 == scale


def test_jitter_unitary():
    a = np.zeros((11, 11))
    a[5, 5] = 1

    j = lentil.convolvable.Jitter()
    scale = 2
    b = j(a, scale)
    assert np.isclose(np.sum(a), np.sum(b))


def test_smear_unitary():
    a = np.zeros((11, 11))
    a[5, 5] = 1

    s = lentil.convolvable.Smear()
    distance = 2
    b = s(a, distance)
    assert np.isclose(np.sum(a), np.sum(b))


def test_pixel_unitary():
    a = np.random.uniform(low=0, high=1, size=(100, 100))
    pixel_sampling = lentil.convolvable.Pixel()
    b = pixel_sampling(a, oversample=1)
    assert np.isclose(np.sum(a), np.sum(b))
