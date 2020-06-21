import numpy as np
import lentil


def test_jitter():
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
