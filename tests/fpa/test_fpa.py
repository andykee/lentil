import numpy as np

import lentil as le


def test_pixelate_unitary():

    a = np.random.rand(100, 100)
    b = le.FPA.pixelate(a, oversample=1)
    assert np.isclose(np.sum(a), np.sum(b))
