import numpy as np

import lentil

def test_expc():
    a = np.random.uniform(low=-1, high=1, size=(5,5))
    assert np.allclose(np.exp(1j*a), lentil.mathtools.expc(a))