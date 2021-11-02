import numpy as np

import lentil
import lentil.util
import lentil.helper


def test_slice_offset():
    shift = np.random.uniform(low=-50, high=50, size=2).astype(int)

    a = lentil.util.circlemask((256, 256), 256//4, shift=shift)
    slc = lentil.util.boundary_slice(a)
    offset = lentil.helper.slice_offset(slc, shape=a.shape)

    assert np.all(offset == shift)
