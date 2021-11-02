import numpy as np

import lentil
import lentil.util
import lentil.helper


def test_mesh_nonsquare():
    xx, yy = lentil.helper.mesh((256, 512))
    assert xx.shape == (256, 512)


def test_boundary_slice():
    a = np.zeros((10, 10))
    r, c = np.floor(np.random.uniform(low=0, high=10, size=2)).astype(int)
    a[r:r+3, c:c+3] = 1

    rmin = np.max((r, 0))
    rmax = np.min((r+3, a.shape[0]))
    cmin = np.max((c, 0))
    cmax = np.min((c+3, a.shape[1]))

    slc = lentil.helper.boundary_slice(a)

    assert (slice(rmin, rmax), slice(cmin, cmax)) == slc


def test_slice_offset():
    shift = np.random.uniform(low=-50, high=50, size=2).astype(int)

    a = lentil.circlemask((256, 256), 256//4, shift=shift)
    slc = lentil.helper.boundary_slice(a)
    offset = lentil.helper.slice_offset(slc, shape=a.shape)

    assert np.all(offset == shift)
