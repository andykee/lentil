import numpy as np

import lentil


def test_boundary():
    mask = lentil.util.hexagon((256, 256), 100, rotate=True)
    bounds = lentil.util.boundary(mask)
    assert np.array_equal(bounds, [42, 214, 28, 228])


def test_rebin():
    img = np.array([[1, 1, 2, 2], [1, 1, 2, 2], [3, 3, 4, 4], [3, 3, 4, 4]])
    factor = 2
    img_rebinned = lentil.util.rebin(img, factor)
    assert np.array_equal(img_rebinned, np.array([[4, 8], [12, 16]]))


def test_rebin_cube():
    img = np.zeros((2, 4, 4))
    img[0] = np.array([[1, 1, 2, 2], [1, 1, 2, 2], [3, 3, 4, 4], [3, 3, 4, 4]])
    img[1] = np.array([[2, 2, 4, 4], [2, 2, 4, 4], [6, 6, 8, 8], [6, 6, 8, 8]])
    factor = 2
    img_rebinned = lentil.util.rebin(img, factor)
    img_rebinned_expected = np.zeros((2, 2, 2))
    img_rebinned_expected[0] = np.array([[4, 8], [12, 16]])
    img_rebinned_expected[1] = np.array([[8, 16], [24, 32]])
    assert np.array_equal(img_rebinned, img_rebinned_expected)


def test_rescale_unitary():
    a = np.random.uniform(low=0, high=1, size=(100, 100))
    b = lentil.util.rescale(a, scale=0.5, order=3, mode='nearest', unitary=True)
    assert np.isclose(np.sum(a), np.sum(b))


def test_mesh_nonsquare():
    xx, yy = lentil.util.mesh((256, 512))
    assert xx.shape == (256, 512)


def test_pixelscle_nyquist():
    wave = np.random.uniform(low=350e-9, high=750e-9)
    f_number = np.random.uniform(low=5, high=30)
    assert lentil.util.pixelscale_nyquist(wave, f_number) == wave * f_number / 2


def test_make_mask():
    mask = lentil.util.circlemask((256, 256), 128)
    index = lentil.util.make_index(mask)
    assert np.array_equal(lentil.util.make_mask(index), mask)


def test_sparse():
    mask = lentil.util.circlemask((256, 256), 128)
    index = lentil.util.make_index(mask)
    vec = lentil.util.m2v(mask, index)
    mat = lentil.util.v2m(vec, index)
    assert np.array_equal(mat, mask)
