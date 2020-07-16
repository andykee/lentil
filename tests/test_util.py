import pytest
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


def test_pad():
    img = np.ones((3, 3))
    out = lentil.util.pad(img, (5, 5))
    truth = np.zeros((5, 5))
    truth[1:4, 1:4] = 1

    assert np.array_equal(out, truth)


def test_pad_cube():
    img = np.ones((2, 3, 3))
    out = lentil.util.pad(img, (5, 5))
    truth = np.zeros((2, 5, 5))
    truth[:, 1:4, 1:4] = 1

    assert np.array_equal(out, truth)


def test_pad_shrink():
    img = np.random.uniform(size=(5, 5))
    out = lentil.util.pad(img, (3, 3))
    assert np.array_equal(out, img[1:4, 1:4])


win_data = np.diag(np.arange(1, 11))
for i in np.diag(win_data):
    win_data[i-1, i-1:] = win_data[i-1, i-1]
    win_data[i-1:, i-1] = win_data[i-1, i-1]


def test_window_none():
    assert np.array_equal(lentil.util.window(win_data, shape=None, slice=None), win_data)


def test_window_single():
    assert np.array_equal(lentil.util.window(5), np.array(5))


def test_even_window_shape():
    assert np.array_equal(lentil.util.window(win_data, shape=(2, 2)), np.array([[5, 5], [5, 6]]))


def test_even_window_window():
    assert np.array_equal(lentil.util.window(win_data, slice=(0, 2, 0, 2)), np.array([[1, 1], [1, 2]]))


def test_odd_window_shape():
    assert np.array_equal(lentil.util.window(win_data[:-1, :-1], shape=(2, 2)), np.array([[4, 4], [4, 5]]))


def test_odd_window_window():
    assert np.array_equal(lentil.util.window(win_data[:-1, :-1], slice=(0, 2, 0, 2)), np.array([[1, 1], [1, 2]]))


def test_window_shape_window():
    assert np.array_equal(lentil.util.window(win_data, shape=(2, 2), slice=(8, 10, 8, 10)), np.array([[9, 9], [9, 10]]))


def test_window_shape_mismatch():
    with pytest.raises(AssertionError):
        lentil.util.window(win_data, shape=(3, 3), slice=(0, 2, 0, 2))


def test_mesh_nonsquare():
    xx, yy = lentil.util.mesh((256, 512))
    assert xx.shape == (256, 512)


def test_pixelscle_nyquist():
    wave = np.random.uniform(low=350e-9, high=750e-9)
    f_number = np.random.uniform(low=5, high=30)
    assert lentil.util.pixelscale_nyquist(wave, f_number) == wave * f_number / 2


def test_slit():
    slit = lentil.util.slit((5, 5), 1)
    assert np.array_equal(slit[:, 0], np.array([0,0,1,0,0]))


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
