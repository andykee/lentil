import pytest
import numpy as np

import lentil


def test_boundary():
    mask = lentil.hexagon((256, 256), 100, rotate=True)
    bounds = lentil.boundary(mask)
    assert np.array_equal(bounds, [28, 228, 41, 215])


@pytest.mark.parametrize('shape, oversample, offset, result', [
    ((2, 2), 1, 0, [[0, 1], [2, 3]]),
    ((2, 2), 2, 0, [3]),
    ((2, 2), 2, (-1, 0), [1]),
    ((2, 2), 2, (0, -1), [2]),
    ((2, 2), 2, (-1, -1), [0]),
    ((3, 3), 1, 0, [[0, 1, 2], [3, 4, 5], [6, 7, 8]]),
    ((3, 3), 3, 0, [4]),
    ((3, 3), 3, (-1, 0), [1]),
    ((3, 3), 3, (1, 0), [7]),
    ((3, 3), 3, (0, -1), [3]),
    ((3, 3), 3, (0, 1), [5]),
    ((3, 3), 3, (-1, -1), [0]),
    ((3, 3), 3, (-1, 1), [2]),
    ((3, 3), 3, (1, -1), [6]),
    ((3, 3), 3, (1, 1), [8]),
    ((4, 4), 1, 0, [[0, 1, 2, 3], [4, 5, 6, 7],
                    [8, 9, 10, 11], [12, 13, 14, 15]]),
    ((4, 4), 2, 0, [[5, 7], [13, 15]]),
    ])
def test_sample(shape, oversample, offset, result):
    x = np.arange(np.prod(shape)).reshape(shape)
    assert np.all(lentil.sample(x, oversample, offset) == result)


@pytest.mark.parametrize('shape, oversample, offset', [
    ((2, 2), 1, 1),
    ((2, 2), 1, -2),
    ((2, 2), 2, 1),
    ((3, 3), 3, -2),
    ((3, 3), 3, 2),
    ((4, 4), 2, 1),
    ((4, 4), 2, -2)
    ])
def test_sample_value_error(shape, oversample, offset):
    x = np.arange(np.prod(shape)).reshape(shape)
    with pytest.raises(ValueError):
        lentil.sample(x, oversample, offset)


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


def test_pixelscle_nyquist():
    wave = np.random.uniform(low=350e-9, high=750e-9)
    f_number = np.random.uniform(low=5, high=30)
    assert lentil.util.pixelscale_nyquist(wave, f_number) == wave * f_number / 2


amp = lentil.circle((512, 512), 256)


def test_amp_norm():
    # make sure we're not starting with an already normalized array!
    assert np.sum(np.abs(amp)**2) != 1.0


def test_normalize_power():
    amp_norm = lentil.normalize_power(amp)
    assert np.isclose(np.sum(np.abs(amp_norm)**2), 1.0)


def test_normalize_power_2():
    amp_norm = lentil.normalize_power(amp, power=2)
    assert np.isclose(np.sum(np.abs(amp_norm)**2), 2.0)


@pytest.mark.parametrize('mode, indexing, result', [
    ('absolute', 'ij', (44, 94)),
    ('absolute', 'xy', (94, 44)),
    ('center', 'ij', (-20, 30)),
    ('center', 'xy', (30, 20)),
])
def test_centroid(mode, indexing, result):
    # circle shifted 20 rows up and 30 columns right of center
    circ = lentil.circle((128, 128), 25, shift=(-20, 30), antialias=False)
    assert np.allclose(lentil.centroid(circ, mode=mode, indexing=indexing), result)


@pytest.mark.parametrize('n', [64, 65])
@pytest.mark.parametrize('indexing', ['ij', 'xy'])
def test_centroid_center_mode_centered(n, indexing):
    circ = lentil.circle((n, n), 20, antialias=False)
    assert np.allclose(lentil.centroid(circ, mode='center', indexing=indexing), (0, 0))
