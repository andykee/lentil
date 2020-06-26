import pytest
import numpy as np

from lentil.detector import Windowable

data = np.diag(np.arange(1, 11))
for i in np.diag(data):
    data[i-1, i-1:] = data[i-1, i-1]
    data[i-1:, i-1] = data[i-1, i-1]

even = Windowable(data)
odd = Windowable(data[:-1, :-1])


def test_even_window_shape():
    assert np.array_equal(even.window(shape=(2, 2)), np.array([[5, 5], [5, 6]]))


def test_even_window_window():
    assert np.array_equal(even.window(window=(0, 2, 0, 2)), np.array([[1, 1], [1, 2]]))


def test_odd_window_shape():
    assert np.array_equal(odd.window(shape=(2, 2)), np.array([[4, 4], [4, 5]]))


def test_odd_window_window():
    assert np.array_equal(odd.window(window=(0, 2, 0, 2)), np.array([[1, 1], [1, 2]]))


def test_window_shape_window():
    assert np.array_equal(even.window(shape=(2, 2), window=(8, 10, 8, 10)), np.array([[9, 9], [9, 10]]))


def test_window_shape_mismatch():
    with pytest.raises(AssertionError):
        even.window(shape=(3, 3), window=(0, 2, 0, 2))
