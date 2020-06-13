import numpy as np
import lentil


def test_saturation_capacity():
    gain = lentil.detector.Gain(gain=1, saturation_capacity=1)
    img = 10 * np.ones((10, 10))
    assert np.array_equal(gain(img), np.ones((10, 10)))


def test_float_gain():
    g = np.random.uniform(size=1)
    gain = lentil.detector.Gain(g, saturation_capacity=None)
    img = 10*np.random.uniform(size=(10, 10))
    assert np.array_equal(gain(img), np.floor(g*img))


def test_array_gain():
    g = np.random.uniform(size=(10, 10))
    gain = lentil.detector.Gain(g, saturation_capacity=None)
    img = 10*np.random.uniform(size=(10, 10))
    assert np.array_equal(gain(img), np.floor(g*img))


def test_nonlinear_gain():
    g = np.random.uniform(low=0.9, high=1.5, size=4)
    gain = lentil.detector.PolynomialGain(g, saturation_capacity=None)
    img = 10*np.ones((10, 10))
    assert np.isclose(gain(img)[0, 0], np.polyval(g, 10))


def test_nonlinear_gain_cube():
    g = np.random.uniform(low=0.9, high=1.5, size=(10, 10, 4))
    gain = lentil.detector.PolynomialGain(g, saturation_capacity=None)
    img = 10*np.ones((10, 10))
    index = np.random.uniform(low=0, high=9, size=2)
    r = int(index[0])
    c = int(index[1])
    assert np.isclose(gain(img)[r, c], np.polyval(g[r, c, :], 10))

