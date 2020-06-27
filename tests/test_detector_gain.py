import numpy as np
import lentil


def test_saturation_capacity():
    img = 10 * np.ones((10, 10))
    out = lentil.detector.adc(img, gain=1, saturation_capacity=1)
    assert np.array_equal(out, np.ones((10, 10)))


def test_float_gain():
    g = np.random.uniform(size=1)
    img = 10*np.random.uniform(size=(10, 10))
    out = lentil.detector.adc(img, g, saturation_capacity=None)
    assert np.array_equal(out, np.floor(g*img))


def test_array_gain():
    g = np.random.uniform(size=(10, 10))
    img = 10*np.random.uniform(size=(10, 10))
    out = lentil.detector.adc(img, g, saturation_capacity=None)
    assert np.array_equal(out, np.floor(g*img))


def test_nonlinear_gain():
    g = np.random.uniform(low=0.9, high=1.5, size=4)
    img = 10*np.ones((10, 10))
    out = lentil.detector.adc(img, g, saturation_capacity=None)
    assert np.isclose(out[0, 0], np.floor(np.polyval(np.append(g, 0), 10)))


def test_nonlinear_gain_cube():
    g = np.random.uniform(low=0.9, high=1.5, size=(4, 10, 10))
    img = 10*np.ones((10, 10))
    out = lentil.detector.adc(img, g, saturation_capacity=None)
    index = np.random.uniform(low=0, high=9, size=2)
    r = int(index[0])
    c = int(index[1])
    assert np.isclose(out[r, c], np.floor(np.polyval(np.append(g[:, r, c], 0), 10)))

