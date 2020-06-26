import pytest
import numpy as np

from lentil.radiometry import Spectrum


def test_spectrum_wave_value_shape_match():
    wave = [1, 2, 3]
    value = [4, 5, 6, 7]
    with pytest.raises(ValueError):
        s = Spectrum(wave, value)


def test_spectrum_add_scalar():
    s = Spectrum([1, 2, 3], [1, 1, 1])
    r = s + 5
    assert np.array_equal(r.wave, s.wave)
    assert np.array_equal(r.value, [6, 6, 6])


def test_spectrum_add_spectrum():
    s1 = Spectrum([1, 2, 3], [1, 1, 1])
    s2 = Spectrum([1, 2, 3], [2, 3, 4])
    r = s1 + s2
    assert np.array_equal(r.wave, s1.wave)
    assert np.array_equal(r.value, [3, 4, 5])


def test_spectrum_add_array():
    s = Spectrum([1, 2, 3], [1, 1, 1])
    r = s + [5, 6, 7]
    assert np.array_equal(r.wave, s.wave)
    assert np.array_equal(r.value, [6, 7, 8])


def test_spectrum_subtract_scalar():
    s = Spectrum([1, 2, 3], [3, 3, 3])
    r = s - 2
    assert np.array_equal(r.wave, s.wave)
    assert np.array_equal(r.value, [1, 1, 1])


def test_spectrum_subtract_spectrum():
    s1 = Spectrum([1, 2, 3], [2, 3, 4])
    s2 = Spectrum([1, 2, 3], [1, 1, 1])
    r = s1 - s2
    assert np.array_equal(r.wave, s1.wave)
    assert np.array_equal(r.value, [1, 2, 3])


def test_spectrum_subtract_array():
    s = Spectrum([1, 2, 3], [5, 6, 7])
    r = s - [1, 1, 1]
    assert np.array_equal(r.wave, s.wave)
    assert np.array_equal(r.value, [4, 5, 6])


def test_spectrum_mul_scalar():
    s = Spectrum([1, 2, 3], [1, 2, 3])
    r = s * 5
    assert np.array_equal(r.wave, s.wave)
    assert np.array_equal(r.value, [5, 10, 15])


def test_spectrum_mul_spectrum():
    s1 = Spectrum([1, 2, 3], [2, 2, 2])
    s2 = Spectrum([1, 2, 3], [2, 3, 4])
    r = s1 * s2
    assert np.array_equal(r.wave, s1.wave)
    assert np.array_equal(r.value, [4, 6, 8])


def test_spectrum_mul_array():
    s = Spectrum([1, 2, 3], [2, 2, 2])
    r = s * [5, 6, 7]
    assert np.array_equal(r.wave, s.wave)
    assert np.array_equal(r.value, [10, 12, 14])


def test_spectrum_truediv_scalar():
    s = Spectrum([5, 10, 15], [5, 10, 15])
    r = s / 5
    assert np.array_equal(r.wave, s.wave)
    assert np.array_equal(r.value, [1, 2, 3])


def test_spectrum_truediv_spectrum():
    s1 = Spectrum([1, 2, 3], [4, 6, 8])
    s2 = Spectrum([1, 2, 3], [2, 3, 4])
    r = s1 / s2
    assert np.array_equal(r.wave, s1.wave)
    assert np.array_equal(r.value, [2, 2, 2])


def test_spectrum_truediv_array():
    s = Spectrum([1, 2, 3], [10, 12, 14])
    r = s / [5, 6, 7]
    assert np.array_equal(r.wave, s.wave)
    assert np.array_equal(r.value, [2, 2, 2])


def test_spectrum_rmul_scalar():
    s = Spectrum([1, 2, 3], [1, 2, 3])
    r = 5 * s
    assert np.array_equal(r.wave, s.wave)
    assert np.array_equal(r.value, [5, 10, 15])


def test_spectrum_pow():
    s = Spectrum([1, 2, 3], [2, 3, 4])
    r = s**2
    assert(np.array_equal(r.wave, s.wave))
    assert(np.array_equal(r.value, s.value**2))


def test_spectrum_sample():
    s = Spectrum([1, 2, 3], [2, 3, 4])
    r = s.sample([1.5, 2.5])
    assert np.array_equal(r, [2.5, 3.5])


def test_spectrum_resample():
    s = Spectrum([1, 2, 3], [2, 3, 4])
    s.resample([1.5, 2.5])
    assert np.array_equal(s.wave, [1.5, 2.5])
    assert np.array_equal(s.value, [2.5, 3.5])


def test_spectrum_asarray():
    s = Spectrum([1, 2, 3], [2, 3, 4])
    r = s.asarray()
    assert np.array_equal(r, np.array([[1, 2, 3], [2, 3, 4]]))


def test_spectrum_bin_preserve_power_true():
    s = Spectrum([1, 2, 3, 4, 5], [2, 2, 2, 2, 2])
    assert np.array_equal(s.bin([2, 3]), [0.5, 1.5])


def test_spectrum_bin_preserve_power_false():
    s = Spectrum([1, 2, 3, 4, 5], [2, 2, 2, 2, 2])
    assert np.array_equal(s.bin([2, 3], preserve_power=False), [1, 3])


def test_spectrum_integrate_no_bounds():
    s = Spectrum([1, 2, 3, 4, 5], [2, 2, 2, 2, 2])
    assert s.integrate() == 8.0


def test_spectrum_integrate_bounds():
    s = Spectrum([1, 2, 3, 4, 5], [2, 2, 2, 2, 2])
    assert s.integrate(2, 4) == 4.0


def test_spectrum_ends():
    s = Spectrum([1, 2, 3, 4, 5], [0.1, 1, 100, 1, 0.1])
    assert np.array_equal(s.ends(tol=0.001), [1, 3])


def test_spectrum_trim():
    s = Spectrum([1, 2, 3, 4, 5], [0.1, 1, 100, 1, 0.1])
    s.trim(tol=0.001)
    assert np.array_equal(s.wave, [2, 3, 4])
    assert np.array_equal(s.value, [1, 100, 1])


def test_spectrum_trim_value_all_zeros():
    s = Spectrum([1, 2, 3, 4, 5], [0, 0, 0, 0, 0])
    s.trim(tol=0.0001)
    assert np.array_equal(s.wave, [1, 2, 3, 4, 5])
    assert np.array_equal(s.value, [0, 0, 0, 0, 0])


def test_spectrum_crop():
    s = Spectrum([1, 2, 3, 4, 5], [0.1, 1, 100, 1, 0.1])
    s.crop(2, 5)
    assert np.array_equal(s.wave, [2, 3, 4, 5])
    assert np.array_equal(s.value, [1, 100, 1, 0.1])


def test_spectrum_crop_outside_wave_limits():
    s = Spectrum([1, 2, 3, 4, 5], [0.1, 1, 100, 1, 0.1])
    s.crop(2, 10)
    assert np.array_equal(s.wave, [2, 3, 4, 5])
    assert np.array_equal(s.value, [1, 100, 1, 0.1])
