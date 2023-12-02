import pytest
import numpy as np
import lentil


def test_collect_charge_scalar():
    img = np.random.uniform(low=0, high=100, size=(10, 10))
    qe = np.random.uniform()
    out = lentil.detector.collect_charge(img, 1, qe)
    assert np.array_equal(out, img*qe)


def test_collect_charge_array():
    img = np.random.uniform(low=0, high=100, size=(5, 10, 10))
    qe = np.random.uniform(size=5)
    out = lentil.detector.collect_charge(img, np.ones(5), qe)
    assert np.array_equal(out, np.einsum('ijk,i->jk', img, qe))


def test_collect_charge_spectrum():
    img = np.random.uniform(low=0, high=100, size=(5, 10, 10))
    wave = np.arange(1, 6)
    value = np.random.uniform(size=5)
    qe = lentil.radiometry.Spectrum(wave, value)
    out = lentil.detector.collect_charge(img, wave, qe)
    assert np.array_equal(out, np.einsum('ijk,i->jk', img, value))


def test_collect_charge_bayer_even():
    img = np.random.uniform(low=0, high=100, size=(5, 2, 2))
    qe_red = np.random.uniform(size=5)
    qe_green = np.random.uniform(size=5)
    qe_blue = np.random.uniform(size=5)
    bayer_pattern = 'RGGB'
    out = lentil.detector.collect_charge_bayer(img, np.ones(5), qe_red,
                                               qe_green, qe_blue, bayer_pattern)

    ul = np.sum(img[:, 0, 0]*qe_red)
    ur = np.sum(img[:, 0, 1]*qe_green)
    ll = np.sum(img[:, 1, 0]*qe_green)
    lr = np.sum(img[:, 1, 1]*qe_blue)

    assert np.allclose(out, np.array([[ul, ur], [ll, lr]]))


def collect_charge_bayer_odd():
    img = np.random.uniform(low=0, high=100, size=(5, 3, 3))
    qe_red = np.random.uniform(size=5)
    qe_green = np.random.uniform(size=5)
    qe_blue = np.random.uniform(size=5)
    bayer_pattern = 'RGGB'
    out = lentil.detector.collect_charge_bayer(img, np.ones(5), qe_red,
                                               qe_green, qe_blue, bayer_pattern)

    a = np.sum(img[:, 0, 0]*qe_red)
    b = np.sum(img[:, 0, 1]*qe_green)
    c = np.sum(img[:, 0, 2]*qe_red)
    d = np.sum(img[:, 1, 0]*qe_green)
    e = np.sum(img[:, 1, 1]*qe_blue)
    f = np.sum(img[:, 1, 2]*qe_green)
    g = np.sum(img[:, 2, 0]*qe_red)
    h = np.sum(img[:, 2, 1]*qe_green)
    i = np.sum(img[:, 2, 2]*qe_red)

    assert np.allclose(out, np.array([[a, b, c], [d, e, f], [g, h, i]]))


def test_collect_charge_bayer_oversample():
    img = np.random.uniform(low=0, high=100, size=(5, 4, 4))
    qe_red = np.random.uniform(size=5)
    qe_green = np.random.uniform(size=5)
    qe_blue = np.random.uniform(size=5)
    bayer_pattern = 'RGGB'
    out = lentil.detector.collect_charge_bayer(img, np.ones(5), qe_red,
                                               qe_green, qe_blue, bayer_pattern,
                                               oversample=2)

    ul = np.einsum('ijk,i->jk', img[:, 0:2, 0:2], qe_red)
    ur = np.einsum('ijk,i->jk', img[:, 0:2, 2:4], qe_green)
    ll = np.einsum('ijk,i->jk', img[:, 2:4, 0:2], qe_green)
    lr = np.einsum('ijk,i->jk', img[:, 2:4, 2:4], qe_blue)

    assert np.allclose(out, np.block([[ul, ur], [ll, lr]]))


def test_pixel_unitary():
    a = np.random.uniform(low=0, high=1, size=(100, 100))
    b = lentil.detector.pixel(a, oversample=1)
    assert np.isclose(np.sum(a), np.sum(b))


def test_pixelate():
    # we've already tested rescale and the pixel MTF separately, so we're just
    # going to ensure the returned image has the right shape here
    img = np.ones((10, 10))
    out = lentil.detector.pixelate(img, 2)
    assert out.shape == (5, 5)


def test_adc_saturation_capacity():
    img = 10 * np.ones((10, 10))
    out = lentil.detector.adc(img, gain=1, saturation_capacity=1)
    assert np.array_equal(out, np.ones((10, 10)))


def test_adc_warn_saturate():
    with pytest.warns(UserWarning):
        img = 10 * np.ones((10, 10))
        out = lentil.detector.adc(img, gain=1, saturation_capacity=1,
                                  warn_saturate=True)


def test_adc_float_gain():
    g = np.random.uniform(size=1)
    img = 10*np.random.uniform(size=(10, 10))
    out = lentil.detector.adc(img, g, saturation_capacity=None)
    assert np.array_equal(out, np.floor(g*img))


def test_adc_array_gain():
    g = np.random.uniform(size=(10, 10))
    img = 10*np.random.uniform(size=(10, 10))
    out = lentil.detector.adc(img, g, saturation_capacity=None)
    assert np.array_equal(out, np.floor(g*img))


def test_adc_nonlinear_gain():
    g = np.random.uniform(low=0.9, high=1.5, size=4)
    img = 10*np.ones((10, 10))
    out = lentil.detector.adc(img, g, saturation_capacity=None)
    assert np.isclose(out[0, 0], np.floor(np.polyval(np.append(g, 0), 10)))


def test_adc_nonlinear_gain_cube():
    g = np.random.uniform(low=0.9, high=1.5, size=(4, 10, 10))
    img = 10*np.ones((10, 10))
    out = lentil.detector.adc(img, g, saturation_capacity=None)
    index = np.random.uniform(low=0, high=9, size=2)
    r = int(index[0])
    c = int(index[1])
    assert np.isclose(out[r, c], np.floor(np.polyval(np.append(g[:, r, c], 0), 10)))


def test_adc_invalid_gain():
    img = 10*np.ones((10, 10))
    with pytest.raises(ValueError):
        out = lentil.detector.adc(img, gain=np.ones((1, 1, 1, 1)))


def test_adc_dtype():
    g = np.random.uniform(size=1)
    img = 10*np.random.uniform(size=(10, 10))
    out = lentil.detector.adc(img, g, dtype=np.int8)
    assert out.dtype == np.int8


def test_shot_noise_poisson():
    img = np.random.uniform(low=0, high=1e5, size=(2, 2))
    shot1 = lentil.detector.shot_noise(img, method='poisson')
    shot2 = lentil.detector.shot_noise(img, method='poisson')
    assert np.all(shot1 != shot2)
    assert shot1.shape == (2, 2)


def test_shot_noise_poisson_seed():
    img = np.random.uniform(low=0, high=1e5, size=(2, 2))
    shot1 = lentil.detector.shot_noise(img, method='poisson', seed=12345)
    shot2 = lentil.detector.shot_noise(img, method='poisson', seed=12345)
    assert np.array_equal(shot1, shot2)


def test_shot_noise_poisson_negative():
    with pytest.raises(ValueError):
        lentil.detector.shot_noise(-10, method='poisson')


def test_shot_noise_poisson_big():
    with pytest.raises(ValueError):
        lentil.detector.shot_noise(10e100, method='poisson')


def test_shot_noise_gaussian():
    img = np.random.uniform(low=0, high=1e5, size=(2, 2))
    shot1 = lentil.detector.shot_noise(img, method='gaussian')
    shot2 = lentil.detector.shot_noise(img, method='gaussian')
    assert np.all(shot1 != shot2)
    assert shot1.shape == (2, 2)


def test_shot_noise_gaussian_seed():
    img = np.random.uniform(low=0, high=1e5, size=(2, 2))
    shot1 = lentil.detector.shot_noise(img, method='gaussian', seed=12345)
    shot2 = lentil.detector.shot_noise(img, method='gaussian', seed=12345)
    assert np.array_equal(shot1, shot2)


@pytest.mark.filterwarnings("ignore:invalid value encountered in sqrt")
def test_shot_noise_gaussian_negative():
    with pytest.raises(ValueError):
        lentil.detector.shot_noise(-10, method='gaussian')


def test_read_noise_gaussian():
    img = np.random.uniform(low=0, high=1e5, size=(2, 2))
    read1 = lentil.detector.read_noise(img, electrons=100)
    read2 = lentil.detector.read_noise(img, electrons=100)
    assert np.all(read1 != read2)
    assert read1.shape == (2, 2)


def test_read_noise_gaussian_seed():
    img = np.random.uniform(low=0, high=1e5, size=(2, 2))
    read1 = lentil.detector.read_noise(img, electrons=100, seed=12345)
    read2 = lentil.detector.read_noise(img, electrons=100, seed=12345)
    assert np.array_equal(read1, read2)


def test_charge_diffusion_shape():
    img = np.random.uniform(low=0, high=1e5, size=(10, 10))
    out = lentil.detector.charge_diffusion(img, sigma=0.5, oversample=1)
    assert out.shape == img.shape


def test_charge_diffusion_shape_oversample():
    img = np.random.uniform(low=0, high=1e5, size=(10, 10))
    out = lentil.detector.charge_diffusion(img, sigma=0.5, oversample=2)
    assert out.shape == img.shape


def test_dark_current_fpn():
    dark1 = lentil.detector.dark_current(rate=100, shape=(10, 10), fpn_factor=0.1, seed=12345)
    dark2 = lentil.detector.dark_current(rate=100, shape=(10, 10), fpn_factor=0.1, seed=12345)
    assert np.array_equal(dark1, dark2)


def test_rule07_dark_current_shape():
    out = lentil.detector.rule07_dark_current(temperature=25,
                                              cutoff_wavelength=6e-6,
                                              pixelscale=18e-6,
                                              shape=(2, 2))

    assert out.shape == (2, 2)


def test_rule07_dark_current_low_cutoff():
    out = lentil.detector.rule07_dark_current(temperature=110,
                                              cutoff_wavelength=4e-6,
                                              pixelscale=18e-6,
                                              shape=1)

    assert np.isclose(out, 65)


def test_rule07_dark_current_high_cutoff():
    out = lentil.detector.rule07_dark_current(temperature=110,
                                              cutoff_wavelength=5e-6,
                                              pixelscale=18e-6,
                                              shape=1)

    assert np.isclose(out, 10403)


def test_cosmic_rays():
    # we're going to punt on this one for now
    frame = lentil.detector.cosmic_rays((10, 10), 0.1*np.ones(3), 1, rate=1)
    assert frame.shape == (10, 10)
