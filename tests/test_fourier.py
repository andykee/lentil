import numpy as np

import lentil.fourier


def test_dft2_even():
    n = 10
    f = np.random.rand(n, n) + 1j * np.random.rand(n, n)

    F_dft = lentil.fourier.dft2(f, 1/n, unitary=False)
    F_fft = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(f)))

    assert np.allclose(F_dft, F_fft)


def test_idft2_even():
    n = 10
    F = np.random.rand(n, n) + 1j * np.random.rand(n, n)

    f_dft = lentil.fourier.idft2(F, 1/n, unitary=False)
    f_fft = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(F)))

    assert np.allclose(f_dft, f_fft)


def test_dft2_odd():
    n = 11
    f = np.random.rand(n, n) + 1j * np.random.rand(n, n)

    F_dft = lentil.fourier.dft2(f, 1/n, unitary=False)
    F_fft = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(f)))

    assert np.allclose(F_dft, F_fft)


def test_idft2_odd():
    n = 11
    F = np.random.rand(n, n) + 1j * np.random.rand(n, n)

    f_dft = lentil.fourier.idft2(F, 1/n, unitary=False)
    f_fft = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(F)))

    assert np.allclose(f_dft, f_fft)


def test_dft2_inverse():
    n = 10
    f = np.random.rand(n, n) + 1j * np.random.rand(n, n)

    F_dft = lentil.fourier.dft2(f, 1/n, unitary=False)
    f_dft = lentil.fourier.idft2(F_dft, 1/n, unitary=False)

    assert np.allclose(f, f_dft)


def test_dft2_unitary():
    n = 10
    f = np.random.rand(n, n) + 1j * np.random.rand(n, n)
    f_power = np.sum(np.abs(f)**2)

    F = lentil.fourier.dft2(f, 1/n, unitary=True)
    F_power = np.sum(np.abs(F)**2)

    assert np.allclose(f_power, F_power)


def test_dft2_shift():
    n = 100
    shift = np.round(np.random.uniform(low=-25, high=25, size=2))
    f = np.ones((n, n))

    F = lentil.fourier.dft2(f, 1/n, shift=shift)

    (xc, yc) = (np.floor(n/2), np.floor(n/2))
    (x, y) = lentil.util.centroid(np.abs(F)**2)
    observed_shift = (x-xc, y-yc)
    assert np.array_equal(shift, observed_shift)


def test_dft2_out():
    n = 10
    f = np.random.normal(size=(n,n))
    F = lentil.fourier.dft2(f, 1/n, out=f)
    assert f is F
