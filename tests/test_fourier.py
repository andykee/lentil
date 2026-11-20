import numpy as np
import pytest

import lentil.fourier
import lentil.helper


# every test runs against both the matrix triple product DFT and the chirp
# Z-transform, which share the same interface and must produce the same result
TRANSFORMS = [
    pytest.param((lentil.fourier.dft2, lentil.fourier.idft2), id='dft2'),
    pytest.param((lentil.fourier.czt2, lentil.fourier.iczt2), id='czt2'),
]


@pytest.fixture(params=TRANSFORMS)
def transform(request):
    return request.param


@pytest.fixture
def fwd(transform):
    return transform[0]


@pytest.fixture
def inv(transform):
    return transform[1]


def test_fwd_even(fwd):
    n = 10
    f = np.random.rand(n, n) + 1j * np.random.rand(n, n)

    F_dft = fwd(f, 1/n, unitary=False)
    F_fft = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(f)))

    assert np.allclose(F_dft, F_fft)


def test_inv_even(inv):
    n = 10
    F = np.random.rand(n, n) + 1j * np.random.rand(n, n)

    f_dft = inv(F, 1/n, unitary=False)
    f_fft = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(F)))

    assert np.allclose(f_dft, f_fft)


def test_fwd_odd(fwd):
    n = 11
    f = np.random.rand(n, n) + 1j * np.random.rand(n, n)

    F_dft = fwd(f, 1/n, unitary=False)
    F_fft = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(f)))

    assert np.allclose(F_dft, F_fft)


def test_inv_odd(inv):
    n = 11
    F = np.random.rand(n, n) + 1j * np.random.rand(n, n)

    f_dft = inv(F, 1/n, unitary=False)
    f_fft = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(F)))

    assert np.allclose(f_dft, f_fft)


def test_inverse(fwd, inv):
    n = 10
    f = np.random.rand(n, n) + 1j * np.random.rand(n, n)

    F_dft = fwd(f, 1/n, unitary=False)
    f_dft = inv(F_dft, 1/n, unitary=False)

    assert np.allclose(f, f_dft)


@pytest.mark.parametrize('shape', [(10, 10), (11, 11), (10, 13)])
def test_inverse_unitary(fwd, inv, shape):
    m, n = shape
    f = np.random.rand(m, n) + 1j * np.random.rand(m, n)

    F_dft = fwd(f, [1/m, 1/n], unitary=True)
    f_dft = inv(F_dft, [1/m, 1/n], unitary=True)

    assert np.allclose(f, f_dft)


def test_inverse_zoom(fwd, inv):
    # m x n input zero padded to an M x N period at alpha = 1/(M, N)
    m, n = 8, 6
    M, N = 20, 15
    f = np.random.rand(m, n) + 1j * np.random.rand(m, n)

    F = fwd(f, [1/M, 1/N], shape=(M, N), unitary=False)
    f_inv = inv(F, [1/M, 1/N], shape=(m, n), unitary=False)

    assert np.allclose(f, f_inv)


def test_inv_cropped_spectrum(inv):
    # inverting a cropped spectrum must give the same result as inverting the
    # full spectrum with the discarded frequencies set to zero
    n, k = 16, 8
    F = np.random.rand(n, n) + 1j * np.random.rand(n, n)
    crop = slice((n-k)//2, (n+k)//2)

    F_masked = np.zeros_like(F)
    F_masked[crop, crop] = F[crop, crop]

    f_cropped = inv(F[crop, crop], 1/n, shape=n, unitary=False)
    f_masked = inv(F_masked, 1/n, unitary=False)
    f_fft = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(F_masked)))

    assert np.allclose(f_cropped, f_masked)
    assert np.allclose(f_cropped, f_fft)


def test_fwd_unitary(fwd):
    n = 10
    f = np.random.rand(n, n) + 1j * np.random.rand(n, n)
    f_power = np.sum(np.abs(f)**2)

    F = fwd(f, 1/n, unitary=True)
    F_power = np.sum(np.abs(F)**2)

    assert np.allclose(f_power, F_power)


def test_inv_unitary(inv):
    n = 10
    F = np.random.rand(n, n) + 1j * np.random.rand(n, n)
    F_power = np.sum(np.abs(F)**2)

    f = inv(F, 1/n, unitary=True)
    f_power = np.sum(np.abs(f)**2)

    assert np.allclose(f_power, F_power)


def test_fwd_shift(fwd):
    n = 100
    shift = np.round(np.random.uniform(low=-25, high=25, size=2))
    f = np.ones((n, n))

    F = fwd(f, 1/n, shift=shift)

    (rc, cc) = (np.floor(n/2), np.floor(n/2))
    (r, c) = lentil.util.centroid(np.abs(F)**2, indexing='ij')
    observed_shift = (r-rc, c-cc)
    assert np.array_equal(shift, observed_shift)


def test_fwd_offset(fwd):
    n = 10
    m = 3

    f = np.zeros((n,n), dtype=complex)
    r,c = np.random.randint(0, n-m, size=2)
    f[r:r+m, c:c+m] = np.random.rand(m, m) + 1j * np.random.rand(m,m)
    slc = lentil.helper.boundary_slice(f)
    offset = lentil.helper.slice_offset(slc, f.shape)

    F = fwd(f, alpha=1/m, shape=10)
    FF = fwd(f[slc], alpha=1/m, shape=10, offset=offset)

    assert np.allclose(F, FF)


def test_fwd_out(fwd):
    n = 10
    f = np.random.normal(size=(n,n)).astype(complex)
    F_expected = fwd(f.copy(), 1/n)
    F = fwd(f, 1/n, out=f)
    assert f is F
    assert np.allclose(F, F_expected)


def test_fwd_out_bad_dtype(fwd):
    n = 10
    f = np.random.normal(size=(n,n))
    with pytest.raises(TypeError):
        fwd(f, 1/n, out=np.empty((n, n)))


def test_fwd_rect(fwd):
    m, n = 10, 11
    f = np.random.rand(m, n) + 1j * np.random.rand(m, n)

    F_dft = fwd(f, [1/m, 1/n], unitary=False)
    F_fft = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(f)))

    assert np.allclose(F_dft, F_fft)


def _random_case(rng):
    m, n = rng.integers(1, 64, size=2)
    M, N = rng.integers(1, 96, size=2)
    alpha = rng.uniform(0.001, 1.0, size=2) * rng.choice([-1, 1], size=2)
    shift = rng.uniform(-20, 20, size=2)
    offset = rng.uniform(-20, 20, size=2)
    f = rng.standard_normal((m, n)) + 1j * rng.standard_normal((m, n))
    return f, dict(alpha=alpha, shape=(M, N), shift=shift, offset=offset)


@pytest.mark.parametrize('seed', range(50))
def test_czt2_matches_dft2(seed):
    rng = np.random.default_rng(seed)
    f, kwargs = _random_case(rng)

    F_dft = lentil.fourier.dft2(f, **kwargs)
    F_czt = lentil.fourier.czt2(f, **kwargs)

    assert F_czt.shape == F_dft.shape
    assert np.max(np.abs(F_czt - F_dft)) <= 1e-10 * np.max(np.abs(F_dft))


@pytest.mark.parametrize('seed', range(50))
def test_iczt2_matches_idft2(seed):
    rng = np.random.default_rng(seed)
    F, kwargs = _random_case(rng)
    kwargs.pop('offset')  # idft2 has no offset argument

    f_dft = lentil.fourier.idft2(F, **kwargs)
    f_czt = lentil.fourier.iczt2(F, **kwargs)

    assert f_czt.shape == f_dft.shape
    assert np.max(np.abs(f_czt - f_dft)) <= 1e-10 * np.max(np.abs(f_dft))


@pytest.mark.parametrize('shape_in, shape_out', [
    ((1, 1), (1, 1)),
    ((1, 7), (5, 1)),
    ((2, 2), (9, 4)),
    ((64, 8), (4, 128)),   # exercises both pass orderings
    ((8, 64), (128, 4)),
])
def test_czt2_matches_dft2_shapes(shape_in, shape_out):
    rng = np.random.default_rng(0)
    f = rng.standard_normal(shape_in) + 1j * rng.standard_normal(shape_in)
    kwargs = dict(alpha=[0.07, 0.03], shape=shape_out, shift=(1.5, -0.25),
                  offset=(0.5, 2))

    F_dft = lentil.fourier.dft2(f, **kwargs)
    F_czt = lentil.fourier.czt2(f, **kwargs)

    assert np.max(np.abs(F_czt - F_dft)) <= 1e-10 * np.max(np.abs(F_dft))


def test_czt2_workers():
    n = 32
    f = np.random.rand(n, n) + 1j * np.random.rand(n, n)
    assert np.allclose(lentil.fourier.czt2(f, 1/n, workers=2),
                       lentil.fourier.czt2(f, 1/n))
