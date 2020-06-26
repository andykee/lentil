import pytest
import numpy as np
import lentil


def test_zernike_rho_theta():
    with pytest.raises(ValueError):
        lentil.zernike.zernike(mask=1, index=1, normalize=True, rho=1, theta=None)


def test_zernike_index():
    with pytest.raises(ValueError):
        lentil.zernike.zernike_index(-1)


def test_zernike_basis():
    basis = lentil.zernike.zernike_basis(mask=np.ones((3, 3)), modes=1, vectorize=False)
    assert np.array_equal(basis, np.ones((1, 3, 3)))


def test_zernike_basis_vectorize():
    basis = lentil.zernike.zernike_basis(mask=np.ones((3, 3)), modes=1, vectorize=True)
    assert np.array_equal(basis, np.ones((1, 9)))


def test_zernike_fit():
    mask = lentil.util.circlemask((256, 256), 128)
    coeffs = np.random.rand(4)*100e-9
    phase = lentil.zernike.zernike_compose(mask, coeffs)
    fit_coeffs = lentil.zernike.zernike_fit(phase, mask, np.arange(1, 5))
    assert np.all(np.isclose(coeffs, fit_coeffs))


def test_zernike_remove():
    mask = lentil.util.circlemask((256, 256), 128)
    coeffs = np.random.rand(4)*100e-9
    phase = lentil.zernike.zernike_compose(mask, coeffs)
    residual = lentil.zernike.zernike_remove(phase, mask, np.arange(1, 5))
    assert np.all(np.isclose(residual, np.zeros_like(residual)))
