import pytest
import numpy as np
import lentil


def test_zernike_rho_theta():
    with pytest.raises(ValueError):
        lentil.zernike(mask=1, index=1, normalize=True, rho=1, theta=None)


def test_zernike_basis():
    basis = lentil.zernike_basis(mask=np.ones((3, 3)), modes=1, vectorize=False)
    assert np.array_equal(basis, np.ones((1, 3, 3)))


def test_zernike_basis_vectorize():
    basis = lentil.zernike_basis(mask=np.ones((3, 3)), modes=1, vectorize=True)
    assert np.array_equal(basis, np.ones((1, 9)))


def test_zernike_fit():
    mask = lentil.circlemask((256, 256), 128)
    coeffs = np.random.rand(4)*100e-9
    opd = lentil.zernike_compose(mask, coeffs)
    fit_coeffs = lentil.zernike_fit(opd, mask, np.arange(1, 5))
    assert np.all(np.isclose(coeffs, fit_coeffs))


def test_zernike_remove():
    mask = lentil.circlemask((256, 256), 128)
    coeffs = np.random.rand(4)*100e-9
    opd = lentil.zernike_compose(mask, coeffs)
    residual = lentil.zernike_remove(opd, mask, np.arange(1, 5))
    assert np.all(np.isclose(residual, np.zeros_like(residual)))

def test_zernike_random_mask():
    mask = lentil.circlemask((256,256), 128)
    random_mask = np.random.uniform(low=-1, high=1, size=mask.shape) * mask

    assert np.array_equal(lentil.zernike(mask, 1), lentil.zernike(random_mask, 1))
