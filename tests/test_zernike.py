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
    mask = lentil.circle((256, 256), 128, antialias=False)
    coeffs = np.random.rand(4)*100e-9
    opd = lentil.zernike_compose(mask, coeffs)
    fit_coeffs = lentil.zernike_fit(opd, mask, np.arange(1, 5))
    assert np.all(np.isclose(coeffs, fit_coeffs))


def test_zernike_remove():
    mask = lentil.circle((256, 256), 128, antialias=False)
    coeffs = np.random.rand(4)*100e-9
    opd = lentil.zernike_compose(mask, coeffs)
    residual = lentil.zernike_remove(opd, mask, np.arange(1, 5))
    assert np.all(np.isclose(residual, np.zeros_like(residual)))

def test_zernike_random_mask():
    mask = lentil.circle((256, 256), 128, antialias=False)
    random_mask = np.random.uniform(low=-1, high=1, size=mask.shape) * mask

    assert np.array_equal(lentil.zernike(mask, 1), lentil.zernike(random_mask, 1))


# Sample points on a 65x65 circular mask centered at (32, 32):
# right edge, top edge, and upper-right diagonal (45 deg)
RIGHT, TOP, UPPER_RIGHT = (32, 60), (4, 32), (12, 52)


@pytest.mark.parametrize('index, positive_at', [
    (2, RIGHT),         # x tilt
    (3, TOP),           # y tilt
    (5, UPPER_RIGHT),   # 45 deg astigmatism
    (7, TOP),           # y coma
    (8, RIGHT),         # x coma
])
def test_zernike_orientation(index, positive_at):
    # theta is measured CCW from +x (increasing column) toward +y
    # (decreasing row)
    mask = lentil.circle((65, 65), 32, antialias=False)
    z = lentil.zernike(mask, index)
    assert z[positive_at] > 0.5


def test_zernike_coordinates_origin_at_mask_centroid():
    mask = lentil.circle((128, 128), 30, shift=(-20, 30), antialias=False)

    rho_xy, theta_xy = lentil.zernike_coordinates(mask, indexing='xy')
    rho_ij, theta_ij = lentil.zernike_coordinates(mask, indexing='ij')

    assert np.unravel_index(np.argmin(rho_xy), rho_xy.shape) == (44, 94)
    assert np.allclose(rho_xy, rho_ij)
    assert np.allclose(theta_xy, theta_ij)


def test_zernike_coordinates_shift_indexing():
    mask = lentil.circle((128, 128), 30, shift=(-20, 30), antialias=False)

    rho_auto, theta_auto = lentil.zernike_coordinates(mask)
    rho_xy, theta_xy = lentil.zernike_coordinates(mask, shift=(30, 20), indexing='xy')
    rho_ij, theta_ij = lentil.zernike_coordinates(mask, shift=(-20, 30), indexing='ij')

    assert np.allclose(rho_xy, rho_auto)
    assert np.allclose(theta_xy, theta_auto)
    assert np.allclose(rho_ij, rho_auto)
    assert np.allclose(theta_ij, theta_auto)
