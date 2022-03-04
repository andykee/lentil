from math import factorial

import numpy as np

import lentil

__all__ = ['zernike', 'zernike_compose', 'zernike_fit', 'zernike_remove',
           'zernike_basis', 'zernike_coordinates']


def zernike(mask, index, normalize=True, rho=None, theta=None):
    """Compute the circular Zernike polynomial for a given mask.

    Parameters
    ----------
    mask : array_like
        Mask defining the extent to compute the Zernike polynomial over. All
        nonzero entries are included in the result.

    index : int
        Noll Zernike index as defined in [1]

    normalize : bool, optional
        If True (default), the output is normalized according to [1]. If False,
        the output value ranges [-1, 1] over the mask.

    rho : array_like, optional
        Radial coordinates of the mask array. :attr:`rho` should be 0 at the
        origin and 1 at the edge of the circle.

    theta : array_like, optional
        Angular coordinates of the mask array in radians.

    Returns
    -------
    out
        Circular Zernike polynomial computed over the given mask.

    Warning
    -------
    Zernike polynomials are defined to be orthogonal on the unit circle. If
    the supplied mask is non-circular, the Zernike polynomial is computed on an
    outscribing circle and then cropped by the mask. Note that this operation
    breaks the orthogonality of the Zernike polynomial. When working with a
    non-circular mask, care must be taken to understand any side-effects of
    using Zernike polynomials constructed in this manner. Any undesirable
    side-effects can be mitigated by using Zernike polynomials that have
    undergone Gram-Schmidt orthogonalization over the supplied mask.

    References
    ----------
    [1] Noll, RJ. Zernike polynomials and atmospheric turbulence. J Opt Soc Am 66, 207-211  (1976).

    """
    mask = np.asarray(mask, dtype=bool)

    #out = np.zeros_like(mask)
    #mask_slice = util.boundary_slice(mask, pad=0)
    #mask = mask[mask_slice]

    if rho is None:
        rho, theta = zernike_coordinates(mask)
    else:
        if theta is None:
            raise ValueError("Both rho and theta must be specified")

    m, n = zernike_index(index)

    if m == 0:
        if n == 0:
            Z = mask
        else:
            if normalize:
                Z = np.sqrt(n+1) * R(m, n, rho) * mask
            else:
                Z = R(m, n, rho) * mask

    elif m > 0:
        if normalize:
            Z = np.sqrt(2) * np.sqrt(n+1) * R(m, n, rho) * np.cos(m*theta) * mask
        else:
            Z = R(m, n, rho) * np.cos(m*theta) * mask

    else:
        if normalize:
            Z = np.sqrt(2) * np.sqrt(n+1) * R(m, n, rho) * np.sin(m*theta) * mask
        else:
            Z = R(m, n, rho) * np.sin(m*theta) * mask

    #out[mask_slice] = Z
    out = Z
    return out


def R(m, n, rho):

    m = int(np.abs(m))
    n = int(np.abs(n))

    if (n - m) & 1:  # odd
        return 0
    else:
        R = np.zeros(rho.shape)
        for k in range(int(n-m)//2 + 1):
            Rk = ((-1) ** k * factorial(n-k) /
                  (factorial(k) * factorial((n+m)//2-k) * factorial((n-m)//2-k)))
            R += Rk * rho ** (n-2*k)
        return R


def zernike_compose(mask, coeffs, normalize=True, rho=None, theta=None):
    """Create an OPD based on the supplied Zernike coefficients.

    Parameters
    ----------
    mask : array_like
        Mask defining the extent to compute the Zernike polynomial over. All
        nonzero entries are included in the result.

    coeffs : array_like
        List of coefficients corresponding to Zernike indices (Noll ordering)
        used to create the OPD.

    normalize : bool, optional
        If True (default), the output is normalized according to [1]. If False,
        the output value ranges [-1, 1] over the mask.

    rho : array_like, optional
        Radial coordinates of the mask array. :attr:`rho` should be 0 at the
        origin and 1 at the edge of the circle.

    theta : array_like, optional
        Angular coordinates of the mask array in radians.

    Returns
    -------
    ndarray
        OPD

    Examples
    --------
    Compute a random OPD using the first ten Zernikes:

    .. plot::
        :scale: 50
        :include-source:
        :context: reset

        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>> import lentil
        >>> mask = lentil.circlemask((256,256), 120)
        >>> coeffs = np.random.rand(10)*1e-8
        >>> opd = lentil.zernike_compose(mask, coeffs)
        >>> plt.imshow(opd, origin='lower')

    Using the same mask, compute an OPD representing 200 nm focus error (Z4) and -100 nm
    astigmatism error (Z6):

    .. plot::
        :context: close-figs
        :include-source:
        :scale: 50

        >>> opd = lentil.zernike_compose(mask, [0, 0, 0, 200e-9, 0, -100e-9])
        >>> plt.imshow(opd, origin='lower')

    References
    ----------
    [1] Noll, RJ. Zernike polynomials and atmospheric turbulence. J Opt Soc Am 66, 207-211  (1976).

    """
    mask = np.asarray(mask)
    coeffs = np.asarray(coeffs)
    opd = np.zeros(mask.shape)

    for index, coeff in np.ndenumerate(coeffs):
        opd += coeff * zernike(mask, index[0]+1, normalize, rho, theta)

    return opd


def zernike_basis(mask, modes, vectorize=False, normalize=True, rho=None, theta=None):
    """Compute a Zernike basis set for a given mask.

    Parameters
    ----------
    mask : array_like
        Mask defining the extent to compute the Zernike polynomial over. All
        nonzero entries are included in the result.

    modes : array_like
        List of modes (Noll ordering) to return.

    vectorize : bool, optional
        If True, the output is returned as a
        ``(length(modes), modes.shape[0]*modes.shape[1])`` array If False
        (default), the output is returned as a
        ``(length(terms), mask.shape[0], mask.shape[1])`` cube.

    normalize : bool, optional
        If True (default), the output is normalized according to [1]. If False,
        the output value ranges [-1, 1] over the mask.

    rho : array_like, optional
        Radial coordinates of the mask array. :attr:`rho` should be 0 at the
        origin and 1 at the edge of the circle.

    theta : array_like, optional
        Angular coordinates of the mask array in radians.

    Returns
    -------
    ndarray
        Zernike basis set

    References
    ----------
    [1] Noll, RJ. Zernike polynomials and atmospheric turbulence. J Opt Soc Am 66, 207-211  (1976).

    """
    modes = np.asarray(modes)
    if modes.shape == ():
        modes = modes[..., np.newaxis]

    mask = np.asarray(mask)
    basis = np.zeros(modes.shape + mask.shape)

    for index, mode in np.ndenumerate(modes):
        basis[index] = zernike(mask, mode, normalize, rho, theta)

    if vectorize:
        # reshape basis from cube to matrix
        return basis.reshape(basis.shape[0], -1)
    else:
        return basis


def zernike_fit(opd, mask, modes, normalize=True, rho=None, theta=None):
    """Fit a Zernike basis set to an OPD.

    Parameters
    ----------
    opd : array_like
        OPD to fit.

    mask : array_like
        Mask defining the extent to compute the Zernike polynomial over. All
        nonzero entries are included in the result.

    modes : array_like
        List of modes (Noll ordering) to fit.

    normalize : bool, optional
        If True (default), the output is normalized according to [1]. If False,
        the output value ranges [-1, 1] over the mask.

    rho : array_like, optional
        Radial coordinates of the mask array. :attr:`rho` should be 0 at the
        origin and 1 at the edge of the circle.

    theta : array_like, optional
        Angular coordinates of the mask array in radians.

    Returns
    -------
    ndarray
        List of coefficients fit to the supplied OPD over the specified number
        of Zernike modes.

    Example
    -------
    .. code:: pycon

        >>> import numpy as np
        >>> import lentil
        >>> mask = lentil.circlemask((256,256),128)
        >>> coeffs = np.random.rand(4)*100e-9
        >>> opd = lentil.zernike_compose(mask, coeffs)
        >>> fit_coeffs = lentil.zernike_fit(opd, mask, np.arange(2,4))
        >>> print('Tip/tilt coefficients:', coeffs[1:3])
        >>> print('Fit tip/tilt coefficients:', fit_coeffs)

        Tip/tilt coefficients: [9.69097470e-08 9.94332699e-08]
        Fit tip/tilt coefficients: [9.69545890e-08 9.94781119e-08]

    See Also
    --------
    :func:`zernike_remove` Fit and remove a Zernike basis set from an OPD.

    References
    ----------
    [1] Noll, RJ. Zernike polynomials and atmospheric turbulence. J Opt Soc Am 66, 207-211  (1976).

    """
    opd = np.asarray(opd)
    mask = np.asarray(mask)

    basis = zernike_basis(mask, modes, True, normalize, rho, theta)

    basis = np.linalg.pinv(basis)

    return np.einsum('ij,i->j', basis, opd.ravel())


def zernike_remove(opd, mask, modes, rho=None, theta=None):
    """Fit and remove a Zernike basis set from an OPD.

    Parameters
    ----------
    opd : array_like
        OPD to fit.

    mask : array_like
        Mask defining the extent to compute the Zernike polynomial over. All
        nonzero entries are included in the result.

    modes : array_like
        List of modes (Noll ordering) to remove.

    rho : array_like, optional
        Radial coordinates of the mask array. :attr:`rho` should be 0 at the
        origin and 1 at the edge of the circle.

    theta : array_like, optional
        Angular coordinates of the mask array in radians.

    Returns
    -------
    ndarray
        Residual OPD after the specified Zernike modes have been fit and
        removed.

    See Also
    --------
    :func:`zernike_fit`

    References
    ----------
    [1] Noll, RJ. Zernike polynomials and atmospheric turbulence. J Opt Soc Am 66, 207-211  (1976).

    """
    opd = np.asarray(opd)
    mask = np.asarray(mask)

    coeffs = zernike_fit(opd, mask, modes, rho, theta)
    fit_opd = zernike_compose(mask, coeffs, rho, theta)

    residual = opd - fit_opd

    return residual


def zernike_index(j):
    """Convert a 1-D Noll index to a 2-D radiul and azimuthal index.

    Parameters
    ----------
    j : int
        Noll Zernike index.

    Returns
    -------
    n : int
        Radial Zernike index.

    m : int
        Azimuthal Zernike index.

    """
    if j < 1:
        raise ValueError('Zernike index j must be a positive integer')

    # find the row j is in
    n = int(np.ceil((-1 + np.sqrt(1 + 8*j)) / 2) - 1)
    if n == 0:
        m = 0
    else:
        # find where in the row j is
        k = (n + 1) * (n + 2) / 2
        r = int(j - k - 1)

        if j & 1:  # odd
            sign = -1
        else:
            sign = 1

        if n & 1:  # odd
            row_m = [1, 1]
        else:
            row_m = [0]

        for i in range(int(np.floor(n / 2))):
            row_m.append(row_m[-1] + 2)
            row_m.append(row_m[-1])

        m = row_m[r] * sign

    return m, n


def zernike_coordinates(mask, shift=None, rotate=0):
    """Compute the Zernike coordinate system for a given mask.

    Parameters
    ----------
    mask : array_like
        Mask defining the extent to compute the Zernike polynomial over. All
        nonzero entries are included in the result.

    shift : (2,) array_like or None, optional
        x, y shift of the coordinate system origin in pixels. If None
        (default), shift is computed automatically to locate the origin at the
        mask centroid.

    rotate: float, optional
        Angle to rotate coordinate system by in degrees. rotate is specified
        relative to the x-axis. Default is 0.

    """

    mask = np.asarray(mask, dtype=bool)

    if shift is None:
        center = np.asarray(mask.shape)/2  # center in (r, c)
        centroid = lentil.centroid(mask)     # centroid in (r, c)
        shift = (centroid[0]-center[0], centroid[1]-center[1])

    rr, cc = lentil.helper.mesh(mask.shape, shift)

    r = np.abs(rr+1j*cc)
    rho = r/np.max(r*mask)  # rho is defined to be 1 on the edge of the aperture

    # there is a 90 degree offset since np.angle places 0 degrees up
    angle = (90-rotate) * np.pi/180

    # rr is negative here since +y is in the -row direction
    theta = np.angle(-rr*np.exp(1j*angle) + 1j*cc*np.exp(1j*angle))

    return rho, theta
