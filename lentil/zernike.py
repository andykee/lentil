from math import factorial

import numpy as np

import lentil
import lentil._mesh

def zernike(mask, index, normalize=True, order='noll', rho=None, theta=None):
    """Compute the circular Zernike polynomial for a given mask.

    Parameters
    ----------
    mask : array_like
        Binary mask defining the extent to compute the Zernike polynomial
        over.
    index : int
        Zernike mode index
    normalize : bool, optional
        If True (default), the output is normalized according to [1]. If False,
        the output value ranges [-1, 1] over the mask.
    order : {'noll', 'fringe', 'ansi'}, optional
        Zernike ordering scheme. Default is 'noll'.
    rho : array_like, optional
        Radial coordinates of the mask array. :attr:`rho` should be 0 at the
        origin and 1 at the edge of the circle.
    theta : array_like, optional
        Angular coordinates of the mask array in radians.

    Returns
    -------
    ndarray

    Notes
    -----
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

    if order == 'noll':
        n, m = noll_to_nm(index)
    elif order == 'fringe':
        n, m = fringe_to_nm(index)
    elif order =='ansi':
        n, m = ansi_to_nm(index)
    else:
        raise ValueError('Unknown order', order)
    
    if rho is None:
        rho, theta = zernike_coordinates(mask)
    else:
        if theta is None:
            raise ValueError("Both rho and theta must be specified")

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
            Z = np.sqrt(2) * np.sqrt(n+1) * R(m, n, rho) * np.sin(np.abs(m)*theta) * mask
        else:
            Z = R(m, n, rho) * np.sin(np.abs(m)*theta) * mask

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
    

def zernike_compose(mask, coeffs, normalize=True, order='noll', rho=None, theta=None):
    """Create an OPD based on the supplied Zernike coefficients.

    Parameters
    ----------
    mask : array_like
        Binary mask defining the extent to compute the Zernike polynomial over.
    coeffs : array_like
        List of coefficients corresponding to Zernike indices used to create
        the OPD.
    normalize : bool, optional
        If True (default), the output is normalized according to [1]. If False,
        the output value ranges [-1, 1] over the mask.
    order : {'noll', 'fringe', 'ansi'}, optional
        Zernike ordering scheme. Default is 'noll'.
    rho : array_like, optional
        Radial coordinates of the mask array. :attr:`rho` should be 0 at the
        origin and 1 at the edge of the circle.
    theta : array_like, optional
        Angular coordinates of the mask array in radians.

    Returns
    -------
    ndarray

    Examples
    --------
    Compute a random OPD using the first ten Zernikes:

    .. plot::
        :scale: 50
        :include-source:
        :context: reset

        >>> mask = prtools.circle((256,256), 120, antialias=False)
        >>> coeffs = np.random.rand(10)*1e-8
        >>> opd = prtools.zernike_compose(mask, coeffs)
        >>> plt.imshow(opd)

    Using the same mask, compute an OPD representing 200 nm focus error (Z4) and -100 nm
    astigmatism error (Z6):

    .. plot::
        :context: close-figs
        :include-source:
        :scale: 50

        >>> opd = prtools.zernike_compose(mask, [0, 0, 0, 200e-9, 0, -100e-9])
        >>> plt.imshow(opd)

    References
    ----------
    [1] Noll, RJ. Zernike polynomials and atmospheric turbulence. J Opt Soc Am 66, 207-211  (1976).

    """
    mask = np.asarray(mask)
    coeffs = np.asarray(coeffs)
    opd = np.zeros(mask.shape)

    for index, coeff in enumerate(coeffs):
        opd += coeff * zernike(mask, index+1, normalize, order, rho, theta)

    return opd

def zernike_basis(mask, modes, vectorize=False, normalize=True, order='noll',
                  rho=None, theta=None):
    """Compute a Zernike basis set for a given mask.

    Parameters
    ----------
    mask : array_like
        Binary mask defining the extent to compute the Zernike polynomial over.
    modes : array_like
        List of modes to return.
    vectorize : bool, optional
        If True, the output is returned as a
        ``(length(modes), modes.shape[0]*modes.shape[1])`` array If False
        (default), the output is returned as a
        ``(length(terms), mask.shape[0], mask.shape[1])`` cube.
    normalize : bool, optional
        If True (default), the output is normalized according to [1]. If False,
        the output value ranges [-1, 1] over the mask.
    order : {'noll', 'fringe', 'ansi'}, optional
        Zernike ordering scheme. Default is 'noll'.
    rho : array_like, optional
        Radial coordinates of the mask array. :attr:`rho` should be 0 at the
        origin and 1 at the edge of the circle.
    theta : array_like, optional
        Angular coordinates of the mask array in radians.

    Returns
    -------
    ndarray

    References
    ----------
    [1] Noll, RJ. Zernike polynomials and atmospheric turbulence. J Opt Soc Am 66, 207-211  (1976).

    """
    modes = np.asarray(modes)
    if modes.shape == ():
        modes = modes[..., np.newaxis]

    mask = np.asarray(mask)

    basis = []
    for index, mode in enumerate(modes):
        basis.append(zernike(mask, mode, normalize, order, rho, theta))
    basis = np.asarray(basis)

    if vectorize:
        # reshape basis from cube to matrix
        return basis.reshape(basis.shape[0], -1)
    else:
        return basis


def zernike_fit(opd, mask, modes, normalize=True, order='noll', rho=None,
                theta=None):
    """Fit a Zernike basis set to an OPD.

    Parameters
    ----------
    opd : array_like
        OPD to fit. NaN values in ``opd`` are excluded from the fit.
    mask : array_like
        Binary mask defining the extent to compute the Zernike basis over.
    modes : array_like
        List of modes to fit.
    normalize : bool, optional
        If True (default), the output is normalized according to [1]. If False,
        the output value ranges [-1, 1] over the mask.
    order : {'noll', 'fringe', 'ansi'}, optional
        Zernike ordering scheme. Default is 'noll'.
    rho : array_like, optional
        Radial coordinates of the mask array. :attr:`rho` should be 0 at the
        origin and 1 at the edge of the circle.
    theta : array_like, optional
        Angular coordinates of the mask array in radians.

    Returns
    -------
    ndarray

    Example
    -------
    .. code:: pycon

        >>> mask = prtools.circle((256,256),128)
        >>> coeffs = np.random.rand(4)*100e-9
        >>> opd = prtools.zernike_compose(mask, coeffs)
        >>> fit_coeffs = prtools.zernike_fit(opd, mask, np.arange(2,4))
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
    opd = np.asarray(opd).ravel()
    mask = np.asarray(mask)
    basis = zernike_basis(mask, modes, True, normalize, order, rho, theta)
    
    # valid data includes all non-NaN values inside the mask
    valid = np.isfinite(opd) & mask.ravel().astype(bool)
    
    coeff, *_ = np.linalg.lstsq(basis[:, valid].T, opd[valid], rcond=None)

    return coeff


def zernike_remove(opd, mask, modes, order='noll', rho=None, theta=None):
    """Fit and remove a Zernike basis set from an OPD.

    Parameters
    ----------
    opd : array_like
        OPD to fit.
    mask : array_like
        Binary mask defining the extent to compute the Zernike basis over.
    modes : array_like
        List of modes to remove.
    order : {'noll', 'fringe', 'ansi'}, optional
        Zernike ordering scheme. Default is 'noll'.
    rho : array_like, optional
        Radial coordinates of the mask array. :attr:`rho` should be 0 at the
        origin and 1 at the edge of the circle.
    theta : array_like, optional
        Angular coordinates of the mask array in radians.

    Returns
    -------
    ndarray

    See Also
    --------
    :func:`zernike_fit`

    References
    ----------
    [1] Noll, RJ. Zernike polynomials and atmospheric turbulence. J Opt Soc Am 66, 207-211  (1976).

    """
    opd = np.asarray(opd)
    mask = np.asarray(mask)

    coeffs = zernike_fit(opd, mask, modes, True, order, rho, theta)
    fit_opd = zernike_compose(mask, coeffs, True, order, rho, theta)

    residual = opd - fit_opd

    return residual


def noll_to_nm(j):
    """Convert Noll index to (n, m) index.

    Parameters
    ----------
    j : int
        Noll Zernike index

    Returns
    -------
    n : int
        Radial degree.
    m : int
        Azimuthal degree.
    """
    if j < 1:
        raise ValueError('Noll index must be greater than zero')

    # find the row j is in
    n = int(np.ceil((-1 + np.sqrt(1 + 8 * j)) / 2) - 1)
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

        for _ in range(n // 2):
            row_m.append(row_m[-1] + 2)
            row_m.append(row_m[-1])

        m = row_m[r] * sign

    return int(n), int(m)


def ansi_to_nm(j):
    """Convert ANSI index to (n, m) index.

    Parameters
    ----------
    j : int
        ANSI index

    Returns
    -------
    n : int
        Radial degree
    m : int
        Azimuthal degree
    """
    if j < 0:
        raise ValueError('ANSI index must be zero or larger')
    
    # https://en.wikipedia.org/wiki/Zernike_polynomials
    n = np.floor((np.sqrt(8 * j + 1) - 1) / 2)
    m = 2*j - n*(n + 2)
    return int(n), int(m)


def fringe_to_nm(j):
    """Convert Fringe index to (n, m) index.

    Parameters
    ----------
    j : int
        Fringe index

    Returns
    -------
    n : int
        Radial degree
    m : int
        Azimuthal degree
    """
    if j < 1:
        raise ValueError('Fringe index must be greater than zero')
    
    # https://github.com/arielmission-space/PAOS/blob/main/paos/classes/zernike.py
    m_n = 2 * (np.ceil(np.sqrt(j)) - 1)
    g_s = (m_n / 2) ** 2 + 1
    n = m_n / 2 + np.floor((j - g_s) / 2)
    m = (m_n - n) * (1 - np.mod(j - g_s, 2) * 2)
    return int(n), int(m)


def zernike_coordinates(mask, shift=None, angle=0, indexing='xy'):
    """Compute the Zernike coordinate system for a given mask.
    
    Parameters
    ----------
    mask : array_like
        Binary mask defining the extent to compute the Zernike polynomial 
        over.
    shift : (2,) array_like, optional
        How far to shift center in float (rows, cols). If None (default),
        shift is computed automatically to locate the origin at the mask
        centroid.
    angle: float, optional
        Angle to rotate coordinate system by in degrees. rotate is specified 
        relative to the x-axis. Default is 0.
    indexing : {'ij', 'xy'}, optional
        Cartesian ('xy', default) or matrix ('ij') or indexing of shift.
    
    Returns
    -------
    rho, theta : ndarrays
    
    """
    mask = np.asarray(mask, dtype=bool)

    if shift is None:
        shift = lentil.centroid(mask, mode='center', indexing=indexing)

    # since the Zernike coordinates are always computed in xy indexing, we
    # need to transform shift into xy indexing if it was specified in ij
    if indexing == 'ij':
        shift = (shift[1], -shift[0])

    x1,x2 = lentil._mesh.mesh(shape=mask.shape, shift=shift, angle=angle, indexing='xy')

    r = np.abs(x1+1j*x2)
    rho = r/np.max(r*mask)  # rho is defined to be 1 on the edge of the aperture
    theta = np.angle(x1 + 1j*x2)
    
    rho[np.where(rho==0)] = 1e-99
    theta[np.where(theta==0)] = 1e-99

    return rho, theta