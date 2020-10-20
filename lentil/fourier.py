import functools
import numpy as np

__all__ = ['dft2', 'idft2']


def dft2(f, alpha, npix=None, shift=(0, 0), unitary=True, out=None):
    """Compute the 2-dimensional discrete Fourier Transform.

    Parameters
    ----------
    f : array_like
        2D array to Fourier Transform

    alpha : float or array_like
        Output plane sampling interval (frequency). If :attr:`alpha` is an
        array, ``alpha[1]`` represents row-wise sampling and ``alpha[2]``
        represents column-wise sampling. If :attr:`alpha` is a scalar,
        ``alpha[1] = alpha[2] = alpha`` gives uniform sampling across the rows
        and columns of the output plane.

    npix : int or array_like, optional
        Size of the output array :attr:`F`. If :attr:`npix` is an array,
        ``F.shape = (npix[1], npix[2])``. If :attr:`npi`` is a scalar,
        ``F.shape = (npix, npix)``. Default is ``f.shape``.

    shift : array_like, optional
        Number of pixels in (x,y) to shift the DC pixel in the output plane
        with the origin centrally located in the plane. Default is ``[0,0]``.

    unitary : bool, optional
        Normalization flag. If ``True``, a normalization is performed on the
        output such that the DFT operation is unitary and energy is conserved
        through the Fourier transform operation (Parseval's theorem). In this
        way, the energy in in a limited-area DFT is a fraction of the total
        energy corresponding to the limited area. Default is ``True``.

    out : ndarray or None
        A location into which the result is stored. If provided, it must have
        shape = npix and dtype = np.complex. If not provided or None, a 
        freshly-allocated array is returned.

    Returns
    -------
    F : complex ndarray

    Notes
    -----
    * Setting ``alpha = 1/f.shape`` and ``npix = f.shape`` is equivalent to

      ::

          F = np.fft.ifftshift(np.fft.fft2(np.fft.ifftshift(f)))

    * ``dft2()`` is designed to place the DC pixel in the same location as a
      well formed call to any standard FFT for both even and odd sized input
      arrays. The DC pixel is located at ``np.floor(npix/2) + 1``, which is
      consistent with calls to Numpy's FFT method where the input and output
      are correctly shifted:
      ``np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(f)))``.

    * If the y-axis shift behavior is not what you are expecting, you most
      likely have your plotting axes flipped (matplotlib's default behavior is
      to place [0,0] in the upper left corner of the axes). This may be resolved
      by either flipping the sign of the y component of ``shift`` or by passing
      ``origin = 'lower'`` to ``imshow()``.

    """

    alpha = np.asarray(alpha)
    if alpha.size == 1:
        ay, ax = float(alpha), float(alpha)
    else:
        ay, ax = float(alpha[0]), float(alpha[1])

    f = np.asarray(f)
    m, n = f.shape

    if npix is None:
        npix = [m, n]

    npix = np.asarray(npix)
    if npix.size == 1:
        M, N = int(npix), int(npix)
    else:
        M, N = int(npix[0]), int(npix[1])

    if isinstance(shift, np.ndarray):
        sx, sy = float(shift[0]), float(shift[1])
    else:
        sx, sy = shift

    if out is not None:
        if not np.can_cast(complex, out.dtype):
            raise TypeError(f"Cannot cast complex output to dtype('{out.dtype}')")

    E1, E2 = _dft2_matrices(m, n, M, N, ax, ay, sx, sy)
    F = np.dot(E1.dot(f), E2, out=out)

    # now calculate the answer, without reallocating memory
    if unitary:
        np.multiply(F, np.sqrt(np.abs(ax * ay)), out=F)
    
    return F


@functools.lru_cache(maxsize=32)
def _dft2_matrices(m, n, M, N, ax, ay, sx, sy):
    X, Y, U, V = _dft2_coords(m, n, M, N)
    E1 = np.exp(-2.0 * np.pi * 1j * ay * np.outer(Y-sy, V-sy)).T
    E2 = np.exp(-2.0 * np.pi * 1j * ax * np.outer(X-sx, U-sx))
    return E1, E2


@functools.lru_cache(maxsize=32)
def _dft2_coords(m, n, M, N):
    # Y and X are (r,c) coordinates in the (m x n) input plane f
    # V and U are (r,c) coordinates in the (M x N) ourput plane F

    X = np.arange(n) - np.floor(n/2.0)
    Y = np.arange(m) - np.floor(m/2.0)
    U = np.arange(N) - np.floor(N/2.0)
    V = np.arange(M) - np.floor(M/2.0)
    
    return X, Y, U, V


def idft2(F, alpha, npix=None, shift=(0, 0), unitary=True, out=None):
    """Compute the 2-dimensional inverse discrete Fourier Transform.

    Parameters
    ----------
    F : array_like
        2D array to Fourier Transform

    alpha : float or array_like
        Input plane sampling interval (frequency). If :attr:`alpha` is an array,
        ``alpha[1]`` represents row-wise sampling and ``alpha[2]`` represents
        column-wise sampling. If :attr:`alpha` is a scalar,
        ``alpha[1] = alpha[2] = alpha`` represents uniform sampling across the
        rows and columns of the input plane.

    npix : int or array_like, optional
        Size of the output array :attr:`F`. If :attr:`npix` is an array,
        ``F.shape = (npix[1], npix[2])``. If :attr:`npix` is a scalar,
        ``F.shape = (npix, npix)``. Default is ``F.shape``

    shift : array_like, optional
        Number of pixels in (x,y) to shift the DC pixel in the output plane with
        the origin centrally located in the plane. Default is `[0,0]`.

    unitary : bool, optional
        Normalization flag. If ``True``, a normalization is performed on the
        output such that the DFT operation is unitary and energy is conserved
        through the Fourier transform operation (Parseval's theorem). In this
        way, the energy in in a limited-area DFT is a fraction of the total
        energy corresponding to the limited area. Default is ``True``.

    Returns
    -------
    f : complex ndarray

    Notes
    -----
    * Setting ``alpha = 1/F.shape`` and ``npix = F.shape`` is equivalent to

      ::

          F = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(F)))

    * ``idft2()`` is designed to place the DC pixel in the same location as a
      well formed call to any standard FFT for both even and odd sized input
      arrays. The DC pixel is located at ``np.floor(npix/2) + 1``, which is
      consistent with calls to Numpy's IFFT method where the input and output
      are correctly shifted:
      ``np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(f)))``.

    * If the y-axis shift behavior is not what you are expecting, you most
      likely have your plotting axes flipped (matplotlib's default behavior is
      to place [0,0] in the upper left corner of the axes). This may be resolved
      by either flipping the sign of the y component of ``shift`` or by passing
      ``origin = 'lower'`` to ``imshow()``.

    References
    ----------
    * `Expressing the inverse DFT in terms of the DFT <https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Expressing_the_inverse_DFT_in_terms_of_the_DFT>`_.  

    """
    F = np.asarray(F)
    N = F.size
    # will allocate memory for F if out == None
    F = dft2(np.conj(F), alpha, npix, shift, unitary=unitary, out=out)
    np.conj(F, out = F)
    return np.divide(F, F.size,out=F)
    # return np.conj(dft2(np.conj(F), alpha, npix, shift, unitary))/N
