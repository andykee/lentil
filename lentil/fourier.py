import functools
import numpy as np


def dft2(f, alpha, npix=None, shift=(0, 0), offset=(0, 0), unitary=True, out=None):
    r"""Compute the 2-dimensional discrete Fourier Transform.

    The DFT is defined in one dimension as

    .. math::

        A_k = \sum_{m=0}^{n-1} a_n \exp \left(-2\pi i\frac{mk}{n}\right) \hspace{3em} k=0, \dots, n-1

    This function allows independent control over input shape, output shape,
    and output sampling by implementing the matrix triple product algorithm
    described in [1].

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
        Number of pixels in (r,c) to shift the DC pixel in the output plane
        with the origin centrally located in the plane. Default is ``(0,0)``.
    offset : array_like, optional
        Number of pixels in (r,c) that the input plane is shifted relative to
        the origin. Default is ``(0,0)``.
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

          F = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(f)))

    * ``dft2()`` is designed to place the DC pixel in the same location as a
      well formed call to any standard FFT for both even and odd sized input
      arrays. The DC pixel is located at ``np.floor(npix/2) + 1``, which is
      consistent with calls to Numpy's FFT method where the input and output
      are correctly shifted:
      ``np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(f)))``.

    * If the y-axis shift behavior is not what you are expecting, you most
      likely have your plotting axes flipped (matplotlib's default behavior is
      to place [0,0] in the upper left corner of the axes). This may be resolved
      by either flipping the sign of the y component of ``shift`` or by passing
      ``origin = 'lower'`` to ``imshow()``.

    References
    ----------
    [1] Soummer, et. al. Fast computation of Lyot-style coronagraph propagation (2007)

    """
    alpha_row, alpha_col = np.broadcast_to(alpha, (2,))

    f = np.asarray(f)
    m, n = f.shape

    if npix is None:
        npix = [m, n]
    M, N = np.broadcast_to(npix, (2,))

    shift_row, shift_col = np.broadcast_to(shift, (2,))
    offset_row, offset_col = np.broadcast_to(offset, (2,))

    if out is not None:
        if not np.can_cast(complex, out.dtype):
            raise TypeError(f"Cannot cast complex output to dtype('{out.dtype}')")

    E1, E2 = _dft2_matrices(m, n, M, N, alpha_row, alpha_col, shift_row, shift_col,
                            offset_row, offset_col)
    F = np.dot(E1.dot(f), E2, out=out)

    # now calculate the answer, without reallocating memory
    if unitary:
        np.multiply(F, np.sqrt(np.abs(alpha_row * alpha_col)), out=F)

    return F


def _dft2_matrices(m, n, M, N, alphar, alphac, shiftr, shiftc, offsetr, offsetc):
    R, S, U, V = _dft2_coords(m, n, M, N)
    E1 = np.exp(-2.0 * 1j * np.pi * alphar * np.outer(R-shiftr+offsetr, U-shiftr)).T
    E2 = np.exp(-2.0 * 1j * np.pi * alphac * np.outer(S-shiftc+offsetc, V-shiftc))
    return E1, E2


@functools.lru_cache(maxsize=32)
def _dft2_coords(m, n, M, N):
    # R and S are (r,c) coordinates in the (m x n) input plane f
    # V and U are (r,c) coordinates in the (M x N) output plane F

    R = np.arange(m) - np.floor(m/2.0)
    S = np.arange(n) - np.floor(n/2.0)
    U = np.arange(M) - np.floor(M/2.0)
    V = np.arange(N) - np.floor(N/2.0)

    return R, S, U, V


def idft2(F, alpha, npix=None, shift=(0,0), unitary=True, out=None):
    r"""Compute the 2-dimensional inverse discrete Fourier Transform.

    The IDFT is defined in one dimension as

    .. math::

        a_m = \frac{1}{n}\sum_{k=0}^{n-1} A_k \exp \left(2\pi i\frac{mk}{n}\right) \hspace{3em} m=0, \dots, n-1

    This function allows independent control over input shape, output shape,
    and output sampling by implementing the matrix triple product algorithm
    described in [1].

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

          F = np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(F)))

    * ``idft2()`` is designed to place the DC pixel in the same location as a
      well formed call to any standard FFT for both even and odd sized input
      arrays. The DC pixel is located at ``np.floor(npix/2) + 1``, which is
      consistent with calls to Numpy's IFFT method where the input and output
      are correctly shifted:
      ``np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(f)))``.

    * If the y-axis shift behavior is not what you are expecting, you most
      likely have your plotting axes flipped (matplotlib's default behavior is
      to place [0,0] in the upper left corner of the axes). This may be resolved
      by either flipping the sign of the y component of ``shift`` or by passing
      ``origin = 'lower'`` to ``imshow()``.

    References
    ----------
    [1] Soummer, et. al. Fast computation of Lyot-style coronagraph propagation (2007)

    [2] `Expressing the inverse DFT in terms of the DFT <https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Expressing_the_inverse_DFT_in_terms_of_the_DFT>`_.

    """
    F = np.asarray(F)
    N = F.size
    # will allocate memory for F if out == None
    F = dft2(np.conj(F), alpha, npix, shift, unitary=unitary, out=out)
    np.conj(F, out=F)
    return np.divide(F, N, out=F)
