import functools
import numpy as np
import scipy.fft


def dft2(f, alpha, shape=None, shift=(0, 0), offset=(0, 0), unitary=True, out=None):
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
    shape : int or array_like, optional
        Size of the output array :attr:`F`. If :attr:`shape` is an array,
        ``F.shape = (shape[0], shape[1])``. If :attr:`shape`` is a scalar,
        ``F.shape = (shape, shape)``. Default is ``f.shape``.
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
        A location into which the result is stored. If provided, out.shape ==
        shape and out.dtype == np.complex. If not provided or None, a
        freshly-allocated array is returned.

    Returns
    -------
    F : complex ndarray

    Notes
    -----
    * Setting ``alpha = 1/f.shape`` and ``shape = f.shape`` is equivalent to

      ::

          F = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(f)))

    * ``dft2()`` is designed to place the DC pixel in the same location as a
      well formed call to any standard FFT for both even and odd sized input
      arrays. The DC pixel is located at ``np.floor(shape/2) + 1``, which is
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

    if shape is None:
        shape = [m, n]
    M, N = np.broadcast_to(shape, (2,))

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
    E1 = np.exp(-2.0 * 1j * np.pi * alphar * np.outer(R+offsetr, U-shiftr)).T
    E2 = np.exp(-2.0 * 1j * np.pi * alphac * np.outer(S+offsetc, V-shiftc))
    return E1, E2


@functools.lru_cache(maxsize=32)
def _dft2_coords(m, n, M, N):
    # R and S are (r,c) coordinates in the (m x n) input plane f
    # U and V are (r,c) coordinates in the (M x N) output plane F
    R = np.arange(m) - np.floor(m/2.0)
    S = np.arange(n) - np.floor(n/2.0)
    U = np.arange(M) - np.floor(M/2.0)
    V = np.arange(N) - np.floor(N/2.0)
    return R, S, U, V


def idft2(F, alpha, shape=None, shift=(0,0), unitary=True, out=None):
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
    shape : int or array_like, optional
        Size of the output array :attr:`F`. If :attr:`npshapeix` is an array,
        ``F.shape = (shape[0], shape[1])``. If :attr:`shape` is a scalar,
        ``F.shape = (shape, shape)``. Default is ``F.shape``
    shift : array_like, optional
        Number of pixels in (x,y) to shift the DC pixel in the output plane with
        the origin centrally located in the plane. Default is `[0,0]`.
    unitary : bool, optional
        Normalization flag. If ``True``, the output is scaled by
        ``sqrt(abs(alpha[0]*alpha[1]))`` such that the DFT operation is
        unitary and energy is conserved through the Fourier transform
        operation (Parseval's theorem). In this way, the energy in in a
        limited-area DFT is a fraction of the total energy corresponding to
        the limited area. If ``False``, the output is scaled by
        ``abs(alpha[0]*alpha[1])`` so that ``idft2()`` inverts
        ``dft2(..., unitary=False)``. When ``alpha = 1/F.shape`` this is
        equivalent to the ``1/F.size`` normalization used by
        :func:`numpy.fft.ifft2`. Default is ``True``.

    Returns
    -------
    f : complex ndarray

    Notes
    -----
    * Setting ``alpha = 1/F.shape`` and ``shape = F.shape`` is equivalent to

      ::

          F = np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(F)))

    * ``idft2()`` is designed to place the DC pixel in the same location as a
      well formed call to any standard FFT for both even and odd sized input
      arrays. The DC pixel is located at ``np.floor(shape/2) + 1``, which is
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
    # will allocate memory for F if out == None
    F = dft2(np.conj(F), alpha, shape, shift, unitary=unitary, out=out)
    np.conj(F, out=F)
    if unitary:
        # unitary normalization already handled by dft2
        return F
    return _inverse_normalize(F, alpha)


def _inverse_normalize(F, alpha):
    # non-unitary inverse normalization. This reduces to the conventional
    # 1/F.size when alpha = 1/F.shape but remains correct when F does not span
    # the full period (e.g. a cropped spectrum)
    alpha_row, alpha_col = np.broadcast_to(alpha, (2,))
    return np.multiply(F, np.abs(alpha_row * alpha_col), out=F)


def czt2(f, alpha, shape=None, shift=(0, 0), offset=(0, 0), unitary=True,
         out=None, *, workers=None):
    r"""Compute the 2-dimensional discrete Fourier Transform using the chirp
    Z-transform.

    The DFT is defined in one dimension as

    .. math::

        A_k = \sum_{m=0}^{n-1} a_n \exp \left(-2\pi i\frac{mk}{n}\right) \hspace{3em} k=0, \dots, n-1

    ``czt2()`` is a drop-in replacement for :func:`dft2`: given the same
    arguments, it returns the same result to within floating point precision.
    It evaluates the transform using Bluestein's chirp Z-transform algorithm
    [1]_, [2]_, which replaces the matrix products in :func:`dft2` with
    FFT-based convolutions and is typically faster for large arrays.

    Parameters
    ----------
    f : array_like
        2D array to Fourier Transform
    alpha : float or array_like
        Output plane sampling interval (frequency). If :attr:`alpha` is an
        array, ``alpha[0]`` represents row-wise sampling and ``alpha[1]``
        represents column-wise sampling. If :attr:`alpha` is a scalar,
        ``alpha[0] = alpha[1] = alpha`` gives uniform sampling across the rows
        and columns of the output plane.
    shape : int or array_like, optional
        Size of the output array :attr:`F`. If :attr:`shape` is an array,
        ``F.shape = (shape[0], shape[1])``. If :attr:`shape` is a scalar,
        ``F.shape = (shape, shape)``. Default is ``f.shape``.
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
        A location into which the result is stored. If provided, out.shape ==
        shape and out.dtype == np.complex. If not provided or None, a
        freshly-allocated array is returned.
    workers : int, optional
        Maximum number of workers to use for parallel computation of the
        underlying FFTs. See :func:`scipy.fft.fft` for details. Default is
        ``None``.

    Returns
    -------
    F : complex ndarray

    See Also
    --------
    dft2 : Matrix triple product DFT with the same interface.
    iczt2 : Inverse transform.

    Notes
    -----
    * The DFT kernel is evaluated on the same centered coordinates used by
      :func:`dft2`. In each dimension, with input coordinates
      :math:`x = r - \lfloor m/2 \rfloor + \mathrm{offset}` and output
      coordinates :math:`u = k - \lfloor M/2 \rfloor - \mathrm{shift}`, the
      kernel factors as

      .. math::

          \exp\left(-2\pi i \alpha x u\right) =
          \exp\left(-i\pi\alpha x^2\right)
          \exp\left(-i\pi\alpha u^2\right)
          \exp\left(i\pi\alpha (u-x)^2\right)

      Because :math:`x` and :math:`u` are both unit-spaced, the sum over
      :math:`x` becomes a pre-multiplication by a chirp, a linear
      convolution with a chirp (evaluated with zero-padded FFTs), and a
      post-multiplication by a chirp. Fractional values of ``shift`` and
      ``offset`` are handled exactly.

    * The 2D transform is computed as two separable 1D transforms. The order
      in which the rows and columns are transformed is chosen to minimize
      the total FFT cost.

    * Setting ``alpha = 1/f.shape`` and ``shape = f.shape`` is equivalent to

      ::

          F = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(f)))

    References
    ----------
    .. [1] Bluestein, L. A linear filtering approach to the computation of
       discrete Fourier transform. IEEE Transactions on Audio and
       Electroacoustics 18, 451-455 (1970).

    .. [2] Rabiner, L., Schafer, R., and Rader, C. The chirp z-transform
       algorithm. IEEE Transactions on Audio and Electroacoustics 17, 86-92
       (1969).

    """
    alpha_row, alpha_col = np.broadcast_to(alpha, (2,))

    f = np.asarray(f)
    m, n = f.shape

    if shape is None:
        shape = [m, n]
    M, N = np.broadcast_to(shape, (2,))

    shift_row, shift_col = np.broadcast_to(shift, (2,))
    offset_row, offset_col = np.broadcast_to(offset, (2,))

    if out is not None:
        if not np.can_cast(complex, out.dtype):
            raise TypeError(f"Cannot cast complex output to dtype('{out.dtype}')")

    row_args = (alpha_row, M, shift_row, offset_row, 0)
    col_args = (alpha_col, N, shift_col, offset_col, 1)
    if _czt2_rows_first(m, n, M, N):
        first, second = row_args, col_args
    else:
        first, second = col_args, row_args

    # f is fully consumed by the first pass, so out may alias f
    F = _czt(f, *first, workers=workers)
    F = _czt(F, *second, workers=workers, out=out)

    if unitary:
        np.multiply(F, np.sqrt(np.abs(alpha_row * alpha_col)), out=F)

    return F


def _czt(f, alpha, M, shift, offset, axis, workers=None, out=None):
    # 1D chirp Z-transform of f along axis, evaluating
    #   F[u] = sum_x f[x] exp(-2 pi i alpha x u)
    # with x = R + offset and u = U - shift (see czt2 Notes)
    m = f.shape[axis]

    x = np.arange(m) - np.floor(m/2.0) + offset
    u = np.arange(M) - np.floor(M/2.0) - shift
    # every difference u - x, from smallest to largest
    d = (u[0] - x[-1]) + np.arange(m + M - 1)

    # chirp phases are computed in float64 before exponentiating to preserve
    # precision when pi*alpha*x**2 is large
    bcast = [1, 1]
    bcast[axis] = -1
    pre = np.exp(-1j * np.pi * alpha * x**2).reshape(bcast)
    post = np.exp(-1j * np.pi * alpha * u**2).reshape(bcast)

    L = scipy.fft.next_fast_len(m + M - 1)
    kernel = scipy.fft.fft(np.exp(1j * np.pi * alpha * d**2), L, workers=workers)

    g = scipy.fft.fft(f * pre, L, axis=axis, workers=workers)
    g *= kernel.reshape(bcast)
    g = scipy.fft.ifft(g, axis=axis, workers=workers, overwrite_x=True)

    idx = [slice(None), slice(None)]
    idx[axis] = slice(m - 1, m - 1 + M)
    return np.multiply(g[tuple(idx)], post, out=out)


def _czt2_rows_first(m, n, M, N):
    # Estimate the FFT cost of transforming along rows first vs. columns
    # first. Returns True if rows first is cheaper
    def cost(length, count):
        L = scipy.fft.next_fast_len(length)
        return count * L * np.log2(max(L, 2))

    rows_first = cost(m + M - 1, n) + cost(n + N - 1, M)
    cols_first = cost(n + N - 1, m) + cost(m + M - 1, N)
    return rows_first <= cols_first


def iczt2(F, alpha, shape=None, shift=(0, 0), unitary=True, out=None, *,
          workers=None):
    r"""Compute the 2-dimensional inverse discrete Fourier Transform using the
    chirp Z-transform.

    The IDFT is defined in one dimension as

    .. math::

        a_m = \frac{1}{n}\sum_{k=0}^{n-1} A_k \exp \left(2\pi i\frac{mk}{n}\right) \hspace{3em} m=0, \dots, n-1

    ``iczt2()`` is a drop-in replacement for :func:`idft2`: given the same
    arguments, it returns the same result to within floating point precision.
    See :func:`czt2` for details of the algorithm.

    Parameters
    ----------
    F : array_like
        2D array to Fourier Transform
    alpha : float or array_like
        Input plane sampling interval (frequency). If :attr:`alpha` is an
        array, ``alpha[0]`` represents row-wise sampling and ``alpha[1]``
        represents column-wise sampling. If :attr:`alpha` is a scalar,
        ``alpha[0] = alpha[1] = alpha`` represents uniform sampling across the
        rows and columns of the input plane.
    shape : int or array_like, optional
        Size of the output array :attr:`f`. If :attr:`shape` is an array,
        ``f.shape = (shape[0], shape[1])``. If :attr:`shape` is a scalar,
        ``f.shape = (shape, shape)``. Default is ``F.shape``
    shift : array_like, optional
        Number of pixels in (r,c) to shift the DC pixel in the output plane
        with the origin centrally located in the plane. Default is ``(0,0)``.
    unitary : bool, optional
        Normalization flag. If ``True``, the output is scaled by
        ``sqrt(abs(alpha[0]*alpha[1]))`` such that the DFT operation is
        unitary and energy is conserved through the Fourier transform
        operation (Parseval's theorem). If ``False``, the output is scaled by
        ``abs(alpha[0]*alpha[1])`` so that ``iczt2()`` inverts
        ``czt2(..., unitary=False)``. When ``alpha = 1/F.shape`` this is
        equivalent to the ``1/F.size`` normalization used by
        :func:`numpy.fft.ifft2`. Default is ``True``.
    out : ndarray or None
        A location into which the result is stored. If provided, out.shape ==
        shape and out.dtype == np.complex. If not provided or None, a
        freshly-allocated array is returned.
    workers : int, optional
        Maximum number of workers to use for parallel computation of the
        underlying FFTs. See :func:`scipy.fft.fft` for details. Default is
        ``None``.

    Returns
    -------
    f : complex ndarray

    See Also
    --------
    idft2 : Matrix triple product IDFT with the same interface
    czt2 : Forward chirp Z-transform

    Notes
    -----
    * Setting ``alpha = 1/F.shape`` and ``shape = F.shape`` is equivalent to

      ::

          f = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(F)))

    """
    F = np.asarray(F)
    # will allocate memory for F if out == None
    F = czt2(np.conj(F), alpha, shape, shift, unitary=unitary, out=out,
             workers=workers)
    np.conj(F, out=F)
    if unitary:
        # unitary normalization already handled by czt2
        return F
    return _inverse_normalize(F, alpha)
