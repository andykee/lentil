import functools
import numpy as np

from lentil import util

__all__ = ['dft2', 'idft2', 'czt2', 'iczt2']


FFTW_WISDOM_FILENAME = '.lentil_fftw_wisdom'
_fftw_wisdom_file = lambda: os.path.join(os.path.expanduser('~'), FFTW_WISDOM_FILENAME)

def _fftw_save_wisdom():
    import pyfftw, pickle, os
    pickle.dump(pyfftw.export_wisdom(), open(_fftw_wisdom_file(), 'wb'))

try:
    import pyfftw, atexit, pickle, os
    pyfftw.interfaces.cache.enable()
    pyfftw.interfaces.cache.set_keepalive_time(5000)
    if os.path.isfile(_fftw_wisdom_file()):
        pyfftw.import_wisdom(pickle.load(open(_fftw_wisdom_file(), 'rb')))
    atexit.register(_fftw_save_wisdom)
except ImportError:
    pyfftw = None


def dft2(f, alpha, npix=None, shift=(0, 0), offset=(0, 0), unitary=True, out=None):
    """Compute the 2-dimensional discrete Fourier Transform.

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
        Number of pixels in (x,y) to shift the DC pixel in the output plane
        with the origin centrally located in the plane. Default is ``(0,0)``.

    offset : array_like, optional
        Number of pixels in (x,y) that the input plane is shifted relative to
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

    ax, ay = _sanitize_alpha(alpha)

    f = np.asarray(f)
    m, n = f.shape

    if npix is None:
        npix = [m, n]
    M, N = _sanitize_npix(npix)

    if offset is None:
        offset = (0, 0)

    sx, sy = _sanitize_shift(shift)
    ox, oy = _sanitize_npix(offset)

    if out is not None:
        if not np.can_cast(complex, out.dtype):
            raise TypeError(f"Cannot cast complex output to dtype('{out.dtype}')")

    E1, E2 = _dft2_matrices(m, n, M, N, ax, ay, sx, sy, ox, oy)
    F = np.dot(E1.dot(f), E2, out=out)

    # now calculate the answer, without reallocating memory
    if unitary:
        np.multiply(F, np.sqrt(np.abs(ax * ay)), out=F)

    return F


@functools.lru_cache(maxsize=32)
def _dft2_matrices(m, n, M, N, ax, ay, sx, sy, ox, oy):
    X, Y, U, V = _dft2_coords(m, n, M, N)
    E1 = util.expc(-2.0 * np.pi * ay * np.outer(Y-sy+oy, V-sy)).T
    E2 = util.expc(-2.0 * np.pi * ax * np.outer(X-sx+ox, U-sx))
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


def idft2(F, alpha, npix=None, shift=(0,0), unitary=True, out=None):
    """Compute the 2-dimensional inverse discrete Fourier Transform.

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



def czt2(f, alpha, npix=None, shift=(0, 0), unitary=True):
    """Compute the 2-dimensional discrete Fourier Transform using the chirp
    z-transform.

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

    Returns
    -------
    F : complex ndarray

    Notes
    -----
    * If ``pyfftw`` is available it is used to increase FFT performance

    * Setting ``alpha = 1/f.shape`` and ``npix = f.shape`` is equivalent to

      ::

          F = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(f)))

    * ``czt2()`` is designed to place the DC pixel in the same location as a
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
    [1] Rabiner, et. al. The Chirp z-Transform Algorithm. Iransactions on Audio and Electroacoustics. (1969)

    """
    ax, ay = _sanitize_alpha(alpha)

    f = np.asarray(f)
    m, n = f.shape

    if npix is None:
        npix = [m, n]
    M, N = _sanitize_npix(npix)

    sx, sy = _sanitize_shift(shift)

    L = _czt2_optimal_L((m,n), (M,N))

    z_row, beta_row, gamma_row = _czt_coeffs(m, M, L[0], sy, ay, unitary)
    z_col, beta_col, gamma_col = _czt_coeffs(n, N, L[1], sx, ax, unitary)

    f_phased = np.multiply(f, beta_col)
    f_phased = np.multiply(f_phased, beta_row[:, np.newaxis], out=f_phased)

    if pyfftw:
        H = pyfftw.interfaces._utils._Xfftn(f_phased, L, (-2,-1), False, 'FFTW_MEASURE',
                                           48, True, True, 'fft2', normalise_idft=False,
                                           ortho=False)
        H = np.multiply(H, z_col, out=H)
        H = np.multiply(z_row[:, np.newaxis], H, out=H)
        X = pyfftw.interfaces._utils._Xfftn(H, L, (-2,-1), False, 'FFTW_MEASURE', 48,
                                           True, True, 'ifft2', normalise_idft=False,
                                           ortho=False)
    else:
        H = np.fft.fft2(f_phased, L)
        H = np.multiply(H, z_col, out=H)
        H = np.multiply(z_row[:, np.newaxis], H, out=H)
        X = np.fft.ifft2(H)

    X = np.multiply(X[0:M, 0:N], gamma_col)
    X = np.multiply(gamma_row[:, np.newaxis], X, out=X)
    return X


def iczt2(F, alpha, npix=None, shift=(0, 0), unitary=True):
    """Compute the 2-dimensional inverse discrete Fourier Transform using the
    chirp z-transform.

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

    Returns
    -------
    F : complex ndarray

    Notes
    -----
    * If ``pyfftw`` is available it is used to increase FFT performance

    * Setting ``alpha = 1/F.shape`` and ``npix = F.shape`` is equivalent to

      ::

          f = np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(F)))

    * ``czt2()`` is designed to place the DC pixel in the same location as a
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
    [1] Rabiner, et. al. The Chirp z-Transform Algorithm. Iransactions on Audio and Electroacoustics. (1969)

    [2] `Expressing the inverse DFT in terms of the DFT <https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Expressing_the_inverse_DFT_in_terms_of_the_DFT>`_.

    """
    F = np.asarray(F)
    N = F.size
    # will allocate memory for F if out == None
    F = czt2(np.conj(F), alpha, npix, shift, unitary=unitary, out=out)
    np.conj(F, out = F)
    return np.divide(F, F.size,out=F)


@functools.lru_cache(maxsize=32, typed=True)
def _czt_coeffs(N, M, L, shift, alpha, unitary):

    dN = N // 2
    dM = M // 2 + shift

    m_hat = np.linspace(start=-dM, stop=M-dM-1, num=M)

    beta = _czt_beta(N, dN, L, alpha, unitary)

    gamma = util.expc(-np.pi * np.square(m_hat) * alpha)


    z_hat_i = np.zeros(L)
    start = dN - dM

    p = np.linspace(start=start, stop=start+M-1, num=M)
    z_hat_i[0:M] = np.pi * np.square(p)

    p = np.linspace(start=start-N+1, stop=start-1, num=N-1)
    z_hat_i[L-N+1:L] = np.pi * np.square(p)

    z_hat = util.expc(z_hat_i * alpha)
    z_hat[M:L-N+1] = 0
    z = np.fft.fft(z_hat)

    return z, beta, gamma

@functools.lru_cache(maxsize=32, typed=True)
def _czt_beta(N, dN, L, alpha, unitary):
    n_hat = np.linspace(start=-dN, stop=N-dN-1, num=N)

    if unitary:
        scale_factor = np.sqrt(alpha)
    else:
        scale_factor = 1

    if pyfftw:
        scale_factor /= L

    return scale_factor * util.expc(-np.pi * np.square(n_hat) * alpha)

@functools.lru_cache(maxsize=32, typed=True)
def _czt2_optimal_L(N, M):
    divisors = set([2, 3, 5, 7]) # this is what FFTW is fast at
    L = []
    for k in range(2):
        min_L = N[k] + M[k] - 1

        for candidate_L in range(min_L, 2 * min_L + 1):
            current_L = candidate_L
            rem = 0
            while rem == 0:
                for div in divisors:
                    rem = current_L % div
                    if not rem: # div evenly divides current_L
                        current_L = current_L//div
                        if current_L in divisors:
                            L.append(int(candidate_L))
                            rem = 1 # exit the while loop
                        break
            if len(L) == k + 1:
                break

    return tuple(L)


def _sanitize_alpha(x):
    """Return consistent representation of alpha as ax, ay"""
    x = np.asarray(x)
    if x.size == 1:
        ay, ax = float(x), float(x)
    else:
        ay, ax = float(x[0]), float(x[1])
    return ax, ay


def _sanitize_npix(x):
    """Return consistent representation of npix as M, N"""
    x = np.asarray(x)
    if x.size == 1:
        M, N = int(x), int(x)
    else:
        M, N = int(x[0]), int(x[1])
    return M, N


def _sanitize_shift(x):
    """Return consistent representation of shift as sx, sy"""
    if isinstance(x, np.ndarray):
        sx, sy = float(x[0]), float(x[1])
    else:
        sx, sy = x
    return sx, sy
