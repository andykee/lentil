import numpy as np

__all__ = ['dft2', 'idft2']


def dft2(f, alpha, npix=None, shift=(0, 0), unitary=True):
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
        ay, ax = alpha, alpha
    else:
        ay, ax = alpha

    f = np.asarray(f)
    m, n = f.shape

    if npix is None:
        npix = [m, n]

    npix = np.asarray(npix)
    if npix.size == 1:
        M, N = npix, npix
    else:
        M, N = npix

    sx, sy = np.asarray(shift)

    # Y and X are (r,c) coordinates in the (m x n) input plane f
    # V and U are (r,c) coordinates in the (M x N) ourput plane F

    X = np.arange(n) - np.floor(n/2.0) - sx
    Y = np.arange(m) - np.floor(m/2.0) - sy
    U = np.arange(N) - np.floor(N/2.0) - sx
    V = np.arange(M) - np.floor(M/2.0) - sy

    E1 = np.exp(-2.0 * np.pi * 1j * ay * np.outer(Y, V)).T
    E2 = np.exp(-2.0 * np.pi * 1j * ax * np.outer(X, U))

    F = E1.dot(f).dot(E2)

    if unitary is True:
        return F * np.sqrt(np.abs(ax * ay))
    else:
        return F


def idft2(F, alpha, npix=None, shift=(0, 0), unitary=True):
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
    * `Expressing the inverse DFT in terms of the DFT <https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Expressing_the_inverse_DFT_in_terms_of_the_DFT>`_.  # NOQA

    """
    F = np.asarray(F)
    N = F.size
    return np.conj(dft2(np.conj(F), alpha, npix, shift, unitary))/N
