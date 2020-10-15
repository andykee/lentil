import numpy as np

__all__ = ['dft2', 'idft2']

_cacheE1 = {}
_cacheE2 = {}
_cacheMbyn = {}
_cacheArange = {}

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

    # for profiling 
    # logFile = open('/Users/shanti/Documents/proj/optiix/fourierlog.txt','a')
    # logFile.write("%d %d %d %d %g %g\n" % (M, m, n, N, sx, sy))

    oldWay = False
    # Y and X are (r,c) coordinates in the (m x n) input plane f
    # V and U are (r,c) coordinates in the (M x N) ourput plane F
    if oldWay:
        Y = np.arange(m, dtype = np.float64) - np.floor(m/2.0) - sy
        V = np.arange(M, dtype = np.float64) - np.floor(M/2.0) - sy
        X = np.arange(n, dtype = np.float64) - np.floor(n/2.0) - sx
        U = np.arange(N, dtype = np.float64) - np.floor(N/2.0) - sx
    else:
        # allocate and precompute some of the arrays
        if m in _cacheArange:
            Y = _cacheArange[m]
        else:
            Y = _cacheArange[m] = np.arange(m, dtype = np.float64) - np.floor(m/2.0)
        if M in _cacheArange:
            V = _cacheArange[M]
        else:
            V = _cacheArange[M] = np.arange(M, dtype = np.float64) - np.floor(M/2.0)
        if n in _cacheArange:
            X = _cacheArange[n]
        else:
            X = _cacheArange[n] = np.arange(n, dtype = np.float64) - np.floor(n/2.0)
        if N in _cacheArange:
            U = _cacheArange[N]
        else:
            U = _cacheArange[N] = np.arange(N, dtype = np.float64) - np.floor(N/2.0)
    
    # Can't generally precompute E1 and E2 because the shift term varies in each call, and 
    # storing all possibilities would lead to tremendous memory use. But we can store
    # and reuse the preallocated empty matrices, since those tend to be the same sizes

    if oldWay:
        E1 = np.exp(-2.0 * np.pi * 1j * ay * np.outer(Y, V)).T
    else:
        key = "%d %d" % (M, m)
        if key in _cacheE1:
            E1 = _cacheE1[key]
        else:
            E1 = _cacheE1[key] = np.empty((M, m), dtype = np.complex128)
            # logFile.write("alloc E1 %s\n" % (key))

        np.outer(V - sy, Y - sy, out=E1) # same as transpose(outer(Y,V))
        np.multiply(E1, -2.0 * np.pi * 1j * ay, out=E1)
        np.exp(E1, out=E1)

    if oldWay:
        E2 = np.exp(-2.0 * np.pi * 1j * ax * np.outer(X, U))
    else:
        key = "%d %d" % (n, N)
        if key in _cacheE2:
            E2 = _cacheE2[key]
        else:
            E2 = _cacheE2[key] = np.empty((n, N), dtype = np.complex128)
            # logFile.write("alloc E2 %s\n" % (key))

        np.outer(X - sx, U - sx, out=E2)
        np.multiply(E2, -2.0 * np.pi * 1j * ax, out=E2)
        np.exp(E2, out=E2)

    # place to store the result of E1.dot(f)
    if oldWay:
        F  = E1.dot(f).dot(E2) 
    else:
        key = "%d %d" % (M, n)
        if key in _cacheMbyn:
            Mbyn = _cacheMbyn[key]
        else:
            Mbyn = _cacheMbyn[key] = np.empty((M, n), dtype = np.complex128)
            # logFile.write("alloc Mxn %s\n" % (key))

        np.dot(E1,f,out = Mbyn)
        F = np.dot(Mbyn,E2,out = out) 

    # now calculate the answer, without allocating memory

    if oldWay:
        if unitary is True:
            return F * np.sqrt(np.abs(ax * ay))
        else:
            return F
    else:
        if unitary is True:
            # return F * np.sqrt(np.abs(ax * ay))
            np.multiply(F,np.sqrt(np.abs(ax * ay)),out=F)
    
        # logFile.close()
        return F


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
    * `Expressing the inverse DFT in terms of the DFT <https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Expressing_the_inverse_DFT_in_terms_of_the_DFT>`_.  # NOQA

    """
    F = np.asarray(F)
    N = F.size
    # will allocate memory for F if out == None
    F = dft2(np.conj(F), alpha, npix, shift, unitary=unitary, out=out)
    np.conj(F, out = F)
    return np.divide(F, F.size,out=F)
    # return np.conj(dft2(np.conj(F), alpha, npix, shift, unitary))/N
