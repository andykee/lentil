*****************************
Caching in ``fourier.dft2()``
*****************************

:func:`~lentil.fourier.dft2` uses the matrix triple product (MTP) approach to 
compute Fourier transforms with variable input and output plane sampling and 
sizes [1]_. 

The Fourier transform of an array :math:`f` is given by

.. math::

    F = E_1 f E_2

where the Fourier kernels :math:`E_1` and :math:`E_2` are constructed as 

.. math::

    E_1 = \exp{(-2\pi i \alpha_y Y V)}

    E_2 = \exp{(-2\pi i \alpha_x X U)}

given vectors of input plane coordinates :math:`X` and :math:`Y`, vectors of output
plane coordinates :math:`U` and :math:`V`, and sampling coefficients :math:`\alpha_x` 
and :math:`\alpha_y`.

Fourier kernel caching
----------------------
If the input and output array sizes and sampling are held constant, :math:`E_1` and 
:math:`E_2` can be computed once and reused providing increased performance for all 
subsequent evaluations. This is achieved by employing Python's ``modeltools.lru_cache`` 
decorator on a method which computes these matrices. The maximum cache size is 32 
entries.

Caching with variable output plane shifts
-----------------------------------------
Repeated calls to :func:`~lentil.fourier.dft2` with a static :attr:`shift` term will
hit the Fourier kernel cache described above. In the event that :attr:`shift` varies
on repeated calls to :func:`~lentil.fourier.dft2` we can still achieve some performance
gains by caching the :math:`X`, :math:`Y`, :math:`U`, and :math:`V` vectors.

Implementation
--------------
The MTP described above is imlemented in Lentil using a two stage LRU caching approach
that caches the calculation of ``X``, ``Y``, ``U``, and ``V`` in all cases and caches
the calculation of ``E1`` and ``E2`` for when ``shift`` does not vary.

.. code:: python

    @functools.lru_cache(maxsize=32)
    def _dft2_coords(m, n, M, N):
        # Y and X are (r,c) coordinates in the (m x n) input plane f
        # V and U are (r,c) coordinates in the (M x N) ourput plane F

        X = np.arange(n) - np.floor(n/2.0)
        Y = np.arange(m) - np.floor(m/2.0)
        U = np.arange(N) - np.floor(N/2.0)
        V = np.arange(M) - np.floor(M/2.0)
    
        return X, Y, U, V


    @functools.lru_cache(maxsize=32)
    def _dft2_matrices(m, n, M, N, ax, ay, sx, sy):
        X, Y, U, V = _dft2_coords(m, n, M, N)
        E1 = np.exp(-2.0 * np.pi * 1j * ay * np.outer(Y-sy, V-sy)).T
        E2 = np.exp(-2.0 * np.pi * 1j * ax * np.outer(X-sx, U-sx))
        return E1, E2


    def dft2(f, alpha, npix=None, shift=(0, 0)):

        ...

        E1, E2 = _dft2_matrices(m, n, M, N, ax, ay, sx, sy)
        F = np.dot(E1.dot(f), E2)

        ...

        return F

Benchmarks
----------


References
----------

.. [1] Jurling et. al., *Techniques for arbitrary sampling in two-dimensional Fourier transforms*.