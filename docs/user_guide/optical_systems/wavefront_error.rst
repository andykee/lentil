.. currentmodule:: lentil

*****************************
Representing Wavefront Errors
*****************************

Wavefront error is represented in a :class:`Plane` by specifying its 
:attr:`~Plane.phase` attribute. For static errors, a simple wavefront error map is
sufficient. For more complicated errors that are random or time-varying in nature, a 
more dynamic and/or state-based approach is required. 

Any Python type that is array-like can be provided and will be converted by 
Lentil to a Numpy array. 

.. note::

    For :class:`Pupil` planes, the :attr:`~Pupil.phase` attribute represents the optical
    path difference (OPD) relative to the pupil's reference sphere. 


Zernike Polynomials
===================
Lentil's :ref:`zernike module <api-zernike>` contains methods for working with `Zernike
polynomials <https://en.wikipedia.org/wiki/Zernike_polynomials>`_. Methods are provided 
for creating, combining, fitting, and removing Zernikes.

.. note::

    Lentil uses the Noll indexing scheme for defining Zernike polynomials [1]_.

Wavefront error maps are easily computed using either the :func:`zernike.zernike` or 
:func:`zernike_compose` functions. For example, we can represent 100 nm of focus over a 
circular aperture with :func:`zernike.zernike`:

.. code-block:: pycon

    >>> import matplotlib.pyplot as plt
    >>> import lentil as le
    >>> mask = le.util.circlemask((256,256), 128)
    >>> z4 = 100e-9 * le.zernike.zernike(mask, index=4)
    >>> plt.imshow(z4)

.. image:: /_static/img/circle_focus.png
    :width: 300px

Any combination of Zernike polynomials can be combined by providing a list of coefficients
to the :func:`zernike.zernike_compose` function. For example, we can represent 200 nm of 
focus and -100 nm of astigmatism as:

.. code-block:: pycon

    >>> import matplotlib.pyplot as plt
    >>> import lentil as le
    >>> mask = le.util.circlemask((256,256), 128)
    >>> coefficients = [0, 0, 0, 200e-9, 0, -100e-9]
    >>> z = le.zernike.zernike_compose(mask, coefficients)
    >>> plt.imshow(z)

.. image:: /_static/img/api/zernike/zernike_compose_2.png
    :width: 300px

Note that the coefficients list is ordered according to the Noll indexing scheme so the
first entry in the list represents piston, the second represents, tilt, and so on.

For models requiring many random trials, it may make more sense to pre-compute the 
Zernike modes once and accumulate the error map for each new state. We can do this by 
creating a vectorized basis set using :func:`zernike.zernike_basis` and accumulating
each independent term using Numpy's `einsum <https://numpy.org/doc/stable/reference/generated/numpy.einsum.html>`_ 
function.

Note that in this case we are only computing the Zernike modes we intend to use (Noll 
indices 4 and 6) so now the first entry in ``coefficients`` corresponds to focus and the
second corresponds to astigmatism.

.. code-block:: pycon

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import lentil as le
    >>> mask = le.util.circlemask((256,256), 128)
    >>> coefficients = [200e-9, -100e-9]
    >>> basis = le.zernike.zernike_basis(mask, modes=(4,6))
    >>> z = np.einsum('ijk,i->jk', basis, coefficients)
    >>> plt.imshow(z)

.. image:: /_static/img/api/zernike/zernike_compose_2.png
    :width: 300px

If you don't love ``einsum``, it's possible to achieve the same result with Numpy's 
`tensordot <https://numpy.org/doc/stable/reference/generated/numpy.tensordot.html>`_:

.. code-block:: pycon

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import lentil as le
    >>> mask = le.util.circlemask((256,256), 128)
    >>> coefficients = [200e-9, -100e-9]
    >>> basis = le.zernike.zernike_basis(mask, modes=(4,6))
    >>> z = np.tensordot(basis, coefficients, axes=(0,0))
    >>> plt.imshow(z)

.. image:: /_static/img/api/zernike/zernike_compose_2.png
    :width: 300px

Normalization
-------------
Each of Lentil's Zernike functions accepts a ``normalize`` parameter. If ``normalize``
is flase (the default), the raw Zernike mode is returned. Each mode will approximately
span [-1 1] although this shouldn't be relied upon because of the discrete sampling of
the result. If ``normalize`` is true, the Zernike mode will be normalized so that its 
standard deviation equals 1. 

Normalization becomes important when trying to achieve a specific error magnitude, 
whether it be in terms of RMS or peak to valley. To acihieve a specific error in terms
of RMS, Zernike modes should be computed with ``normalize=True`` before multiplying by
the error magnitude:

.. code-block:: pycon

    >>> import lentil as le
    >>> import numpy as np
    >>> mask = le.util.circlemask((256,256), 128)
    >>> z4 = 100e-9 * le.zernike.zernike(mask, mode=4, normalize=True)
    >>> np.std(z4[np.nonzero(z4)])

    9.986295346152438e-08

To achieve a specific error in terms of peak to valley, Zernike modes should be computed
and normalized separately. The separate normalization step should be performed to ensure
the discretely sampled mode spans [-0.5 0.5] before multiplying by the error magnitude:

.. code-block:: pycon

    >>> import lentil as le
    >>> import numpy as np
    >>> mask = le.util.circlemask((256,256), 128)
    >>> z4 = le.zernike.zernike(mask, mode=4)
    >>> z4 /= np.max(z4) - np.min(z4)
    >>> z4 *= 100e-9
    >>> np.max(z4) - np.min(z4)

    1e-07

Defining custom coordinates
---------------------------
By default, all of Lentil's Zernike functions place the center of the coordinate system
at the centroid of the supplied mask with its axes aligned with Lentil's 
:ref:`user-guide.coordinate-system`. This works as expected for the vast majority of 
needs, but in some cases it may be desirable to manually define the coordinate system. 
This is accomplished by using :func:`zernike.zernike_coordinates` to compute ``rho`` and
``theta``, and providing these definitions to the appropriate Zernike function. For 
example, if we have an off-centered sub-aperture but wish to compute focus relative to 
the center of the defined array:

.. code-block:: pycon

    >>> import matplotlib.pyplot as plt
    >>> import lentil as le
    >>> mask = le.util.circlemask((256,256), radius=50, shift=(0,60))
    >>> rho, theta = le.zernike.zernike_coordinates(mask, shift=(0,0))
    >>> z4 = le.zernike.zernike(mask, 4, rho=rho, theta=theta)
    >>> plt.imshow(z4)

If we wish to align a tilt mode with one side of a hexagon:

.. code-block:: pycon

    >>> import matplotlib.pyplot as plt
    >>> import lentil as le
    >>> mask = le.util.hexagon((256,256), radius=128)
    >>> rho, theta = le.zernike.zernike_coordinates(mask, shift=(0,0), rotate=60)
    >>> z2 = le.zernike.zernike(mask, 2, rho=rho, theta=theta)
    >>> plt.imshow(z2)

Sensitivity Matrices
====================
The effects of optical element rigid body perturbations and surface figure errors in the 
exit pupil of an optical system are commonly captured using linear sensitivity matrices. 
These linearized models can be used in place of a full ray-tracing model for representing
small perturbations and errors. In general, a linear wavefront error model has the form:

.. math::

    \mathbf{w} = \mathbf{S}\Delta\mathbf{x}

where :math:`\textbf{w}` is the wavefront error map, :math:`S` is the sensitivity 
matrix, and :math:`\Delta\mathbf{x}` is a vector of perturbations relative to the system 
state about which linearization occurred. 

The :math:`\mathbf{S}` matrix will have either two or three dimensions. For a three-
dimensional sensitivity matrix, the wavefront error map is computed by multiplying 
:math:`\mathbf{S}`  by the :math:`\Delta\mathbf{x}` vector and summing along the first 
dimension:

.. code-block:: pycon

    >>> w = np.einsum('ijk,i->jk', S, dx)

For a two-dimensional sensitivity matrix, each mode is assumed to have been unraveled 
into a vector. The wavefront error is computed by taking the dot product of 
:math:`\mathbf{S}` and :math:`\Delta\mathbf{x}` and reshaping the resulting vector into a
two-dimensional error map. For a sensitivity matrix representing a 256 x 256 pixel 
wavefront map:

.. code-block:: pycon

    >>> w = np.dot(S, dx)
    >>> w.reshape((256,256))


.. Chromatic Aberrations
.. =====================
.. Chromatic aberrations are wavelength-dependent errors cause by dispersion. These 
.. aberrations can be further classified as either transverse or longitudinal. Transverse
.. chromatic aberration causes a wavelength-dependent focus shift and can be implemented
.. by customizing :class:`~lentil.DispersivePhase`'s :func:`~lentil.DispersivePhase.multiply` 
.. method. For example, if an 
.. optical system produces best focus at 550 nm and each nm of wavelength change causes 1 nm
.. of focus error, we represent the wavelength-dependent focus by:

.. .. math::

..     \mbox{Focus shift} = \lambda - 550 \times 10^{-9}

.. We implement this focus shift as an additional focus :attr:`~lentil.Plane.phase` term 
.. that is applied within the plane's :func:`~lentil.Plane._phasor` method. Note that 
.. wavelength is given in :attr:`Wavefront.wavelength`

.. .. code-block:: python3

..     import lentil as le

..     class TransverseCA(le.Plane):

..         def __init__(self, *args, **kwargs):
..             super().__init__(*args, **kwargs)

..             # Pre-compute defocus map for efficiency
..            self.defocus = le.zernike.zernike(mask=self.amplitude, 
..                                               index=4,
..                                               normalize=True)

..         def _phasor(amplitude, ):


.. Transverse chromatic aberration causes a wavelength-dependent magnification across the
.. field. 


.. Atmospheric Turbulence
.. ======================



.. Static Errors
.. =============

.. File-based errors
.. -----------------

.. Parametric errors
.. -----------------

.. Parameterized Timeseries Errors
.. ===============================


.. [1] Noll, RJ. Zernike polynomials and atmospheric turbulence. J Opt Soc Am 66, 207-211  (1976).
