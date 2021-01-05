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

Defining custom coordinates
---------------------------
By default, all of Lentil's Zernike functions place the center of the coordinate system
at the centroid of the supplied mask. This works as expected for the vast majority of 
needs, but in some cases it may be desirable to manually define the coordinate system. 
This is accomplished by using :func:`zernike.zernike_coordinates` to compute ``rho`` and
``theta``, and providing these definitions to the appropriate Zernike function.

Linear Influence Functions
==========================

Static Errors
=============

.. File-based errors
.. -----------------

.. Parametric errors
.. -----------------


Atmospheric Turbulence
======================


Parameterized Timeseries Errors
===============================


.. [1] Noll, RJ. Zernike polynomials and atmospheric turbulence. J Opt Soc Am 66, 207-211  (1976).
