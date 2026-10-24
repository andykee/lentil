.. _user.planes:

******
Planes
******

Lentil represents optical systems using one or more planes that have some 
influence on a propagating wavefront. The following planes are the core 
building blocks of most models in Lentil:

* A :ref:`Plane <user.planes.plane>` represents a discretely
  sampled optical plane at an arbitrary location in an optical system.
* A :ref:`Pupil <user.planes.pupil>` represents a discretely
  sampled pupil plane.
* An :ref:`Image <user.planes.image>` represents a discretely
  sampled image plane.

In addition, several "utility" planes are provided. These planes don't 
represent physical components of an optical system, but are used to implement 
commonly encountered optical effects:

* The :ref:`Tilt <user.planes.tilt>` plane is used to represent 
  wavefront tilt in terms x and y tilt (in radians).

.. * The :class:`~lentil.Rotate` plane rotates a wavefront by an arbitrary angle.
.. * The :class:`~lentil.Flip` plane flips a wavefront about its x, y, or both x 
..   and y axes.

.. _user.planes.plane:

Plane basics
============
Planes are described by a common set of attributes. The most commonly used 
attributes are

* ``amplitude`` - defines the electric field amplitude transmission through 
  the plane
* ``opd`` - defines the optical path difference that a wavefront experiences 
  when propagating through the plane
* ``mask`` - defines the binary mask over which the plane data is valid. If 
  ``mask`` is 2-dimensional, the plane is assumed to be monolithic. If 
  ``mask`` is 3-dimensional, the plane is assumed to be segmented with the 
  individual segment masks provided along the first dimension. If mask is not 
  provided, it is automatically created from the nonzero values in 
  ``amplitude``.

.. plot:: _img/python/segmask.py
    :scale: 50

* ``pixelscale`` - defines the physical sampling of the above attributes. A 
  simple example of how to calculate the pixelscale for a discretely sampled 
  circular aperture is given below:

  .. image:: /_static/img/pixelscale.png
    :width: 450px
    :align: center

* ``ptype`` - specifies the plane type. Additional details are available in 
  the sections describing :ref:`user.wavefront.plane_wavefront` 
  and :ref:`user.wavefront.ptype_rules`

  Valid types are:
  
  ============= ============================
  ptype         Planes with this type
  ============= ============================
  ``none``      ``Plane``
  ``pupil``     ``Pupil``
  ``image``     ``Image``
  ``tilt``      ``Tilt``, ``DispersiveTilt``
  ``transform`` ``Rotate``, ``Flip``
  ============= ============================
  
.. note::

    All Plane attributes have sensible default values that have no effect on
    propagations when not specified.

Plane creation
--------------
Create a new plane with

.. plot::
    :include-source:
    :scale: 50

    >>> p = lentil.Plane(amplitude=lentil.circle((256,256), 120))
    >>> plt.imshow(p.amplitude, cmap='gray')

Once a plane is defined, its attributes can be modified at any time:

.. code-block:: pycon

    >>> p = lentil.Plane(amplitude=lentil.circle((256,256), 120))
    >>> p.opd = 2e-6 * lentil.zernike(p.mask, index=4)
    >>> plt.imshow(p.opd, cmap=opd_cmap)

.. plot::
    :scale: 50

    >>> import matplotlib
    >>> p = lentil.Plane(amplitude=lentil.circle((256,256), 120))
    >>> p.opd = 2e-6 * lentil.zernike(p.mask, index=4)
    >>> opd = p.opd
    >>> opd[np.where(opd == 0)] = np.nan
    >>> plt.imshow(opd, cmap=opd_cmap)


.. Basic arithmetic operations
.. ---------------------------
.. Planes with the same ``ptype`` can be added or subtracted, returning a new
.. ``plane``:

.. In-place operations are also supported:


Resampling and rescaling
------------------------
It is possible to resample a plane using either the 
:func:`~lentil.Plane.resample` or :func:`~lentil.Plane.rescale` methods. Both 
methods use intrepolation to resample the ``amplitude``, ``opd``, and ``mask`` 
attributes and readjust the ``pixelscale`` attribute as necessary.


.. _user.planes.pupil:

Pupil
=====
Lentil's ``Pupil`` class provides a convenient way to represent a 
generalized pupil function. Pupil planes behave exactly like plane 
objects but introduce an implied spherical phase term defined by the 
``focal_length`` attribute. The spherical phase term is opaque to the 
user but is given by

.. math::

    \frac{1}{2f} \left(x^2 + y^2\right)

where :math:`f` is the focal length and :math:`x` and :math:`y` are pupil 
plane coordinates.

A pupil is defined by the following required parameters:

* ``focal_length`` - The effective focal length (in meters)
  represented by the pupil
* ``pixelscale`` - Defines the physical sampling of each 
  pixel in the discretely sampled attributes described below

Discreetly sampled pupil attributes can also be specified:

* ``ampltiude`` - Defines the relative electric field amplitude transmission 
  through the pupil
* ``opd`` - Defines the optical path difference that a wavefront experiences 
  when propagating through the pupil.
* ``mask`` - Defines the binary mask over which the pupil data is valid. If 
  ``mask`` is 2-dimensional, the pupil is assumed to be monolithic. If ``mask`` 
  is 3-dimensional, the pupil is assumed to be segmented with the segment 
  masks allocated along the first dimension. If mask is not provided, it is 
  automatically created as needed from the nonzero values in ``amplitude``.

.. note::

    All optional Pupil attributes have sensible default values that have no 
    effect on propagations when not defined.

Create a pupil with:

.. code-block:: pycon

    >>> p = lentil.Pupil(focal_length=10, pixelscale=1/100, amplitude=1, opd=0)

.. _user.planes.image:

Image
=====
Lentil's ``Image`` plane is used to either manipulate or view a wavefront at an 
image plane in an optical system. An image plane does not have any required 
parameters although any of the following can be specified:

* ``pixelscale`` - Defines the physical sampling of each pixel in the image 
  plane. If not provided, the sampling will be automatically selected to ensure 
  the results are at least Nyquist sampled.
* ``shape`` - Defines the shape of the image plane. If not provided, the image 
  plane will grow as necessary to capture all data.
* ``amplitude`` - Definers the relative electric field amplitude transmission 
  through the image plane.
* ``opd`` - Defines the optical path difference that a wavefront experiences 
  when propagating through the image plane.

.. _user.planes.tilt:

Tilt
====
The ``Tilt`` plane provides a mechanism for representing wavefront
tilt in terms of angular rotations about the plane's x and y-axes. This 
representation is separate and in addition to any tilt specified in a 
plane's ``opd`` attribute. Tilt planes most useful for representing global 
tilt in an optical system (for example, due to a pointing error).

Given the following pupil plane:

.. plot::
    :include-source:
    :scale: 50

    >>> pupil = lentil.Pupil(amplitude=lentil.circle((256, 256), 120),
    ...                      focal_length=10, pixelscale=1/250)
    >>> w = lentil.Wavefront(650e-9)
    >>> w *= pupil
    >>> w = lentil.propagate_dft(w, pixelscale=5e-6, shape=(64,64), oversample=2)
    >>> plt.imshow(w.intensity, cmap='inferno')

It is simple to see the effect of introducing a tilted wavefront into the 
system:

.. plot::
    :include-source:
    :scale: 50

    >>> pupil = lentil.Pupil(amplitude=lentil.circle((256, 256), 120),
    ...                      focal_length=10, pixelscale=1/250)
    >>> tilt = lentil.Tilt(x=10e-6, y=-5e-6)
    >>> w = lentil.Wavefront(650e-9)
    >>> w *= pupil
    >>> w *= tilt
    >>> w = lentil.propagate_dft(w, pixelscale=5e-6, shape=(64,64), oversample=2)
    >>> plt.imshow(w.intensity, origin='lower', cmap='inferno')

.. note::

  Notice the use of ``origin='lower'`` in the plot above. For an explanation, 
  see the note :ref:`here <user.coordinate_system.origin>`.

.. .. _user_guide.planes.transformations:

.. Plane transformations
.. =====================
.. The plane transformation examples below are used to modify the following image:

.. .. code-block:: pycon
..
..     >>> pupil = lentil.Pupil(amplitude=lentil.util.circle((256, 256), 128),
..     ...                      focal_length=10, pixelscale=1/256)
..     >>> detector = lentil.Detector(pixelscale=5e-6, shape=(1024, 1024))
..     >>> psf = lentil.propagate([pupil, detector], wave=650e-9, npix=(128, 128))
..     >>> plt.imshow(psf, origin='lower')


.. .. image:: /_static/img/psf_coma.png
..     :width: 300px

.. Rotate
.. ------
.. :class:`~lentil.Rotate` can be used to rotate a Wavefront by an arbitrary amount:

.. .. code-block:: pycon

..     >>> rotation = lentil.Rotate(angle=30, unit='degrees')
..     >>> psf = lentil.propagate([pupil, rotation, detector], wave=650e-9, npix=(128, 128))
..     >>> plt.imshow(psf, origin='lower')

.. .. image:: /_static/img/psf_coma_rotate.png
..     :width: 300px

.. Flip
.. ----
.. :class:`~lentil.Flip` can be used to flip a Wavefront about its axes:

.. .. code-block:: pycon

..     >>> flip = lentil.Flip(axis=1)
..     >>> psf = lentil.propagate([pupil, flip, detector], wave=650e-9, npix=(128, 128))
..     >>> plt.imshow(psf, origin='lower')

.. .. image:: /_static/img/psf_coma_flip.png
..     :width: 300px

DispersiveTilt
==============
The :class:`~lentil.DispersiveTilt` class can be used to represent any
wavelength-depentent tilt. Examples include modeling chromatic aberrations, 
prisms, or diffraction gratings, just to name a few. With a dispersive tilt,
light is dispersed along a path called a spectral trace. The location of a 
specific wavelength along this trace is determined by the dispersion function.

The basic geometry of spectral dispersion is illustrated in the figure below:

.. image:: /_static/img/grism_geometry.png
  :align: center
  :width: 400px

The spectral trace is parameterized by a polynomial of the form

.. math::

  y = a_n x^n + \cdots + a_2 x^2 + a_1 x + a_0

and should return units of meters on the focal plane provided an input
in meters on the focal plane.

Similarly, the wavelength along the trace is parameterized by a
polynomial of the form

.. math::

  \lambda = a_n d^n + \cdots + a_2 d^2 + a_1 d + a_0

and should return units of meters of wavelength provided an input distance
along the spectral trace.

Lentil supports trace and dispersion functions with any arbitrary polynomial
order. While a simple analytic solution exists for modeling first-order trace
and/or dispersion, there is no general solution for higher order functions.

As a result, trace and/or dispersion polynomials with order > 1 are evaluated
numerically. Although the effects are small, this approach impacts both the
speed and precision of modeling grisms with higher order trace and/or
dispersion functions. In cases where speed or accuracy are extremely important,
a custom solution may be required.


Grism
-----
.. warning::

    :class:`~lentil.Grism` is deprecated and will be removed in a future 
    version. Use :class:`~lentil.DispersiveTilt` instead.

A grism is a combination of a diffraction grating and a prism that creates a 
dispersed spectrum normal to the optical axis. This is in contrast to a single 
grating or prism, which creates a dispersed spectrum at some angle that 
deviates from the optical axis. Grisms are most commonly used to create 
dispersed spectra for slitless spectroscopy or to create interference fringes 
for dispersed fringe sensing.

Lentil's :class:`~lentil.Grism` plane provides a straightforward mechanism for
efficiently modeling a grism.


.. Lenslet Arrays
.. ==============

