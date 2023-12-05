.. _user.fundamentals.planes:

.. currentmodule:: lentil

.. |Pupil| replace:: :class:`Pupil`
.. |Image| replace:: :class:`Image`
.. |Detector| replace:: :class:`Detector`
.. |Plane| replace:: :class:`Plane`
.. |Wavefront| replace:: :class:`Wavefront`
.. |Tilt| replace:: :class:`Tilt`

******
Planes
******

Lentil plane types
==================
All Lentil planes are derived from the |Plane| class. This base class defines the
interface to represent any discretely sampled plane in an optical model. It can also
be used directly in a model. Planes typically have some influence on the propagation
of a wavefront though this is not strictly required and models may use *dummy* or
*reference* planes as needed.

Lentil provides several general planes that are the building blocks for most optical
models:

* The |Plane| base class provides the core logic for representing and
  working with discretely sampled planes in an optical model.
* The |Pupil| plane provides a convenient way to represent a pupil plane
  in an optical system. There is nothing particularly special about pupil planes, they
  merely provide a convenient location (mathematically-speaking) to enforce limiting
  apertures or stops and include optical aberrations. More detailed discussion of pupil
  planes is available in [1]_.
* The |Image| plane provides a location where the image formed by an
  optical system may be manipulated or viewed.
* The |Detector| plane is conceptually identical to the Image plane but
  is optimized to very efficiently compute intensity from a complex field.

In addition, several "utility" planes are provided. These planes do not represent
physical components of an optical system, but are used to implement commonly encountered
optical effects:

* The :class:`~lentil.Tilt` is used to represent wavefront tilt in terms of radians
  of x and y tilt.
* The :class:`~lentil.Rotate` plane rotates a wavefront by an arbitrary angle.
* The :class:`~lentil.Flip` plane flips a wavefront about its x, y, or both x and y
  axes.

Plane
=====
Lentil's |Plane| class represents a discretely sampled plane in an optical model. Planes
have attributes for representing the sampled complex amplitude of the plane as well as
additional metadata that may influence how a propagating wavefront interacts with the
plane. A plane is defined by the following parameters:

* :attr:`~lentil.Plane.amplitude` - Defines the relative electric field amplitude
  transmission through the plane
* :attr:`~lentil.Plane.opd` - Defines the optical path difference that a wavefront
  experiences when propagating through the plane
* :attr:`~lentil.Plane.mask` - Defines the binary mask over which the plane data is
  valid. If `mask` is 2-dimensional, the plane is assumed to be monolithic. If `mask`
  is 3-dimensional, the plane is assumed to be segmented with the individual segment
  masks inserted along the first dimension. If mask is not provided, it is automatically
  created as needed from the nonzero values in :attr:`~lentil.Plane.amplitude`.

.. plot:: _img/python/segmask.py
    :scale: 50

* :attr:`~lentil.Plane.pixelscale` - Defines the physical sampling of each pixel in
  the above attributes. A simple example of how to calculate the pixelscale for a
  discretely sampled circular aperture is given below:

  .. image:: /_static/img/pixelscale.png
    :width: 450px
    :align: center

.. note::

    All Plane attributes have sensible default values that have no effect on
    propagations when not specified.


Create a new Plane with

.. plot::
    :include-source:
    :scale: 50

    >>> p = lentil.Plane(amplitude=lentil.util.circle((256,256), 120))
    >>> plt.imshow(p.amplitude, origin='lower')

Once a Plane is defined, its attributes can be modified at any time:

.. plot::
    :include-source:
    :scale: 50

    >>> p = lentil.Plane(amplitude=lentil.util.circle((256,256), 120))
    >>> p.opd = 2e-6 * lentil.zernike(p.mask, index=4)
    >>> plt.imshow(p.opd, origin='lower')


Resampling or rescaling a Plane
-------------------------------
It is possible to resample a plane using either the :func:`~lentil.Plane.resample`
or :func:`~lentil.Plane.rescale` methods. Both methods use intrepolation to
resample the amplitude, opd, and mask attributes and readjust the pixelscale
attribute as necessary. The default behavior is to perform this interpolation
on a copy of the plane, but it is possible to operate in-place by setting
``inplace=True``.

.. _user_guide.planes.fit_tilt:

Fitting and removing Plane tilt
-------------------------------
The plane's :func:`~lentil.Plane.fit_tilt` method performs a least squares fit to
estimate and remove tilt from the opd attribute. The tilt removed from the opd
attribute is accounted for by appending an equivalent :class:`~lentil.Tilt` object
to the plane's :attr:`~lentil.Plane.tilt` attribute. The default behavior is to
perform this operation on a copy of the plane, but it is possible to operate
in-place by setting ``inplace=True``.

See :ref:`user.diffraction.tilt` for additional information on when to use
this method.

Segmented optical systems
=========================
Creating a model of a segmented aperture optical system in Lentil doesn't require any
special treatment. The |Plane| object works the same with sparse or
segmented amplitude, opd, and mask attributes as with monolithic ones.

That being said, it is advantageous from a performance point of view to supply a
3-dimensional `segment mask` when specifying a Plane's :attr:`~lentil.Plane.mask`
attribute rather than a flattened 2-dimensional `global mask` when working
with a segmented aperture, as depicted below:

.. plot:: _img/python/segmask.py
    :scale: 50

This modification is not necessary to achieve accurate propagations, but can
greatly improve performance. For additional details, see
:ref:`user.diffraction.segmented`.

.. _user.planes.pupil:

Pupil
=====
Lentil's |Pupil| class provides a convenient way to represent a generalized pupil
function. |Pupil| planes behave exactly like |Plane| objects but introduce an implied
spherical phase term defined by the :attr:`~lentil.Pupil.focal_length` attribute. The
spherical phase term is opaque to the user but is given by

.. math::

    \frac{1}{2f} \left(x^2 + y^2\right)

where :math:`f` is the focal length and :math:`x` and :math:`y` are pupil plane
coordinates.

A pupil is defined by the following required parameters:

* :attr:`~lentil.Pupil.focal_length` - The effective focal length (in meters)
  represented by the pupil
* :attr:`~lentil.Pupil.pixelscale` - Defines the physical sampling of each pixel in
  the discretely sampled attributes described below

Discretely sampled pupil attributes can also be specified:

* :attr:`~lentil.Pupil.amplitude` - Defines the relative electric field amplitude
  transmission through the pupil
* :attr:`~lentil.Pupil.opd` - Defines the optical path difference that a wavefront
  experiences when propagating through the pupil.
* :attr:`~lentil.Pupil.mask` - Defines the binary mask over which the pupil data is
  valid. If `mask` is 2-dimensional, the pupil is assumed to be monolithic. If `mask`
  is 3-dimensional, the pupil is assumed to be segmented with the segment masks
  allocated along the first dimension. If mask is not provided, it is automatically
  created as needed from the nonzero values in :attr:`~lentil.Pupil.amplitude`.

.. note::

    All optional Pupil attributes have sensible default values that have no effect on
    propagations when not defined.

Create a pupil with:

.. code-block:: pycon

    >>> p = lentil.Pupil(focal_length=10, pixelscale=1/100, amplitude=1, opd=0)

Image
=====
Lentil's |Image| plane is used to either manipulate or view a wavefront at a focal point
in an optical system. An image plane does not have any required parameters although any
of the following can be specified:

* :attr:`~lentil.Image.pixelscale` - Defines the physical sampling of each pixel in
  the image plane. If not provided, the sampling will be automatically selected to
  ensure the results are at least Nyquist sampled.
* :attr:`~lentil.Image.shape` - Defines the shape of the image plane. If not provided,
  the image plane will grow as necessary to capture all data.
* :attr:`~lentil.Image.amplitude` - Definers the relative electric field amplitude
  transmission through the image plane.
* :attr:`~lentil.Image.opd` - Defines the optical path difference that a wavefront
  experiences when propagating through the image plane.

Detector
========
Lentil's |Detector| plane is used to accumulate the intensity in an image plane.
Intensity is computed as the absolute value of the complex amplitude in the image plane
squared:

.. math::

    \mathbf{I} = \left|\mathbf{W}\right|^2

Similar to the |Image| plane, a detector plane does not have any required parameters
although any of the following can be specified:

* :attr:`~lentil.Detector.pixelscale` - Defines the physical sampling of each pixel in
  the image plane. If not provided, the sampling will be automatically selected to
  ensure the results are at least Nyquist sampled.
* :attr:`~lentil.Detector.shape` - Defines the shape of the image plane. If not provided,
  the image plane will grow as necessary to capture all data.

While an |Image| plane can be used to compute intensity, the |Detector| plane implements
an algorithm that greatly reduces the memory footprint and increases the speed of this
operation. Details of this algorithm are available in the :ref:`technical-notes`.

.. note::

  An |Image| plane is interchangeable with a |Detector| plane, but the converse is not
  true. This is because the calculation of the real-valued intensity discards the complex
  field information. Because of this, |Detector| planes can only be used as the final
  plane in a Lentil model.

.. _user.planes.tilt:

Tilt
====
The :class:`~lentil.Tilt` plane provides a mechanism for directly specifying wavefront
tilt outside of the context of a discretely sampled |Plane| object. :class:`~lentil.Tilt`
is most useful for representing global tilt in an optical system (for example, due to a
pointing error).

Given the following |Pupil| and |Detector| planes:

.. plot::
    :include-source:
    :scale: 50

    >>> pupil = lentil.Pupil(amplitude=lentil.util.circle((256, 256), 120),
    ...                      focal_length=10, pixelscale=1/250)
    >>> w = lentil.Wavefront(650e-9)
    >>> w *= pupil
    >>> w = lentil.propagate_dft(w, pixelscale=5e-6, shape=(64,64), oversample=2)
    >>> plt.imshow(w.intensity)

It is simple to see the effect of introducing a tilted wavefront into the system:

.. plot::
    :include-source:
    :scale: 50

    >>> pupil = lentil.Pupil(amplitude=lentil.util.circle((256, 256), 120),
    ...                      focal_length=10, pixelscale=1/250)
    >>> tilt = lentil.Tilt(x=10e-6, y=-5e-6)
    >>> w = lentil.Wavefront(650e-9)
    >>> w *= pupil
    >>> w *= tilt
    >>> w = lentil.propagate_dft(w, pixelscale=5e-6, shape=(64,64), oversample=2)
    >>> plt.imshow(w.intensity, origin='lower')

.. note::

  Notice the use of ``origin='lower'`` in the plot above. For an explanation, see
  the note :ref:`here <user.coordinate_system.origin>`.

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

.. _user_guide.planes.special:

Specialized planes
==================

Grism
-----
A grism is a combination of a diffraction grating and a prism that creates a dispersed
spectrum normal to the optical axis. This is in contrast to a single grating or prism,
which creates a dispersed spectrum at some angle that deviates from the optical axis.
Grisms are most commonly used to create dispersed spectra for slitless spectroscopy or
to create interference fringes for dispersed fringe sensing.

Lentil's :class:`~lentil.Grism` plane provides a straightforward mechanism for
efficiently modeling a grism.


Active optics and deformable mirrors
====================================
Active optics and deformable mirrors are easily represented by defining an OPD that
depends on some parameterized state. Because there is no standard architecture for these
types of optical elements, Lentil does not provide a concrete implementation. Instead,
a custom subclass of either |Plane| or |Pupil| should be defined. The exact
implementation details will vary by application, but a simple example of a tip-tilt
mirror where the plane's OPD is computed dynamically based on the state `x` is
provided below. Additional examples can be found in Model Patterns under
:ref:`patterns.planes`.

.. code-block:: python3

    import lentil
    import numpy as np

    class TipTiltMirror(lentil.Plane):

        def __init__(self):
            self.amplitude = lentil.circle((256,256),120)

            self.x = np.zeros(2)

            # Note that we set normalize=False so that each mode spans [-1, 1] and then
            # multiply by 0.5 so that each mode has peak-valley = 1
            self._infl_fn = 0.5 * lentil.zernike_basis(mask=self.amplitude,
                                                       modes=[2,3],
                                                       normalize=False)

        @property
        def opd(self):
            return np.einsum('ijk,i->jk', self._infl_fn, self.x)

.. code-block:: pycon

    >>> tt = TipTiltMirror()
    >>> tt.x = [1e-6, 3e-6]
    >>> plt.imshow(tt.opd)
    >>> plt.colorbar()

.. plot::
    :scale: 50

    import matplotlib.pyplot as plt
    import lentil

    mask = lentil.circlemask((256,256), 120)
    opd = lentil.zernike_compose(mask, [0, 1e-6, 3e-6], normalize=False)

    im = plt.imshow(opd, origin='lower')
    plt.colorbar(im, fraction=0.046, pad=0.04)

Customizing Plane
=================
The Plane class or any of the classes derived from Plane can be subclassed to modify
any of the default behavior. Reasons to do this may include but are not limited to:

* Dynamically computing the :attr:`~lentil.Plane.opd` attribute
* Changing the Plane-Wavefront interaction by redefining the `Plane.multiply()` method
* Modifying the way a Plane is resampled or rescaled

Some general guidance for how to safely subclass Plane is provided below.

.. note::

    Lentil's |Plane| class and its subclasses all use Python's ``__init_subclass__()``
    method to ensure any required default values are set - even if a user-defined
    subclass does not explicitly call ``Plane``'s constructor ``__init__()`` method. For
    this reason, it is not strictly necessary to call ``super().__init__()`` when
    implementing a custom Plane subclass. It also won't hurt, as long as you're careful
    to either call ``super().__init__()`` before defining any static plane attributes or
    passing these attributes along to the ``super().__init__()`` call to ensure they are
    properly set.

Redefining the amplitude, OPD, or mask attributes
---------------------------------------------------
Plane :attr:`~lentil.Plane.amplitude`, :attr:`~lentil.Plane.opd`, and
:attr:`~lentil.Plane.mask` are all defined as properties, but Python allows you to
redefine them as class attributes without issue:

.. code-block:: python3

    import lentil

    class CustomPlane(le.Plane):
        def __init__(self):
            self.amplitude = lentil.circle((256,256), 128)
            self.opd = lentil.zernike(lentil.circlemask((256,256),128), 4)

If more dynamic behavior is required, the property can be redefined. For example, to
return a new random OPD each time the :attr:`~lentil.Plane.opd` attribute is
accessed:

.. code-block:: python3

    import numpy as np
    import lentil

    class CustomPlane(lentil.Plane):
        def __init__(self):
            self.mask = lentil.circlemask((256,256), 128)
            self.amplitude = lentil.circle((256,256), 128)

        @property
        def phaopdse(self):
            return lentil.zernike_compose(self.mask, np.random.random(10))

It is also straightforward to implement a custom :attr:`~lentil.Plane.opd` property to
provide a stateful OPD attribute:

.. code-block:: python3

    import numpy as np
    import lentil

    class CustomPlane(lentil.Plane):
        def __init__(self, x=np.zeros(10)):
            self.mask = lentil.circlemask((256,256), 128)
            self.amplitude = lentil.circle((256,256), 128)
            self.x = x

        @property
        def opd(self):
            return lentil.zernike_compose(self.mask, self.x)

.. note::

    Polychromatic or broadband diffraction propagations access the OPD, amplitude,
    and mask attributes for each propagatioon wavelength. Because these attributes
    remain fixed during a propagation, it is inefficient to repeatedly recompute
    them. To mitigate this, it can be very useful to provide a mechanism for freezing
    these dynamic attributes. There are many ways to do this. One approach is provided
    below:

    .. code-block:: python3

        import copy
        import numpy as np
        import lentil

        class CustomPlane(lentil.Plane):
            def __init__(self):
                self.mask = lentil.circlemask((256,256), 128)
                self.amplitude = lentil.circle((256,256), 128)

            @property
            def opd(self):
                return lentil.zernike_compose(self.mask, np.random.random(10))

            def freeze(self):
                # Return a copy of CustomPlane with the OPD attribute redefined
                # to be a static copy of the OPD when freeze() is called
                out = copy.deepcopy(self)
                out.opd = self.opd.copy()
                return out


Customizing Plane methods
-------------------------
Any of the |Plane| methods can be redefined in a subclass without restriction. Care
should be taken to ensure any redefined methods return data compatible with the
parent method's return type to preserve compatibility within Lentil.


.. Lenslet Arrays
.. ==============


.. [1] Goodman, *Introduction to Fourier Optics*.
