.. _user-guide.planes:

.. |Pupil| replace:: :class:`~lentil.Pupil`
.. |Image| replace:: :class:`~lentil.Image`
.. |Detector| replace:: :class:`~lentil.Detector`
.. |Plane| replace:: :class:`~lentil.Plane`

*************************
Specifying Optical Planes
*************************

Plane
=====
Lentil's |Plane| class represents a discretely sampled plane in an optical model. Planes
have attributes for representing the sampled complex amplitude of the plane as well as
additional metadata that may influence how a propagating wavefront interacts with the
plane. A plane is defined by the following parameters:

* :attr:`~lentil.Plane.amplitude` - Defines the relative electric field amplitude
  transmission through the plane
* :attr:`~lentil.Plane.phase` - Defines the electric field phase shift that a wavefront
  experiences when propagating through the plane
* :attr:`~lentil.Plane.mask` - Defines the binary mask over which the plane data is
  valid
* :attr:`~lentil.Plane.segmask` - Defines the set of masks over which any segments in
  the plane are valid. If None, no segments are defined. See :ref:`user-guide.segmented`
  for more information on how to use this attribute.
* :attr:`~lentil.Plane.pixelscale` - Defines the physical sampling of each pixel in
  the above attributes
* :attr:`~lentil.Plane.z` - Defines the plane location along the optical axis of the
  unfolded optical system. This attribute is required when modeling near-field
  diffraction with :func:`~lentil.propagate_fresnel`. It is ignored when modeling far-
  field diffraction with :func:`~lentil.propagate`.

.. note::

    All Plane attributes have sensible default values that have no effect on
    propagations when not defined.


Creating Planes
---------------

Simple Planes are easy to define from the command line:

.. code-block:: pycon

    >>> import lentil
    >>> import matplotlib.pyplot as plt
    >>> p = lentil.Plane(amplitude=lentil.util.circle((256,256), 128))
    >>> plt.imshow(p.amplitude)

.. image:: /_static/img/circle_amplitude.png
    :width: 300px

Once a Plane is defined, its attributes can be modified at any time:

.. code-block:: pycon

    >>> p.phase = 2e-6 * lentil.zernike.zernike(aperture.mask, index=4)
    >>> plt.imshow(p.phase)

.. image:: /_static/img/circle_focus.png
    :width: 300px

Depending on their use, sometimes it will be more convenient to define Planes in a
module. In this case, you should subclass Plane:

.. code-block:: python3

    import lentil

    class CustomPlane(le.Plane):
        def __init__(self):
            self.amplitude = lentil.util.circle((256,256), 128)
            self.opd = 2e-6 * lentil.zernike.zernike(lentil.util.circlemask((256,256),128), 4)

Any of Plane's attributes can also be redefined as properties if further customization
is needed. This is typically necessary if an attribute is stateful or has some sort of
randomness:

.. code-block:: python3

    import lentil

    class CustomPlane(lentil.Plane):
        def __init__(self, focus = 0):
            self.mask = lentil.util.circlemask((256,256), 128)
            self.amplitude = lentil.util.circle((256,256), 128)
            self.focus = focus

        @property
        def phase(self):
            focus_opd = self.focus * lentil.zernike.zernike(self.mask)
            random_opd = lentil.zernike.zernike_compose(self.mask, 1e-6*np.random.random(10))
            return focus_opd + random_opd

.. note::

    Lentil's |Plane| class and its standard library subclasses all use Python's
    ``__init_subclass__()`` method to ensure any required default values are set - even
    if a user-defined subclass does not explicitly call ``Plane``'s constructor
    ``__init__()`` method. For this reason, it is not strictly necessary to call
    ``super().__init__()`` when implementing a custom Plane subclass. It also won't
    hurt, as long as you're careful to either call ``super().__init__()`` before
    defining any static plane attributes or passing these attributes along to the
    ``super().__init__()`` call to ensure they are properly set.

See :ref:`user-guide.custom-planes` for additional information on creating custom planes
and defining special plane behavior.

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
* :attr:`~lentil.Pupil.diameter` - The outscribing diameter of the pupil (in meters)
* :attr:`~lentil.Pupil.pixelscale` - Defines the physical sampling of each pixel in
  the discretely sampled attributes described below

Discretely sampled pupil attributes can also be specified:

* :attr:`~lentil.Pupil.amplitude` - Defines the relative electric field amplitude
  transmission through the pupil
* :attr:`~lentil.Pupil.phase` - Defines the electric field phase shift that a wavefront
  experiences when propagating through the pupil. This term is commonly known as the
  optical path difference (OPD).
* :attr:`~lentil.Pupil.mask` - Defines the binary mask over which the pupil data is
  valid
* :attr:`~lentil.Pupil.segmask` - Defines the set of masks over which any segments in
  the pupil are valid. If None, no segments are defined. See :ref:`user-guide.segmented`
  for more information on how to use this attribute.

.. note::

    All optional Pupil attributes have sensible default values that have no effect on
    propagations when not defined.

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
* :attr:`~lentil.Image.phase` - Defines the electric field phase shift that a wavefront
  experiences when propagating through the image plane.

Detector
========
Lentil's |Detector| plane is used to accumulate the intensity in an image plane. 
Intensity is computed as the absolute value of the complex amplitude in the image plane 
squared.  Similar to the |Image| plane, a detector plane does not have any required parameters 
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
  true. This is becuse the calculation of the real-valued intensity discards the complex 
  field information. Because of this, |Detector| planes can only be used as the final 
  plane in a Lentil model.

Representing Dispersion
=======================
Dispersion is most commonly seen in an optical system as a wavelength-dependent phase
change. In some cases, like with a grating or prism, dispersion is used to achieve some 
desired optical effect. In other cases, dispersion causes an unwanted chromatic 
aberration.

Lentil provides two classes for representing the effects of dispersion: 
:class:`~lentil.DispersivePhase` and :class:`~lentil.DispersiveShift`.

Active Optics and Deformable Mirrors
====================================
Active optics and deformable mirrors are easily represented by defining a phase that
depends on some parameterized state. Because there is no standard architecture for these
types of optical elements, Lentil does not provide a concrete implementation. Instead,
a custom subclass of either |Plane| or |Pupil| should be defined. The exact
implementation details will vary by application, but a simple example of a tip-tilt 
mirror is provided below and additional examples can be found in Model Patterns under
:ref:`patterns.planes`.

.. code-block:: python3

    import lentil
    import numpy as np

    class TipTiltMirror(lentil.Plane):

        def __init__(self):
            self.amplitude = lentil.util.circle((256,256),128)

            self.x = np.zeros(2)

            # Note that we set normalize=False so that each mode spans [-1, 1] and then
            # multiply by 0.5 so that each mode has peak-valley = 1
            self._infl_fn = 0.5 * lentil.zernike.zernike_basis(mask=self.amplitude,
                                                               modes=[2,3],
                                                               normalize=False)

        @property
        def phase(self):
            return np.einsum('ijk,i->jk', self._infl_fn, self.x)

.. code-block:: pycon

    >>> tt = TipTiltMirror()
    >>> tt.x = [1e-6, 3e-6]
    >>> plt.imshow(tt.phase)
    >>> plt.colorbar()

.. image:: /_static/img/circle_tilt.png
    :width: 350px

Tilt
====
The :class:`~lentil.Tilt` plane provides a mechanism for directly specifying wavefront
tilt outside of the context of a discretely sampled :class:`~lentil.Plane` object.
:class:`~lentil.Tilt` is most useful for representing global tilt in an optical system
(for example, due to a pointing error).

Given the following :class:`~lentil.Pupil` and :class:`~lentil.Image` planes:

.. code-block:: pycon

    >>> import lentil
    >>> import matplotlib.pyplot as plt
    >>> pupil = lentil.Pupil(amplitude=lentil.util.circle((256, 256), 128), diameter=1,
    ...                      focal_length=10, pixelscale=1/256)
    >>> detector = lentil.Image(pixelscale=5e-6, shape=(1024, 1024))
    >>> psf = lentil.propagate([pupil, detector], wave=650e-9, npix=(64, 64))
    >>> plt.imshow(psf, origin='lower')

.. image:: /_static/img/psf_64.png
    :width: 300px

it is simple to see the effect of introducing a tilted wavefront into the system:

.. code-block:: pycon

    >>> tilt = lentil.Tilt(x=10e-6, y=-5e-6)
    >>> psf = lentil.propagate([tilt, pupil, detector], wave=650e-9, npix=(64, 64))
    >>> plt.imshow(psf, origin='lower')

.. image:: /_static/img/psf_64_tilt.png
    :width: 300px

Plane Transformations
=====================
The plane transformation examples below are used to modify the following image:

.. code-block:: pycon

    >>> import lentil
    >>> import matplotlib.pyplot as plt
    >>> pupil = lentil.Pupil(amplitude=lentil.util.circle((256, 256), 128), diameter=1,
    ...                      focal_length=10, pixelscale=1/256)
    >>> detector = lentil.Image(pixelscale=5e-6, shape=(1024, 1024))
    >>> psf = lentil.propagate([pupil, detector], wave=650e-9, npix=(128, 128))
    >>> plt.imshow(psf, origin='lower')


.. image:: /_static/img/psf_coma.png
    :width: 300px

Rotate
------
:class:`~lentil.Rotate` can be used to rotate a Wavefront by an arbitrary amount:

.. code-block:: pycon

    >>> rotation = lentil.Rotate(angle=30, unit='degrees')
    >>> psf = lentil.propagate([pupil, rotation, detector], wave=650e-9, npix=(128, 128))
    >>> plt.imshow(psf, origin='lower')

.. image:: /_static/img/psf_coma_rotate.png
    :width: 300px

Flip
----
:class:`~lentil.Flip` can be used to flip a Wavefront about its axes:

.. code-block:: pycon

    >>> flip = lentil.Flip(axis=1)
    >>> psf = lentil.propagate([pupil, flip, detector], wave=650e-9, npix=(128, 128))
    >>> plt.imshow(psf, origin='lower')

.. image:: /_static/img/psf_coma_flip.png
    :width: 300px





.. Lenslet Arrays
.. ==============


