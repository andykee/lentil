Available planes
================


Plane
-----


Pupil
-----
Lentil's :class:`~lentil.Pupil` plane provides a convenient way to represent a 
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

Image
-----
Lentil's :class:`~lentil.Image` plane is used to manipulate a wavefront at an 
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
