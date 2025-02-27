.. _user.fundamentals.wavefront:

.. .. currentmodule:: lentil

*******************
Wavefront and Field
*******************

Lentil's |Wavefront| class represents a discretely sampled electric field.
Wavefronts in Lentil are monochromatic (i.e. they represent a single 
wavelength). A wavefront is defined by the following parameters:

* :attr:`~lentil.Wavefront.wavelength` - The wavefront wavelength
* :attr:`~lentil.Wavefront.pixelscale` - The physical sampling of the 
  wavefront's electric field
* :attr:`~lentil.Wavefront.diameter` - The outscribing diameter around
  the wavefront
* :attr:`~lentil.Wavefront.focal_length` - The focal length of the 
  wavefront. A plane wave has an infinite focal length denoted by ``np.inf``
* :attr:`~lentil.Wavefront.tilt` - A tuple of tilt angles about the x and y 
  axes
* :attr:`~lentil.Wavefront.ptype` - The plane type of the wavefront

Wavefronts are mutable through interaction with a |Plane|. For more details
see :ref:`user.fundamentals.wavefront.plane_wavefront`.

Wavefonts have several useful attributes for working with their electric
field:

* :attr:`~lentil.Wavefront.field` - returns the wavefront's complex field
* :attr:`~lentil.Wavefront.intensity` - returns the intensity (modulus
  squared) of the wavefront's complex field
* :attr:`~lentil.Wavefront.shape` - returns the shape of the wavefront's
  complex field. If shape is ``()``, the wavefront is considered infinite
  and will inherit any finite plane shape 

ptype
=====
A wavefront's type (:attr:`~lentil.Wavefront.ptype`) defines how it interacts 
with a |Plane|. When a wavefront interacts with a plane, it inherits the plane's
``ptype``. Plane type is set automatically and unexpected behavior may
occur if it is changed.

|Wavefront| support the following ptypes:

==================  ===========================================================
:class:`none`       Wavefront has no specific type
:class:`pupil`      Wavefront is at a pupil plane and has a finite focal length
:class:`image`      Wavefront is at an image plane
==================  ===========================================================

The rules defining when a wavefront is allowed to interact with a plane based
on ``ptype`` are described :ref:`here <user.fundamentals.wavefront.ptype_rules>`.

.. _user.fundamentals.wavefront.plane_wavefront:

How a plane affects a wavefront
===============================
An optical plane generally has some effect on a wavefront as it propagates
through the plane. A plane may change a propagating wavefront's amplitude, phase,
and/or physical extent. This |Plane|-|Wavefront| interaction is performed by the
plane's :func:`~lentil.Plane.multiply` method. A |Plane| and |Wavefront| can be
multiplied in two ways:

* By calling :func:`Plane.multiply` directly:

.. code:: pycon

    >>> w1 = plane.multiply(w0)

* By using the built-in multiplication operator (which in turn calls
  :func:`Plane.multiply`):

.. code:: pycon

    >>> w1 = plane * w0

The :func:`~Plane.multiply` method constructs a complex phasor from the plane's
:attr:`~lentil.Plane.amplitude` and :attr:`~lentil.Plane.opd` attributes and the
|Wavefront| wavelength. The plane complex phasor is then multiplied element-wise with
the wavefront's complex data array.

.. _user.fundamentals.wavefront.ptype_rules:

Multiplication rules
====================
The table below outlines when a wavefront can interact with a plane based on their
``ptype`` and what the wavefront's pytpe is after interaction:

+-----------------+------------+-------------+-------------+
|                 | .. centered:: Wavefront ``ptype``      |
+ Plane ``ptype`` +------------+-------------+-------------+
|                 | ``none``   | ``pupil``   | ``image``   |
+=================+============+=============+=============+
| ``none``        | ``none``   | Not allowed | Not allowed |
+-----------------+------------+-------------+-------------+
| ``pupil``       | ``pupil``  | ``pupil``   | Not allowed |
+-----------------+------------+-------------+-------------+
| ``image``       | ``image``  | Not allowed | ``image``   |
+-----------------+------------+-------------+-------------+ 
| ``tilt``        | ``none``   | ``pupil``   | ``image``   |
+-----------------+------------+-------------+-------------+
| ``transform``   | ``none``   | ``pupil``   | ``image``   |
+-----------------+------------+-------------+-------------+

Field
=====
Lentil's wavefront consists of two major components: electric field data and 
metadata about the wavefront. Internally |Wavefront| uses one or more |Field| 
objects to store the electric field data. Users shouldn't need to interact with 
fields directly, but understanding how they are used helps to understand Lentil 
better.

A field is defined by the following parameters:

* :attr:`~lentil.field.Field.data` - The complex eletric field data
* :attr:`~lentil.field.Field.pixelscale` - The physical sampling of the 
  electric field data
* :attr:`~lentil.field.Field.offset` - The location of the center of the 
  field data relative to the global optical axis at (0,0)
* :attr:`~lentil.field.Field.tilt` - A list of |Tilt| objects [1]_


.. [1] Technically any object implementing Lentil's :ref:`TiltInterface <user.advanced.extend.tiltinterface>` is supported