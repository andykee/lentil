.. _user.wavefront:

.. .. currentmodule:: lentil

*******************
Wavefront and Field
*******************

Lentil's |Wavefront| class represents a discretely sampled electric field.
Wavefronts in Lentil are monochromatic (i.e. they represent a single 
wavelength). A wavefront is defined by the following parameters:

* ``wavelength`` - The wavefront wavelength
* ``pixelscale`` - The physical sampling of the wavefront's electric field
* ``diameter`` - The outscribing diameter around the wavefront
* ``focal_length`` - The focal length of the wavefront. A plane wave has an 
  infinite focal length denoted by ``None``
* ``tilt`` - A tuple of tilt angles about the x and y axes
* ``ptype`` - The plane type of the wavefront

Wavefronts are mutable through interaction with a |Plane|. For more details
see :ref:`user.wavefront.plane_wavefront`.

Wavefonts have several useful attributes for working with their electric
field:

* ``field`` - returns the wavefront's complex field
* ``intensity`` - returns the intensity (modulus
  squared) of the wavefront's complex field
* ``shape`` - returns the shape of the wavefront's
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
on ``ptype`` are described :ref:`here <user.wavefront.ptype_rules>`.

.. _user.wavefront.plane_wavefront:

How a plane affects a wavefront
===============================
An optical plane generally has some effect on a wavefront as it propagates
through the plane. A plane may change a propagating wavefront's amplitude, 
phase, and/or physical extent. This ``Plane``-``Wavefront`` interaction is 
implemented using multiplication:

.. code:: pycon

    >>> w1 = plane * w0

When a plane and wavefront are multiplied, a complex phasor is constructed 
from the plane's ``amplitude`` and ``opd`` attributes and the wavefront's 
wavelength. The plane complex phasor is then multiplied element-wise with
the wavefront's complex data array. Internally, the operation looks something
like the following pseudocode:

.. code:: python

   phs = plane.amplitude * np.exp(2* np.pi*1j * plane.opd / wavefront.wavelength)
   wavefront.field = wavefront.field * phasor

.. _user.wavefront.ptype_rules:

Plane-wavefront multiplication rules
====================================
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
metadata about the wavefront. Internally a wavefront uses one or more |Field| 
objects to store the electric field data. Users shouldn't need to interact with 
fields directly, but understanding how they are used helps to understand Lentil 
better.

A field is defined by the following parameters:

* ``data`` - The complex eletric field data
* ``pixelscale`` - The physical sampling of the electric field data
* ``offset`` - The location of the center of the field data relative to the 
  global optical axis at (0,0)
* ``tilt`` - A list of |Tilt| objects [1]_

.. [1] All objects implementing Lentil's :ref:`TiltInterface 
       <user.extend.tiltinterface>` are supported
