.. _user.fundamentals.wavefront:

.. currentmodule:: lentil

*******************
Wavefront and Field
*******************

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

.. _user.fundamentals.plane_wavefront:

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
the wavefront's complex data array:

.. math::

    \mathbf{W_1} = \mathbf{A} \exp\left(\frac{2\pi j}{\lambda} \mathbf{\phi}\right) \circ \mathbf{W_0}


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


Broadcasting
============