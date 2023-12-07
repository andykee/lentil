.. _user.fundamentals.wavefront:

.. currentmodule:: lentil

*********************
Wavefronts and Fields
*********************

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

Broadcasting
============