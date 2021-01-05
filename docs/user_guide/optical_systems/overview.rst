********
Overview
********

Lentil uses Fourier transform-based algorithms to numerically model the propagation of
a monochromatic wavefront through an optical system represented by a list of planes.
This process is repeated a number of times (once per discrete wavelength being
represented) to simulate the propagation of broadband light through the system. The
exact planes used in a Lentil model depend primarily on the optical design of the system
being represented, but the type of diffraction propagation required to accurately
capture all desired effects also plays a factor.

Lentil plane types
==================
All Lentil planes are derived from its :class:`~lentil.Plane` class. This base class
defines the interface to represent a discretely sampled plane in an optical model. It
can also be used directly in a model. Planes typically have some influence on the
propagation of a wavefront though this is not strictly required and some models may use
*dummy* or *reference* planes as needed.

Lentil provides several general planes that are the building blocks for most optical
models:

* The :class:`~lentil.Pupil` plane provides a convenient way to represent a pupil plane
  in an imaging system. There is nothing particularly special about pupil planes, they
  merely provide a convenient location (mathematically-speaking) to enforce limiting
  apertures or stops and include optical aberrations. More detailed discussion of pupil
  planes is available in [1]_.
* The :class:`~lentil.Image` plane provides a location where the image formed by a
  pupil may be manipulated or viewed.
* The :class:`~lentil.Detector` plane is conceptually identical to the image plane but
  is optimized to very efficiently compute intensity from an image plane complex field.

In addition, several "utility" planes are provided. These planes do not represent
physical components of an optical system, but are used to implement commonly seen 
optical effects:

* The :class:`~lentil.Tilt` is used to represent wavefront tilt about the z-axis in
  terms of radians of x and y tilt.
* The :class:`~lentil.Rotate` plane rotates a wavefront by an arbitrary angle.
* The :class:`~lentil.Flip` plane flips a wavefront about its x, y, or both x and y
  axes.

Finally, a number of specialized optical planes are provided. More discussion of the
types and uses of these planes can be found in :ref:`user-guide.planes`.


``Plane.multiply()``
====================
The :class:`~lentil.Plane`-wavefront interaction is governed by the plane's 
:func:`~lentil.Plane.multiply` method. Generally speaking, this method constructs a
complex phasor from the plane's :attr:`~lentil.Plane.amplitude` and 
:attr:`~lentil.Plane.phase` attributes and the propagation wavelength and performs
an element-wise multiplication with the wavefront's complex amplitude array.


Simple imaging systems
======================
Most imaging systems can be adequately modeled by far-field propagation between a pupil
and image plane. This includes most cameras, telescopes, and imaging instruments. In
these models, all of the optics in a system are represented by a single
:class:`~lentil.Pupil` plane. The results of the diffraction propagation (assuming the
imaging system is operating near focus) can be viewed using either an
:class:`~lentil.Image` or :class:`~lentil.Detector` plane. Because of the optimizations
mentioned above, the :class:`~lentil.Detector` plane should be used to most efficiently
compute image plane intensity.

.. image:: /_static/img/cassegrain.png
    :width: 550px
    :align: center

.. image:: /_static/img/simple_optical_system.png
    :width: 375px
    :align: center

More complicated imaging systems
================================
More complicated imaging systems may contain multiple pupil and image planes. This
includes systems like spectrometers and coronagraphs. With these systems, the
:class:`~lentil.Pupil`, :class:`~lentil.Image`, and :class:`~lentil.Detector` planes are
still used but much more care needs to be taken to ensure each plane is adequately
sampled to avoid the introduction of numerical artifacts in the diffraction propagation.

Fresnel diffraction and intermediate planes
===========================================
If access to an intermediate (non-pupil or image) plane is required or if an imaging 
system is not operating near focus, the near-field (Fresnel) propagation method should 
be used.


.. [1] Goodman, *Introduction to Fourier Optics*.
