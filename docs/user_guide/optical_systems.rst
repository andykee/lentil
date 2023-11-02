.. _user_guide.optical_systems:

.. currentmodule:: lentil

**************************
Describing optical systems
**************************

.. |Pupil| replace:: :class:`Pupil`
.. |Image| replace:: :class:`Image`
.. |Detector| replace:: :class:`Detector`
.. |Plane| replace:: :class:`Plane`
.. |Wavefront| replace:: :class:`Wavefront`
.. |Tilt| replace:: :class:`Tilt`

It is common in physical optics modeling to represent an optical system in its
"unfolded" state where all optical elements are arranged along a straight line (the
optical axis). This approach assumes the following:

* The beam is `paraxial <https://en.wikipedia.org/wiki/Paraxial_approximation>`_
* Powered mirrors are represented as equivalent ideal thin lenses
* The beam does not change dimensions across a lens
* A lens has no thickness so all phase changes occur in a plane

Lentil uses |Plane| objects to represent discretely sampled planes in an optical system
and |Wavefront| objects to represent discretely sampled electromagnetic fields as they
propagate through an optical system.

.. _user_guide.optical_systems.plane_wavefront:

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
:attr:`~lentil.Plane.amplitude` and :attr:`~lentil.Plane.phase` attributes and the
|Wavefront| wavelength. The plane complex phasor is then multiplied element-wise with
the wavefront's complex data array:

.. math::

    \mathbf{W_1} = \mathbf{A} \exp\left(\frac{2\pi j}{\lambda} \mathbf{\theta}\right) \circ \mathbf{W_0}


.. If the |Plane| :attr:`~lentil.Plane.tilt` attribute is not empty, its contents are appended
.. to the |Wavefront|. See :ref:`user_guide.planes.fit_tilt` and :ref:`user_guide.diffraction.tilt`
.. for additional details.

Planes in a simple optical system
=================================
Most optical systems can be adequately modeled by a single far-field propagation
between a :ref:`user_guide.planes.pupil` and image plane. This includes most cameras,
telescopes, and imaging instruments. In these models, all of the optics in a system are
represented by a single |Pupil| plane:

.. plot::
    :scale: 50
    :include-source:

    >>> import matplotlib.pyplot as plt
    >>> import lentil
    >>> amplitude = lentil.circle(shape=(256, 256), radius=120)
    >>> opd = lentil.zernike_compose(mask=amplitude,
    ...                              coeffs=[0, 0, 0, 100e-9, 300e-9, 0, -100e-9])
    >>> pupil = lentil.Pupil(amplitude=amplitude, phase=opd, focal_length=10,
    ...                      pixelscale=1/240)
    >>> plt.imshow(pupil.phase, origin='lower')

Segmented optical systems
=========================
Creating a model of a segmented aperture optical system in Lentil doesn't require any
special treatment. The |Plane| and |Pupil| objects work the same with sparse or
segmented amplitude, phase, and mask attributes as with monolithic ones.

That being said, it is advantageous from a performance point of view to supply a
3-dimensional `segment mask` when specifying a Plane's :attr:`~lentil.Plane.mask`
attribute rather than a flattened 2-dimensional `global mask` when working
with a segmented aperture, as depicted below:

.. plot:: _img/python/segmask.py
    :scale: 50

This modification is not necessary to achieve accurate propagations, but can
greatly improve performance. For additional details, see
:ref:`user_guide.diffraction.segmented`.


More complicated optical systems
================================
More complicated imaging systems may contain multiple pupil and image planes. This
includes systems like spectrometers and coronagraphs. With these systems, the
:class:`~lentil.Pupil`, :class:`~lentil.Image`, and :class:`~lentil.Detector` planes are
still used but much more care needs to be taken to ensure each plane is adequately
sampled to avoid the introduction of numerical artifacts in the diffraction propagation.

If access to an intermediate (non-pupil or image) plane is required or if an imaging
system is not operating near focus, the near-field (Fresnel) propagation methods should
be used instead.
