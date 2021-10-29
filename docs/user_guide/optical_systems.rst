.. _user_guide.optical_systems:

**************************
Describing optical systems
**************************

.. |Pupil| replace:: :class:`~lentil.Pupil`
.. |Image| replace:: :class:`~lentil.Image`
.. |Detector| replace:: :class:`~lentil.Detector`
.. |Plane| replace:: :class:`~lentil.Plane`
.. |Wavefront| replace:: :class:`~lentil.Wavefront`
.. |Tilt| replace:: :class:`~lentil.Tilt`

Lentil uses |Plane| objects to represent discretely sampled planes in an
optical system and Fourier transform-based algorithms to numerically model
the propagation of a monochromatic |Wavefront| from plane to plane through
the optical system. The basic propagation flow is:

    1. Create a new |Wavefront|
    2. Propagate the wavefront through the first plane in the optical system
    3. Propagate the wavefront to the next plane in the optical system
    4. Repeat steps 2 and 3 until reaching the final plane in the optical system

This process may be repeated a number of times (once per discrete wavelength being
represented) to simulate the propagation of broadband light through the system.

.. _user_guide.optical_systems.plane_wavefront:

How a plane affects a wavefront
===============================
An optical plane generally has some effect on a wavefront as it propagates
through the plane. A plane may change a propagating wavefront's amplitude, phase,
or physical extent. This |Plane|-|Wavefront| interaction is governed by the
plane's :func:`~lentil.Plane.multiply` method. This method constructs a complex
phasor from the plane's :attr:`~lentil.Plane.amplitude` and
:attr:`~lentil.Plane.phase` attributes and the |Wavefront| wavelength. The plane
complex phasor is then multiplied element-wise with the wavefront's complex data
array:

.. math::

    \mathbf{W} = \mathbf{A} \exp(\frac{-2\pi j}{\lambda} \mathbf{\theta}) \circ \mathbf{W}

If the |Plane| :attr:`~lentil.Plane.tilt` attribute is not empty, its contents are appended
to the |Wavefront|. See :ref:`user_guide.planes.fit_tilt` and :ref:`user_guide.diffraction.tilt`
for additional details.

Planes in a simple optical system
=================================
Most optical systems can be adequately modeled by far-field propagation between a pupil
and image plane. This includes most cameras, telescopes, and imaging instruments. In
these models, all of the optics in a system are represented by a single
:class:`~lentil.Pupil` plane. The results of the diffraction propagation (assuming the
optical system is operating near focus) can be viewed using either an
:class:`~lentil.Image` or :class:`~lentil.Detector` plane. Because of the optimizations
mentioned above, the :class:`~lentil.Detector` plane should be used to most efficiently
compute image plane intensity.

.. image:: /_static/img/cassegrain.png
    :width: 550px
    :align: center

.. image:: /_static/img/simple_optical_system.png
    :width: 375px
    :align: center

If the optical system depicted above has a 1m diameter primary mirror, a secondary
mirror obscuration of 0.33m centered over the primary, a focal length of 10m, and
a focal plane with 5um pixels. We describe this optical system using a |Pupil| and
|Detector| plane as follows:

.. code-block:: pycon

    >>> amplitude = lentil.circle(shape=(256, 256), radius=128) -
    ...             lentil.circle(shape=(256, 256), radius=128/3)
    >>> pupil = lentil.Pupil(amplitude=amplitude, phase=opd, focal_length=10,
    ...                      pixelscale=1/256)
    >>> detector = lentil.Detector(pixelscale=5e-6)

Segmented optical systems
=========================
Creating a model of a segmented aperture optical system in Lentil doesn't require any
special treatment. The |Plane| and |Pupil| objects work just as well with sparse or
segmented amplitude, phase, and mask attributes as with monolithic ones.

That being said, it is advantageous from a performance point of view to supply a
3-dimensional "segment mask" when specifying a Plane's :attr:`~lentil.Plane.mask`
attribute rather than a flattened 2-dimensional mask when working
with a segmented aperture, as depicted below:

.. image:: /_static/img/segmask.png
    :width: 650px
    :align: center

This modification is not necessary to achieve physically correct propagations, but can
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
