.. _user.getting_started.optical_system:

**************************
Describing optical systems
**************************

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

Planes in a simple optical system
=================================
Planes have attributes for representing 
Most optical systems can be adequately modeled by a single far-field propagation
between a :ref:`user.planes.pupil` and image plane. This includes most cameras,
telescopes, and imaging instruments. In these models, all of the optics in a system are
represented by a single |Pupil| plane:

.. plot::
    :scale: 50
    :include-source:

    >>> amplitude = lentil.circle(shape=(256, 256), radius=120)
    >>> opd = lentil.zernike_compose(mask=amplitude,
    ...                              coeffs=[0, 0, 0, 100e-9, 300e-9, 0, -100e-9])
    >>> pupil = lentil.Pupil(amplitude=amplitude, opd=opd, focal_length=10,
    ...                      pixelscale=1/240)
    >>> plt.imshow(pupil.opd)

More complicated optical systems
================================
More complicated imaging systems can still be moreled using far-field diffraction
but may contain multiple planes. This includes systems like spectrometers and 
coronagraphs which can typically still be represented by a series of |Pupil|, 
|Image|, and |Detector| planes. It is also possible to model and analyze effects 
that cannot be represented by pupil or image planes using |Plane| objects and the 
near-field propagation routines. In both multi-plane far-field and near-field
models, much more care needs to be taken to ensure each plane is adequately 
sampled to avoid the introduction of numerical artifacts in the diffraction 
propagation simulation.


