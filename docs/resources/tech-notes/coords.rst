A consistent coordinate system for wave optics simulations
==========================================================

A choice of harmonic time dependence fixes the sign of all phase factors throughout the diffraction formalism, including wave propagation operators, Fourier transforms, Optical Path Difference (OPD) definitions, and numerical propagation algorithms. Consistent use of a single convention is essential to avoid sign errors in simulated wavefronts and PSFs and when incorporating laboratory measurements or results from other analysis. Key results having implications for wave optics simulation are summarized in Table 1. The remainder of this document describes how the conventions in the table were derived through the specification of a sign convention for harmonic time dependence. 

**Table 1**: Summary of harmonic time dependence sign convention and implications for wave optics simulation and phase retrieval.

+-------------------------------+-----------------------------------------------------+----------------------------------------------------+
| Time convention               | :math:`\exp(-i\omega t)`                            | :math:`\exp(i\omega t)`                            |
+===============================+=====================================================+====================================================+
| Spatial phasor contribution   | :math:`\displaystyle\frac{2\pi i}{\lambda}W`        | :math:`\displaystyle\frac{-2\pi i}{\lambda}W`      |
+-------------------------------+-----------------------------------------------------+----------------------------------------------------+
| Forward propagation           | ``fft()``                                           | ``ifft()``                                         |
+-------------------------------+-----------------------------------------------------+----------------------------------------------------+
| Positive OPD                  | Longer path length                                  | Shorter path length                                |
+-------------------------------+-----------------------------------------------------+----------------------------------------------------+
| Positive wavefront error      | Longer path length                                  | Shorter path length                                |
+-------------------------------+-----------------------------------------------------+----------------------------------------------------+
| Positive focus error          | Longer focal length                                 | Shorter focal length                               |
+-------------------------------+-----------------------------------------------------+----------------------------------------------------+
| Positive tilt about :math:`x` | :math:`-v` translation                              | :math:`+v` translation                             |
+-------------------------------+-----------------------------------------------------+----------------------------------------------------+
| Positive tilt about :math:`y` | :math:`+u` translation                              | :math:`-u` translation                             |
+-------------------------------+-----------------------------------------------------+----------------------------------------------------+
| Converging spherical wave     | :math:`\exp\left[-i\frac{k}{2z} (x^2 + y^2)\right]` | :math:`\exp\left[i\frac{k}{2z} (x^2 + y^2)\right]` |
+-------------------------------+-----------------------------------------------------+----------------------------------------------------+


Coordinate System
-----------------

We begin with the right-handed coordinate frame depicted in Figure 1. A diffracting aperture lies in the :math:`(x,y)` plane with the optical axis :math:`z` pointing in the direction of propagation. The diffracted wave field is observed in the :math:`(u,v)` plane which is parallel to the aperture plane and separated by some distance :math:`+\Delta z`.

.. figure:: /_static/img/coords_3d.png
    :width: 500px
    :align: center

    **Figure 1**: Coordinate system.


Wave Propagation Sign Convention
--------------------------------
The phase of a wave traveling in :math:`+z` takes the form :math:`\cos\left(kz-\omega t\right)` or equivalently :math:`\cos\left(\omega t-kz\right)`. The choice of sign for the time component (harmonic time dependence) has downstream effects for wave optics simulations, including determining the sign convention for representing the complex pupil amplitude, the direction of Fourier transform used to model diffraction propagation, and whether positive OPD corresponds with longer or shorter path length.

Under the convention used most commonly in physics and optics, the harmonic time dependence term is negative. In complex notation, this gives us

.. math::
    :label: 

    E(z,t) = E_0\exp\left(ikz\right) \exp\left(-i\omega t\right) \, ,

corresponding with clockwise rotation in the complex plane. Extending Equation 1 to two dimensions gives

.. math::
    :label: 

    u(x,y,t) = f(x,y)\exp(-i\omega t) \, ,

where

.. math::
    :label: 

    f(x,y) = A(x,y)\exp\left[\frac{2\pi i}{\lambda}W(x,y)\right]

is the complex pupil function with amplitude :math:`A`, wavefront aberration :math:`W`, and wavelength :math:`\lambda`. Propagating this scalar field to an observation plane is described using the scalar diffraction integral derived from the wave equation. When the propagation distance is much larger than the aperture size and wavelength, the far-field (Fraunhofer) approximation may be applied. In the far-field, the quadratic phase terms associated with the propagation become negligible relative to the linear phase variation in the field and may therefore be omitted. The resulting field in the observation plane is given by

.. math::
    :label:

    F(u,v) = \iint_{-\infty}^\infty f(x,y) \exp\left[\frac{-2\pi i}{\lambda z} \left(xu + yv\right)\right] \,dx \,dy \, .

From the expression for the two-dimensional continuous Fourier transform

.. math::
    :label:

    \displaystyle F(u,v) = \iint_{-\infty}^{\infty}f(x,y)\exp\left[-2\pi i\left(ux + vy\right)\right]\,dx \,dy

we see that Equation 4 is simply the Fourier transform of the complex pupil function :math:`f(x,y)` evaluated at spatial frequencies

.. math::
    :label:

    f_x = u/(\lambda z)\, , \hspace{3em} f_y = v/(\lambda z) \, .

Most practical applications of Equation \ref{eq:fraunhofer} use the Fast Fourier Transform (FFT) implementation of the Discrete Fourier Transform (DFT) to propagate from pupil plane to image plane. Across common scientific computing environments (NumPy, SciPy, MATLAB, IDL, FFTW, FFTPACK) the Fourier transform sign convention is consistent: ``fft()`` has a negative complex exponential and ``ifft()`` has a positive complex exponential. If a custom DFT implementation is used, care should be taken to ensure the sign convention used matches the standard Fourier transform sign convention.

Wavefront Error Sign Convention
-------------------------------
The sign convention used in the representation of wavefront error follows directly from the convention chosen for harmonic time dependence in Section \ref{3}. Taking :math:`\exp\left(-i\omega t\right)` as the time-dependent component, we rewrite Eqs. (2) and (3) to represent a wave propagating in the :math:`+z` direction as

.. math::
    :label:

    u(x,y,t) = A(x,y)\exp\left(\frac{2\pi i}{\lambda}z\right) \exp\left(-i\omega t\right)

with amplitude :math:`A`, wavelength :math:`\lambda`, and position along the optical axis :math:`z`. From the spatial component of Equation \ref{eq:full_phasor} we see that the spatial phase accumulating after traveling some distance :math:`z` is 

.. math::
    :label:

    \phi = \frac{2\pi}{\lambda}z

and that a wave traveling farther in :math:`z` will have a more positive phase :math:`\phi`.

We define OPD as the aberrated optical path length minus the reference optical path length. [#f1]_ In two-dimensions, the OPD evaluated at each pupil coordinate is:

.. math::
    :label:

    W(x,y) = z(x,y) - z_{ref}(x,y) \, .

Finally, noting that wavefront error is simply the phase corresponding to a given OPD:

.. math::
    :label:
    
    \phi(x,y) = \frac{2\pi}{\lambda}W(x,y) \, ,

it follows that longer path length or OPD corresponds with positive wavefront error and with a more positive phase.

We now turn our attention to examining the impact of this sign convention on two of the most commonly observed optical aberrations: tilt and focus. Tilt in the pupil plane is represented by a linear phase ramp

.. math::
    :label:

    \phi(x,y) = \frac{2\pi}{\lambda}\left(\theta_x x + \theta_y y\right)

From Equation \ref{eq:pupil}, the corresponding complex pupil function is

.. math::
    :label:

    f(x,y) = A(x,y)\exp\left[\frac{2\pi i}{\lambda} \left(\theta_x x + \theta_y y\right)\right] \, .

Substituting into Equation \ref{eq:fraunhofer} gives

.. math::
    :label:

    F(u,v) = \iint_{-\infty}^\infty A(x,y) \exp\left[\frac{-2\pi i}{\lambda z} \left(x\left(u-\theta_x z\right) + y\left(v-\theta_y z\right)\right)\right] \,dx \,dy \, .

By the Fourier shift theorem, the sign of the observation plane shift follows directly from the sign of the angular deviation of the propagation vector due to the tilt:

.. math::
    :label:

    u = \theta_x z\, , \hspace{3em} v=\theta_y z\, .

Using the right hand rule, a positive rotation about :math:`x` rotates the direction of propagation in :math:`-y`, resulting in a shift in :math:`-v` in the observation plane. For a rotation about :math:`y`, the propagation direction rotates in :math:`+x` and a corresponding shift in :math:`+u` in the observation plane. Note that this asymmetry in sign arises from the geometry of the coordinate system. The Fourier shift theorem itself treats :math:`x` and :math:`y` symmetrically, and the sign of the image plane shift follows directly from the sign of the OPD slope, which is set by the right-hand rule.

For focus, we return to the convention for wavefront error and OPD identified earlier: positive OPD corresponds with a longer path length. From Figure 2 we see that a wavefront with positive focus error has a focal point displaced by :math:`+\Delta z` relative to the nominal focal point (longer focal length). Relative to the reference wavefront, rays near the edge of the pupil travel a longer optical path before reaching focus, producing a positive OPD.

.. figure:: /_static/img/focus_direction.png
    :width: 500px
    :align: center

    **Figure 2**: Positive focus corresponds with a wavefront whose focal length is longer relative to a reference sphere. A wavefront with negative focus has a shorter focal length.

Validation via numerical simulation
-----------------------------------

.. image:: /_static/img/array_indexing.png
        :width: 500px
        :align: center



Tip/tilt
--------
A positive x-tilt rotates the yz plane clockwise about the x-axis resulting
in a shift in the image plane in the positive y direction. A positive y-tilt
rotates the xz plane clockwise about the y-axis resulting in a shift in the
image plane in the negative x direction.

.. plot:: _plot/tilt_images.py
    :scale: 50

Focus
-----
A positive focus has the effect of shifting the focus point beyond the image plane (+z)
while a negative focus has the effect of shifting the focus point in front of the image
plane (-z).

.. image:: /_static/img/focus_direction.png
    :width: 450px
    :align: center

We can test this by observing +/- defocused point spread functions of an
imaging system with an asymmetric aperture (for example, in the shape of the letter
F). We expect the positive focus image to have the same orientation as the aperture
(consistent with observing the image before coming to focus) and the negative focus
image to be flipped about both axes relative to the aperture (consistent with
observing the image after passing through focus). The results of this exercise are
presented below:

.. plot:: _plot/focus_images.py
    :scale: 50



.. rubric:: Footnotes

.. [#f1] Modern raytrace packages like Code V and Zemax use the opposite sign convention.



