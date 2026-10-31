Propagation with large tilts
============================
A tilt in the pupil plane appears as a translation in the image plane. When
tilt is represented by a linear phase ramp across a pupil, it is sampled like
any other phase. When this tilt is large, it causes two problems:

* **Aliasing.** If the phase ramp changes by more than half a wave per sample,
  it will alias to a smaller tilt and the resulting PSF will appear in the
  wrong place in the image plane.
* **Cost.** The output plane must be large enough to contain the shifted PSF.

.. image:: /_static/img/propagate_tilt_phase_wrap.png
    :width: 650px
    :align: center

Lentil avoids both issues by removing tilt from a plane's OPD and applying it
grometrically as a shift in the output plane.

.. image:: /_static/img/propagate_tilt_angle.png
    :width: 600px
    :align: center

How it works
------------
Calling a plane's :func:`~lentil.Plane.fit_tilt` method fits tip and tilt to
the plane's :attr:`~lentil.Plane.opd` via least squares. It removes the fitted
tilt from the OPD and appends an equivalent |Tilt| to the plane's
:attr:`~lentil.Plane.tilt` list.

When the plane is applied to a wavefront, its tilts are attached to the
wavefront's fields. During a propagation, :func:`~lentil.propagate_dft`
converts the attached tilts into a displacement in the output plane. The
whole-pixel part of the displacement moves the computed region
(``prop_shape``). The DFT applies the fractional part exactly, so the
result is not snapped to the output grid.

.. image:: /_static/img/propagate_tilt_angle_steps.png
    :width: 600px
    :align: center

Using ``fit_tilt()``
--------------------
.. code-block:: python

    pupil = lentil.Pupil(amplitude=amp, opd=opd, pixelscale=1/240,
                         focal_length=10)
    pupil = pupil.fit_tilt()

    w = lentil.Wavefront(650e-9) * pupil
    w = lentil.propagate_dft(w, pixelscale=5e-6, shape=600, prop_shape=64)

* :func:`~lentil.Plane.fit_tilt` returns a copy. Pass ``inplace=True`` to
  modify the plane in-place.
* The returned plane is frozen with the tilt-removed OPD cached. Calling
  :func:`~lentil.Plane.thaw` brings the fitted tilt back. See
  :doc:`../guides/performance`.
* On a segmented plane, tilt is computed independently for each segment. See
  :doc:`segmented`.

Other sources of tilt
---------------------
Tilt that is not part of a plane's OPD uses the same mechanism:

* A field angle or other image-plane shift: 
  ``lentil.Wavefront(wavelength, tilt=[rx, ry])``.
* An explicit tilt, such as a pointing offset: multiply by
  ``lentil.Tilt(x, y)``.
* Wavelength-dependent displacement, such as chromatic aberration or a prism:
  multiply by a :class:`~lentil.DispersiveTilt`.

Limitations
-----------
* :func:`~lentil.propagate_fft` does not currently support fitted tilt and
  raises ``NotImplementedError``.
* ``prop_shape`` must be large enough to hold the PSF.
* If a tilt moves the computed region entirely outside the output plane,
  nothing is computed for that field.