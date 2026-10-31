Controlling the output plane
============================
:func:`~lentil.propagate_dft` provides different ways for controlling the
exact region of the output plane that is computed. This can significantly
improve propagation performance in some circumstances.

.. image:: /_static/svg/propagate_prop_shape.svg
    :width: 450px
    :align: center

==================  ===========================================================
``shape``           Size of the output plane. The returned arrays are
                    ``shape * oversample``.
``prop_shape``      Size of the region that is computed, centered on the
                    optical axis. Defaults to ``shape``.
``shift``           Moves the optical axis, and the computed region with it,
                    within the output plane.
``mask``            Restricts computation to the bounding box of the nonzero
                    values in ``mask``.
==================  ===========================================================

Everything outside the computed region is zero. A |Wavefront| stores only
the computed data and the :attr:`~lentil.Wavefront.field` and
:attr:`~lentil.Wavefront.intensity` attributes place it in a ``shape`` array
when they are read.

``prop_shape``
--------------
``prop_shape`` sets the size of the computed region in the output. Decreasing
``prop_shape`` improves performance, but data outside the region is lost:

.. plot:: _plot/prop_shape.py
    :scale: 50

For most pupil to image propagations, a ``prop_shape`` of 128 or 256 pixels
balances speed and coverage. For image to pupil propagations, the computed
region must contain the whole re-imaged pupil, so it is safest to leave
``prop_shape`` at its default.

``prop_shape`` pays off most when modeling large tilts. See :doc:`tilt`.

``shift``
---------
``shift`` moves the optical axis by a number of output pixels. It is given in
native pixels and is scaled by ``oversample`` internally.

.. code-block:: python

    # move the optical axis +20 px in x and +10 px in y
    w = lentil.propagate_dft(w, pixelscale=5e-6, shape=128, shift=(20, 10))

Shifts from fitted tilt are added to ``shift``.

``mask``
--------
``mask`` restricts the computation to part of the output plane, such as a
detector subarray or within a field stop. It must have shape 
``shape * oversample``. Only the rectangular bounding box around its nonzero
values is computed, and the mask values are not applied to the result.

.. code:: python

    # only compute the region covered by an off-axis 40 x 60 px subarray
    mask = np.zeros((256, 256))
    mask[150:190, 160:220] = 1

    w = lentil.propagate_dft(w, pixelscale=5e-6, shape=256, mask=mask)

.. plot:: _plot/prop_mask.py
    :scale: 50
