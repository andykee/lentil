Multi-plane propagation
=======================
Many systems have intermediate focal planes and re-imaged pupils, including
relays, coronagraphs, spatial filters, and instruments behind a telescope that
provides a collimated interface. Lentil models these by chaining propagations
between planes. Each call to :func:`~lentil.propagate_fft` or 
:func:`~lentil.propagate_dft` moves the wavefront between a pupil and an
image plane, and planes multiplied in between act at those conjugates.

Focal length and diameter through a chain
-----------------------------------------
A wavefront carries a ``focal_length`` and a ``diameter`` from plane to
plane:

* Multiplying by a |Pupil| with a ``focal_length`` sets both. A |Pupil|
  without one, such as a Lyot stop, leaves them unchanged.
* :func:`~lentil.propagate_dft` uses the wavefront's ``focal_length`` unless a
  ``focal_length`` argument is given. The returned wavefront carries whichever
  value was used for the propagation.
* An image to pupil propagation with ``focal_length=f2`` scales the
  wavefront's ``diameter`` by :math:`f_2/f_1`.

The ``pixelscale`` of the re-imaged pupil is up to you. To keep the same
number of samples across the pupil, scale it by :math:`f_2/f_1` as well.

Example: a relay with a field stop
----------------------------------
A 1 m, f/10 telescope focuses onto a circular field stop with a 26 µm
(40 :math:`\lambda /D`) radius. A 5 m focal length collimator re-images the
pupil with a 0.5 m diameter:

.. code-block:: python

    pupil = lentil.Pupil(amplitude=amp, pixelscale=1/240, focal_length=10)
    stop = lentil.Image(amplitude=lentil.circle((256, 256), 26),
                        pixelscale=1e-6)

    w = lentil.Wavefront(650e-9) * pupil
    w = lentil.propagate_dft(w, to=stop)    # focus on the stop's grid
    w = w * stop

    w = lentil.propagate_dft(w, pixelscale=1/480, shape=256,
                             focal_length=5)

.. code-block:: pycon

    >>> w.ptype, w.focal_length, w.diameter
    (ptype('pupil'), 5, np.float64(0.5))

.. plot:: _plot/prop_relay.py
    :scale: 50

Without the stop, the pupil is re-imaged with good quality and crisp edges.
Including the stop masks high spatial frequencies in the image plane,
resulting in a blurred re-imaged pupil. Note that only a small portion of the
focal plane was computed, which also acts a s an additional square field stop.
Propagating to a larger image plane would resuce the ringing seen in the
pupil image with no stop. See :doc:`sampling` for details.

Pupil plane optics after the relay are |Pupil| planes. To form a final
image, propagate again:

.. code-block:: python

    lyot = lentil.Pupil(amplitude=lentil.circle((256, 256), 108),
                        pixelscale=1/480)
    w = w * lyot
    w = lentil.propagate_dft(w, pixelscale=2e-6, shape=128)

With no ``focal_length`` argument, the propagation uses the collimator's
5 m focal length as if the light were refocused by the same optic. If the
camera has a different focal length, pass it as ``focal_length`` or give it to
a |Pupil| plane to be multiplied with the propagating wavefront.

.. note::

    A stop placed exactly at focus is modeled here with far-field
    propagations only. A stop away from focus needs a near-field
    propagation. See :doc:`near-field`.