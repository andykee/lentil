Sampling
========
The accuracy of a numerical diffraction simulation depends on adequately 
sampling the input and output planes. This section describes how to select
appropriate sampling fand how to configure propagations to avoid introducing
numerical artifacts.

Output sampling and aliasing
----------------------------
The relationship between spatial sampling in the input plane and output plane 
is described by 

.. math::

    Q = \frac{\lambda \ F/\#}{\Delta u}

where :math:`\lambda` is wavelength, :math:`F/\#=f/D` is focal length 
divided by the pupil diameter, and :math:`\Delta u` is output pixelscale. 

The intensity of a diffraction pattern contains no spatial frequencies above
:math:`1/(\lambda F/\#)`, so the output is Nyquist sampled for intensity when
:math:`Q \ge 2`. Aliasing is obvious below :math:`Q = 1.5`:

.. plot:: _plot/dft_discrete_Q_sweep.py
    :scale: 50

A helper function :func:`~lentil.pixelscale_nyquist` returns pixelscale that
gives :math:`Q = 2`.

Oversampling
~~~~~~~~~~~~
If detector pixel size, wavlength, and F/# are fixed and :math:`Q < 2`, the
output can be computed on a finer grid by specifying ``oversample``. The
result has ``shape * oversample`` pixels sampled at 
``pixelscale / oversample``, and

.. math::

    Q_{\mathrm{os}} = Q \cdot \mathtt{oversample}

Choose ``oversample`` so that :math:`Q_{\mathrm{os}} \ge 2`, then rebin the
intensity back to native sampling:

.. code-block:: python

    # f/10, 650 nm, 5 um pixels -> Q = 1.3
    oversample = int(np.ceil(2 / 1.3))    # 2

    w = lentil.propagate_dft(w, pixelscale=5e-6, shape=64,
                             oversample=oversample)
    img = lentil.rebin(w.intensity, oversample)    # (64, 64)

Input sampling and periodic wraparound
--------------------------------------
The Fourier transform of a sampled pupil is periodic. The output repeats with
period

.. math::

    P = \frac{\lambda z}{\Delta x} = N \lambda F/\#

where :math:`\Delta x` is the input ``pixelscale``, :math:`z` is the focal
length, and :math:`N = D/\Delta x` is the number of samples across the pupil.
Light more than :math:`P/2` from the optical axis wraps around to the other
side of the output plane:

.. plot:: _plot/sampling_period.py
    :scale: 50

Be careful to limit the output field of view to fall within half a period so
that neighboring copies don't contaminate the output:

.. math::

    \mathtt{shape} \cdot \Delta u \le \frac{P}{2}
    \quad\Longleftrightarrow\quad
    \mathtt{shape} \le \frac{N \, Q}{2}

where ``shape`` is in native (not oversampled) pixels. For example, a pupil
240 pixels across imaged at :math:`Q = 1.3` supports an output up to 156
pixels across. To go larger, sample the pupil more finely with
:func:`~lentil.Plane.rescale`.

.. plot:: _plot/dft_q_sweep.py
    :scale: 50

The same limit applies to tilt: a tilt that moves the PSF more than
:math:`P/2` aliases. Lentil avoids this limit with a separate path for
tilt, described in :doc:`tilt`.

Image to pupil propagations
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The same relations hold for image to pupil propagation, although the effects
are slightly different:

* The re-imaged pupil must fit inside one period, :math:`\lambda z / \Delta x`,
  where :math:`\Delta x` is now the image plane ``pixelscale``.
* The image plane field of view, ``shape * pixelscale``, limits the finest
  detail in the re-imaged pupil to roughly :math:`\lambda z / (\mathtt{shape}
  \cdot \Delta x)`. A small image plane acts as a field stop. See
  :doc:`multi`.

Choosing sampling parameters
----------------------------
#. Set ``pixelscale`` to the detector pixel size, or to
   :func:`~lentil.pixelscale_nyquist` if no detector is modeled.
#. Compute :math:`Q` and set ``oversample`` to at least :math:`2/Q`.
#. Set ``shape`` to cover the region of interest. Check that
   :math:`\mathtt{shape} \le N Q / 2`. If it isn't, sample the pupil more
   finely or reduce ``shape``.
#. For broadband models, run these checks at the **shortest** wavelength,
   where both :math:`Q` and :math:`P` will be smallest.