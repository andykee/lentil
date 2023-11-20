***********************************
Discrete Fourier transform sampling
***********************************

.. currentmodule:: lentil

Both the :func:`fourier.dft2` and :func:`fourier.czt2` methods require the user to
specify an output plane sampling in the frequency domain. Lentil calls this parameter
``alpha``. For a propagating wavefront, ``alpha`` is given by

.. math::

    \alpha = \frac{dx \ du}{\lambda \ z}

where

* :math:`dx` is the input plane sampling (spatial)
* :math:`du` is the output plane sampling (spatial)
* :math:`\lambda` is the propagation wavelength
* :math:`z` is the propagation distance

The accuracy of a numerical diffraction simulation depends on adequately sampling the
input and output planes. The sections below describe how to select appropriate sampling
for these planes to avoid the introduction of numerical artifacts.


Ensuring Nyquist-sampled output
===============================
The relationship between spatial sampling in the input plane and output plane is defined
by :math:`Q` and should be at least 2 in numerical simulations to ensure the output plane
is Nyquist-sampled for intensity:

.. math::

    Q = \frac{\lambda \ F/\#}{du}

High-frequency ailiasing is clearly apparent in propagations where :math:`Q < 1.5` and
still visibile to those with a keen eye when :math:`1.5 < Q < 2`:

.. plot:: _img/python/dft_discrete_Q_sweep.py
    :scale: 50


In cases where the system :math:`Q` is less than 2, the propagation fidelity should
be increased by oversampling to avoid ailiasing. For a given imaging system, the system's
:math:`F/\#`, focal plane sampling :math:`du`, and propagation wavelength(s)
:math:`\lambda` will be fixed values. As a result, :math:`Q` can only be increased in a
discrete propagation by introducing an ``oversample`` factor that effectively decreases
the output plane sampling :math:`du`. In order to view the results of a propagation at
native system sampling, the oversampled output plane must be resampled or rebinned.

.. math::

    \alpha_\mbox{os} = \frac{dx \ du}{\lambda \ z \ \texttt{oversample}}

and increases :math:`Q`:

.. math::

    Q_{\mbox{os}} = \frac{\lambda \ F/\# \ \texttt{oversample}}{du}

.. note::

    ``oversample`` should always be chosen to ensure :math:`Q > 2` for accurate
    propagation results.


Avoiding periodic wraparound
============================




.. math::

    q = \frac{1}{\alpha \ \ \texttt{npix} \ \ \texttt{oversample}}

.. math::

    \texttt{npix}_{\mbox{DFT}} = \frac{1}{2 \ \alpha \ \ \texttt{oversample}}

.. plot:: _img/python/dft_q_sweep.py
    :scale: 50
