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

.. image:: /_static/img/discrete_Q_sweep.png
    :width: 800px
    :align: center

In cases where the system :math:`Q` is less than 2, the propagation fidelity should
be increased by oversampling to avoid ailiasing. For a given imaging system, :math:`du`, 
:math:`\lambda`, and :math:`z` are likely to be fixed. Therefore, :math:`Q` can only be 
increased in a discrete propagation by either more finely sampling the input plane 
(decreasing :math:`dx`) or by oversampling the output plane (temporarily decreasing 
:math:`du`) before rebinning the results to the native output sampling.

We more finely sample the output plane by introducing an ``oversample`` factor. This
decreases ``alpha``


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


.. image:: /_static/img/q_sweep.png
    :width: 800px
    :align: center
