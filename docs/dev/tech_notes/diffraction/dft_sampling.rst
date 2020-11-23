***********************************
Discrete Fourier transform sampling
***********************************

.. currentmodule:: lentil

Both the :func:`fourier.dft2` and :func:`fourier.czt2` methods require the user to 
provide an output plane sampling (frequency) parameter ``alpha``:

.. math::

    \alpha = \frac{dx \ du}{\lambda \ z}

where

* :math:`dx` is the input plane sampling (spatial)
* :math:`du` is the output plane sampling (spatial)
* :math:`\lambda` is the propagation wavelength
* :math:`z` is the propagation distance

Ensuring Nyquist-sampled output
===============================
The relationship between input plane sampling and output plane sampling is defined
by :math:`Q` and should be at least 2 to ensure the output plane is Nyquist-sampled
for intensity:

.. math::

    Q = \frac{\lambda \ F/\#}{du}

In cases where the system :math:`Q` is less than 2, the propagation fidelity should 
be increased by oversampling to avoid ailiasing. High-frequency ailiasing is clearly 
apparent in propagations where :math:`Q < 2`:

.. image:: /_static/img/discrete_Q_sweep.png
    :width: 800px
    :align: center

For a given imaging system, :math:`du`, :math:`\lambda`, and :math:`z` are likely to 
be fixed. Therefore, :math:`Q` can only be increased in a discrete propagation by either 
more finely sampling the input plane (decrease :math:`dx`) or by oversampling the output
plane (temporarily decrease :math:`du`) before rebinning the results to the native 
output sampling. 

We more finely sample the output plane by introducing an ``oversample`` factor. The
result is a smaller ``alpha`` term


.. math::

    \alpha_\mbox{os} = \frac{dx \ du}{\lambda \ z \ \texttt{oversample}}

and a larger :math:`Q`:

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