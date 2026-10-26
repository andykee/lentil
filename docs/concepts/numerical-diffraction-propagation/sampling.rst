Sampling
========
The accuracy of a numerical diffraction simulation depends on adequately 
sampling the input and output planes. The sections below describe how to 
select appropriate sampling for these planes and how to configure the
propagations to avoid the introduction of numerical artifacts.

Ensuring Nyquist-sampled output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The relationship between spatial sampling in the input plane and output plane 
is defined by :math:`Q` and should be at least 2 in numerical simulations to 
ensure the output plane is Nyquist-sampled for intensity:

.. math::

    Q = \frac{\lambda \ F/\#}{du}

where :math:`lambda` is propagation wavelength, :math:`F/\#` is focal length 
divided by diameter, and :math:`du` is output plane spatial sampling. 
High-frequency ailiasing is clearly apparent in propagations where 
:math:`Q < 1.5` and still visibile to those with a keen eye when 
:math:`1.5 < Q < 2`:

.. plot:: user/plots/dft_discrete_Q_sweep.py
    :scale: 50

In cases where the optical system :math:`Q` is less than 2, the simulation 
fidelity should be increased by oversampling to avoid ailiasing. For a given 
imaging system, the system's :math:`F/\#`, output plane sampling, and 
propagation wavelength(s) will be fixed values. As a result, :math:`Q` can 
only be increased in a discrete propagation by introducing an ``oversample`` 
factor that effectively decreases the output plane sampling :math:`du`. In 
order to view the results of a propagation at native system sampling, the 
oversampled output plane must be resampled or rebinned.

.. math::

    Q_{\mbox{os}} = \frac{\lambda \ F/\# \ \texttt{oversample}}{du}

.. note::

    ``oversample`` should always be chosen to ensure :math:`Q > 2` for accurate
    propagation results.

Avoiding periodic wraparound
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

    \texttt{npix}_{\mbox{DFT}} = \frac{1}{2 \ \alpha \ \ \texttt{oversample}}

.. plot:: user/plots/dft_q_sweep.py
    :scale: 50