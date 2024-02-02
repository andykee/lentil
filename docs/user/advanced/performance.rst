.. _user.performance:

**********************
Optimizing Performance
**********************

Selecting sensible propagation parameters
=========================================
Performing numerical propagation is almost always the most expensive part of 
any model. Even the most efficient and streamlined models will experience
high memory usage and slow performance when working with large arrays and many
wavelengths. Understanding how the propagation method's parameters influence 
accuracy and speed will allow you to identify appropriate values for your 
specific needs.

``propagate_dft()``
-------------------

.. list-table:: 
   :widths: 20 40 40
   :header-rows: 1

   * - Parameter
     - Usage
     - Performance considerations
   * - ``wave``, ``weight``
     - Arrays representing wavelength and the corresponding weight of each 
       wavelength in the propagation.
     - Propagation time scales linearly with the number of wavelengths. 
   * - ``shape``
     - Shape of output plane.
     - The output plane size has minimal impact on propagation performance 
       unless it is astronomically large. ``shape`` should be set to ensure 
       all data is adequately captured by the output plane.
   * - ``prop_shape``
     - Shape of propagation plane.
     - Propagation time scales quadtatically with ``prop_shape``. 
       ``prop_shape`` should be chosen to capture the data extent necessary to
       accurately model all required diffraction effects.
   * - ``oversample``
     - Number of times to oversample the output and propagation planes.
     - Propagation time scales quadratically with ``oversample``. For accuracy, 
       ``oversample`` should be selected to ensure propagations are Nyquist 
       sampled, but there is typically no benefit in selecting larger values.



Strategies for increasing runtime performance
=============================================

Perform fewer propagations by increasing wavelength sampling
------------------------------------------------------------
Selecting an appropriate wavelength sampling depends on many factors. Any 
opportunity to reduce the number of propagations required will increase 
overall model performance. Two patterns for more coarse wavelength sampling
are provided :ref:`here<examples.bandpass_resampling>`


Use appropriate plane sampling
------------------------------
Planes should be sized to ensure the smallest spatial features of interest are
adequately sampled. Small features (e.g. secondary mirror supports, segment 
gaps) should be represented by at leat two samples. 



.. Faster photon to electron (quantum efficiency) calculations
.. -----------------------------------------------------------
