***********
Signal flow
***********

An image sensor converts photons striking a pixel into a digital signal. The signal
produced by each pixel is subject to a number of physical characteristics and is
degraded by a variety of noise sources. The core components of this process are
represented in the figure below:

.. image:: /_static/img/image_sensor_signal_flow.png
    :width: 75%


Photons striking a pixel during an integration time will result in the accumulation of
electrons in the pixel well. The fraction of photons striking a pixel that are collected
as electrons is the sensor's *quantum efficiency* (QE). Additional electrons "leaking"
into each pixel are represented as *dark signal*. The electrons collected in each pixel
during an integration time are converted to a voltage (up to some saturation point),
amplified, and converted to a digital signal by an ADC. Additional noise generated
during the readout and conversion of the analog signal is represented as *read noise*.

This is a somewhat simplistic model of a digital image sensor but it provides the
framework for capturing the majority of physical effects and noise sources present and
is easily extended to capture more complex and nuanced effects.

Input signal
============

Charge collection
=================


Pixel effects
=============


Analog to digital conversion
============================

