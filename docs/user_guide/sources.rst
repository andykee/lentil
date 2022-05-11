************************
Computational radiometry
************************

   .. currentmodule:: lentil

Computational radiometry is used to model the propagation of radiant energy through an
optical system. It uses geometry and known optical and imaging properties to compute the
irradiance from an observed scene at a detector. Lentil's
:ref:`radiometry<api.radiometry>` module provides a few helpful objects and functions
for working with ratiometric data. 

End to end radiometric modeling is part science and part art and this user guide 
provides limited insight. For a more in-depth treatment of the subject see:

* *Electro-Optical System Analysis and Design: A Radiometry Perspective*, Cornelius J. Willers
* *The Art of Radiometry*, James M. Palmer


The Spectrum Object
===================
:class:`radiometry.Spectrum` is a data structure for working with spectral 
data. A spectrum is defined by the following parameters:

* :attr:`~radiometry.Spectrum.wave` - Defines the wavelengths represented in the
  spectrum
* :attr:`~radiometry.Spectrum.value` - Defines the values corresponding to each
  wavelength
* :attr:`~radiometry.Spectrum.waveunit` - Defines the wavelength units
* :attr:`~radiometry.Spectrum.valueunit` - Defines the value units

Create a new Spectrum with

.. plot::
    :include-source:
    :scale: 50

    >>> import matplotlib.pyplot as plt
    >>> import lentil
    >>> s = lentil.radiometry.Spectrum(wave=np.arange(400,650),
    ...                                value=np.ones(250),
    ...                                waveunit='nm', valueunit=None)
    >>> plt.plot(s.wave, s.value)
    >>> plt.grid()
    >>> plt.xlabel(f'Wavelength [{s.waveunit}]')
    >>> plt.ylabel('Transmission [A.U.]')

Spectrum objects can be created from CSV files using 
:func:`radiometry.Spectrum.from_csv`


Although there is no built-in function, it is also possible to create a Spectrum
from a FITS file using `astropy <https://docs.astropy.org/en/stable/io/fits/index.html>`_



Manipulating Spectrum objects
-----------------------------



Units
-----

=========  ==============  =========  =======  ==============  ===================
Operand 1  Operation       Operand 2  Output   Units           Overlap
=========  ==============  =========  =======  ==============  ===================
Source     :math:`+`       Source     Source   Source 1 units  Full, partial, none
Source     :math:`=`       Source     Source   Source 1 units  Full, partial, none
Source     :math:`\times`  Scalar     Source   Source units    N/A
Source     :math:`\times`  Vector     Source   Source units    Full
=========  ==============  =========  =======  ==============  ===================

Radiometric Units
-----------------

Blackbody Emitters
==================

Thermal self emission for IR systems
------------------------------------

Optical Transmission
====================

The transmission of an optical element is a measure of the element's ability to transmit
light. Transmission models are primarily used in combination with a flux model to
estimate the irradiance of a source at a focal plane as observed through an optical
system.

In the simplest case, a transmission model can be represented by a single
:class:`~lentil.radiometry.Spectrum` defined analytically or loaded from a CSV file.
More complicated models can be built as the product of a number of independent
transmissions. A slight variation of this approach is used to represent transmission
models with selectable elements (like with a filter wheel).

Getting light into an optical system
====================================


Defining Sources
================
Two things to consider:

1. Need to compute the correct number of photons entering the system
    * Irradiance (<quantity>/m^2) * collecting area
    * collecting area is usually pi*(pm_diameter/2)^2 * fill_factor
2. Need to appropriately normalize the amplitude function and ensure the DFT is unitary



