**********
Radiometry
**********
Computational radiometry is used to model the propagation of radiant energy through an
optical system. It uses geometry and known optical and imaging properties to compute the
irradiance from an observed scene at a detector. Lentil's
:ref:`radiometry<api-radiometry>` module provides a
:class:`~lentil.radiometry.Spectrum` object for representing and working with
radiometric quantities but because radiometric models tend to be very application=
specific, not much additional help is provided.

The Spectrum Object
===================
:class:`~lentil.radiometry.Spectrum` is a data structure for holding spectral data.

=========  ==============  =========  =======  ==============  ===================
Operand 1  Operation       Operand 2  Output   Units           Overlap
=========  ==============  =========  =======  ==============  ===================
Source     :math:`+`       Source     Source   Source 1 units  Full, partial, none
Source     :math:`=`       Source     Source   Source 1 units  Full, partial, none
Source     :math:`\times`  Scalar     Source   Source units    N/A
Source     :math:`\times`  Vector     Source   Source units    Full
=========  ==============  =========  =======  ==============  ===================

Blackbody Emitters
==================

Radiometric Units
=================

Defining Sources
================
Two things to consider:

1. Need to compute the correct number of photons entering the system
    * Irradiance (<quantity>/m^2) * collecting area
    * collecting area is usually pi*(pm_diameter/2)^2 * fill_factor
2. Need to appropriately normalize the amplitude function and ensure the DFT is unitary

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

Thermal Emission for IR Systems
===============================

