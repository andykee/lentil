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

Creating Spectrum from CSV
--------------------------
Spectrum objects can be created from CSV files using 
:func:`radiometry.Spectrum.from_csv`. Given a CSV file formatted as

.. code-block::

    Wavelength (nm), Flux
    9.00452000e+01, 1.2757340e-17
    9.01354000e+01, 1.7265205e-17
    9.02258000e+01, 1.8341225e-17
    ... 
    2.98786752e+05, 1.3599320e-19
    2.99086286e+05, 1.3535845e-19
    2.99386120e+05, 1.3385094e-19

a Spectrum object is created with

.. code-block:: pycon

    >>> from lentil.radiometry import Spectrum
    >>> import matplotlib.pyplot as plt
    >>> vega = Spectrum.from_csv('vega.csv', waveunit='nm',
    ...                          valueunit='flam', header_rows=1)
    >>> plt.plot(wave=vega.wave, value=vega.value)
    >>> plt.grid()
    >>> plt.xlabel('Wavelength [nm]')
    >>> plt.ylabel('Flux [erg s^-1 cm^-2]')

.. image:: /_static/img/vega.png
    :width: 350px
    :align: center

Creating Spectrum from FITS
---------------------------
Although there is no built-in function, it is also possible to create a Spectrum
from a FITS file using `astropy <https://docs.astropy.org/en/stable/io/fits/index.html>`_.
Here, we'll create a dictionary of Johnson-Cousins filter transmissions as Spectrum
objects and then plot their transmissions:

.. code-block:: pycon
    
    >>> from lentil.radiometry import Spectrum
    >>> import matplotlib.pyplot as plt
    >>> from astropy.io import fits
    >>> jc = {}
    >>> for f in ('U','B','V','R','I'):
    ...     hdul = fits.open(f'johnson_{f.lower()}.fits')
    ...     jc[f] = (Spectrum(wave=hdul[1].data['WAVELENGTH'], 
    ...                       value=hdul[1].data['THROUGHPUT'],
    ...                       waveunit='nm', valueunit=None))
    >>> for band in jc:
    ...     plt.plot(jc[band].wave, jc[band].value, label=band)
    >>> plt.grid()
    >>> plt.legend()
    >>> plt.xlabel('Wavelength [nm]')
    >>> plt.ylabel('Transmission [A.U.]')

.. image:: /_static/img/johnson.png
    :width: 350px
    :align: center

The exact layout of spectral data within a FITS file may vary, but this example
illustrates a general approach for creating Spectrum objects from FITS data.

Units
-----
Spectrum objects provide support for a limited set of length and flux units and
allows for conversion between units.

The following wavelength units are supported:

===================== =========================
``waveunit``          Unit
===================== =========================
``m``, ``meter``      SI base unit
``um``, ``micron``    :math:`10^{-6}\ \mbox{m}`
``nm``, ``nanometer`` :math:`10^{-9}\ \mbox{m}`
``angstrom``          :math:`10^{-10}\ \mbox{m}`
===================== =========================

The following flux units are supported:

============= ===========================================
``valueunit`` Units
============= ===========================================
``photlam``   :math:`\mbox{photons s}^{-1} \mbox{m}^{-2}`
``wlam``      :math:`\mbox{W m}^{-2}`
``flam``      :math:`\mbox{erg s}^{-1} \mbox{cm}^{-2}`
============= ===========================================

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



