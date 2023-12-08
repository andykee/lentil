.. _user.fundamentals.radiometry:

.. currentmodule:: lentil

************************
Computational radiometry
************************

Computational radiometry is used to model the propagation of radiant energy through an
optical system. It uses geometry and known optical and imaging properties to compute the
irradiance from an observed scene at a detector. Lentil's
:ref:`radiometry<api.radiometry>` submodule provides a few helpful objects and functions
for working with ratiometric data. 

End to end radiometric modeling is part science and part art and this user guide 
provides limited insight. For a more in-depth treatment of the subject see:

* *Electro-Optical System Analysis and Design: A Radiometry Perspective*, Cornelius J. Willers
* *The Art of Radiometry*, James M. Palmer


The Spectrum Object
===================
:class:`radiometry.Spectrum` is a data structure for working with spectral 
data. A spectrum is defined by the following parameters:

* :attr:`~radiometry.Spectrum.wave` - the wavelengths represented in the
  spectrum
* :attr:`~radiometry.Spectrum.value` - the values corresponding to each
  wavelength
* :attr:`~radiometry.Spectrum.waveunit` - the wavelength units
* :attr:`~radiometry.Spectrum.valueunit` - the value units

Create a new Spectrum with

.. plot::
    :include-source:
    :scale: 50

    >>> import lentil
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
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
from a FITS file using `Astropy <https://docs.astropy.org/en/stable/io/fits/index.html>`_.
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
Spectrum objects provide support for a limited set of wavelength and flux units and
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

Converting between units is done with the Spectrum's :func:`~radiometry.Spectrum.to`
method:

.. code-block:: pycon

    >>> import lentil
    >>> import numpy as np
    >>> s = lentil.radiometry.Spectrum(wave=np.arange(400,700,50),
    ...                                value=np.ones(6),
    ...                                waveunit='nm', valueunit=None)
    >>> s.wave
    array([400, 450, 500, 550, 600, 650])
    >>> s.to('m')
    >>> s.wave
    array([4.0e-07, 4.5e-07, 5.0e-07, 5.5e-07, 6.0e-07, 6.5e-07])


Manipulating Spectrum objects
-----------------------------

Basic operations
~~~~~~~~~~~~~~~~
The following arithmetic operations are defined for Spectrum objects:

.. currentmodule:: lentil.radiometry

* :func:`Spectrum.add`
* :func:`Spectrum.subtract`
* :func:`Spectrum.multiply`
* :func:`Spectrum.divide`
* :func:`Spectrum.power`

A new Spectrum object is created and returned for each operation:

.. code-block:: pycon
   
    >>> import lentil
    >>> import numpy as np
    >>> # multiply a spectrum by a scalar
    >>> a = lentil.radiometry.Spectrum(wave=np.arange(400,700,50),
    ...                                value=np.ones(6),
    ...                                waveunit='nm', valueunit=None)
    >>> b = a * 2
    >>> b.value
    array([2, 2, 2, 2, 2, 2])

    >>> # add two spectrum together
    >>> c = lentil.radiometry.Spectrum(wave=np.arange(400,700,50),
    ...                                value=2*np.ones(6),
    ...                                waveunit='nm', valueunit=None)
    >>> d = a + c 
    >>> d.value
    array([3, 3, 3, 3, 3, 3])


Arithmetic operations work on Spectrum objects in the following ways:

* a scalar with a Spectrum elementwise over all Spectrum values
* a vector with a Spectrum elementwise over all Spectrum values (note the vector
  length must match the size of the Spectrum)
* a Spectrum with another Spectrum elementwise (note that all arithmetic 
  operations are supported for fully and partially overlapping data and addition 
  and subtraction are supported for disjoint data)

Standard arithmetic behavior is available using the appropriate overloaded operator 
(``+``, ``-``, ``*``, ``/``, or ``**``) with additional custom behavior defining
wavelength sampling and value interpolation options available by calling the 
arithmetic method directly:

.. code-block:: pycon
   
    >>> import lentil
    >>> import numpy as np
    >>> a = lentil.radiometry.Spectrum(wave=np.arange(400,700,50),
    ...                                value=np.ones(6),
    ...                                waveunit='nm', valueunit=None)
    >>> b = lentil.radiometry.Spectrum(wave=np.arange(500,700,25),
    ...                                value=2*np.ones(8),
    ...                                waveunit='nm', valueunit=None)
    >>> c = a.multiply(b, sampling=100)
    >>> c.value
    array([3, 3, 3])


Cropping, trimming, and joining operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following operations are available for manipulating Spectrum data:

* :func:`Spectrum.append` - Append a Spectrum to the end of another
* :func:`Spectrum.crop` - Crop a Spectrum by wavelength
* :func:`Spectrum.pad` - Pad a Spectrum with additional values
* :func:`Spectrum.trim` - Trim zeros off the ends of a Spectrum

Sampling and binning operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following operations are available for adjusting the sampling of
Spectrum data:

* :func:`Spectrum.bin` - Compute binned data at desired central wavelengths
* :func:`Spectrum.sample` - Sample a spectrum at desired wavelengths
* :func:`Spectrum.resample` - Sample a spectrum at desired wavelengths (in-place)
* :func:`Spectrum.integrate` - Compute integrated value between specified wavelengths

.. currentmodule:: lentil

Blackbody Emitters
==================
Create a :class:`radiometry.Blackbody` object with:

.. plot::
    :context:
    :include-source:
    :scale: 50

    >>> import lentil
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> wave = np.arange(400,4000)
    >>> temp = 5000
    >>> src = lentil.radiometry.Blackbody(wave,temp,waveunit='nm')
    >>> plt.plot(src.wave, src.value), plt.grid()
    >>> plt.xlabel('Wavelength [nm]'), plt.ylabel('Flux [photons/sec/m^2/sr]')

Because Blackbody subclasses Spectrum, all of the Spectrum methods are available:


.. plot::
    :context: close-figs
    :include-source:
    :scale: 50

    >>> src.to('wlam')
    >>> plt.plot(src.wave, src.value), plt.grid()
    >>> plt.xlabel('Wavelength [nm]'), plt.ylabel('Flux [W/m^2/sr]')


Path transmission and emission
==============================
The :class:`radiometry.Material` object is useful for representing an optic with a 
specific transmission (or reflectance), emission, and contamination level. Several
materials can be combined together in a list to compute path transmission or emission
using the :func:`radiometry.path_transmission` and :func:`radiometry.path_emission`
functions:



.. Thermal self emission for IR systems
.. ------------------------------------

.. Optical Transmission
.. ====================

.. The transmission of an optical element is a measure of the element's ability to transmit
.. light. Transmission models are primarily used in combination with a flux model to
.. estimate the irradiance of a source at a focal plane as observed through an optical
.. system.

.. In the simplest case, a transmission model can be represented by a single
.. :class:`~lentil.radiometry.Spectrum` defined analytically or loaded from a CSV file.
.. More complicated models can be built as the product of a number of independent
.. transmissions. A slight variation of this approach is used to represent transmission
.. models with selectable elements (like with a filter wheel).

.. Getting light into an optical system
.. ====================================


.. Defining Sources
.. ================
.. Two things to consider:

.. 1. Need to compute the correct number of photons entering the system
..     * Irradiance (<quantity>/m^2) * collecting area
..     * collecting area is usually pi*(pm_diameter/2)^2 * fill_factor
.. 2. Need to appropriately normalize the amplitude function and ensure the DFT is unitary



