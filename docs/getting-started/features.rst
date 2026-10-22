Features
========

Most of the functionality in Lentil resides in its standard library. This
includes:

* Objects for representing :doc:`elements of an optical system <../concepts/planes/index>`
* Tools for modeling various sources of :doc:`wavefront error <../concepts/planes/wavefront-error>`
  and common :doc:`imaging artifacts <../concepts/imaging-artifacts>`
* A suite of algorithms for :doc:`numerically simulating wave optics propagation <../concepts/numerical-diffraction-propagation/index>`
* Various utility functions for drawing shapes, manipulating array-based data,
  and performing routine optical calculations

Lentil also includes several subpackages that provide domain-specific
functionality. These include:

* ``lentil.detector`` provides tools for modeling :doc:`focal plane arrays <../concepts/image-sensors>`
* ``lentil.radiomatry`` provides tools for working with :doc:`spectral data <../concepts/radiometry>`