.. _user.overview:

****************
Package overview
****************

Lentil is a Python library for modeling the imaging chain of an optical system. It was 
originally developed at NASA's Jet Propulsion Lab by the Wavefront Sensing and Control 
group (383E) to provide an easy to use framework for simulating point spread functions 
of segmented aperture telescopes.

Lentil provides a framework for defining optical planes and combining them to model 
diffraction and other imaging effects. It can be used to create models with varying 
levels of complexity; from simple prototypes to extremely high-fidelity models of ground 
or space telescopes. Lentil was developed to be:

* User friendly - a simple, consistent interface optimized for common use cases.
* Modular - models are made by combining configurable building blocks. It is easy to 
  replace or update individual elements without impacting the rest of the model.
* Easy to extend - core objects are designed to be subclassed, modified, extended, 
  broken apart, and put back together in different ways.

Lentil organization
===================
Lentil is organized into two parts: a standard library and a few subpackages 
providing domain-specific tools. The features of the standard library and each
of the subpackages is described in the table below:

============================================== ===============================================
Namespace                                      Purpose
============================================== ===============================================
:ref:`lentil <api>`                            Standard library
:ref:`lentil.detector <api.detector>`          Model focal planes
:ref:`lentil.radiometry <api.radiometry>`      Work with spectral data
============================================== ===============================================

Lentil's standard library can be imported as follows:

.. code-block:: python3

    >>> import lentil

If you'd prefer to use Lentil *en français*, you can try

.. code-block:: python3

    >>> import lentil as le

Public subpackages are automatically imported with Lentil.


Getting help
============
The best place to ask for help on subjects not covered in this documentation or suggest new 
features/ideas is by opening a ticket on `Github <https://github.com/andykee/lentil/issues>`__.
