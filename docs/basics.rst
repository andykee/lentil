******
Basics
******
Lentil provides a framework for defining optical planes and combining them to model
diffraction and other imaging effects. It can be used to create models with varying
levels of complexity; from simple prototypes to extremely high-fidelity models of ground
or space telescopes. Lentil was developed to be:

* *User friendly* - a simple, consistent interface optimized for common use cases.

* *Modular* - models are made by combining configurable building blocks. It is easy to
  replace or update individual elements without impacting the rest of the model.

* *Easy to extend* - core objects are designed to be subclassed, modified, extended,
  ripped apart, and put back together in different ways.


Importing Lentil
================
The recommended syntax for importing Lentil is

.. code-block:: python3

    >>> import lentil

If you'd prefer to use Lentil *en français*, you can try

.. code-block:: python3

    >>> import lentil as le

The functions and classes imported to the root namespace make up Lentil's core
functionality. Several additional submodules provide additional functionality. The
list below summarizes the general organization of Lentil's public namespace:

* :ref:`lentil<api-lentil>` - core classes and propagation functionality
* :ref:`lentil.convolvable<api-convolvable>` - classes for representing common imaging
  artifacts via convolution
* :ref:`lentil.detector<api-detector>` - functions to help with modeling focal planes
* :ref:`lentil.modeltools<api-modeltools>` - classes, decorators, and functions to
  simplify the creation of rich models
* :ref:`lentil.radiometry<api-radiometry>` - classes for performing computational
  radiometry
* :ref:`lentil.util<api-util>` - utility functions
* :ref:`lentil.wfe<api-wfe>` - classes for representing sources of wavefront error
* :ref:`lentil.zernike<api-zernike>` - classes for creating and working with Zernike
  polynomials

A Simple Example
================
